#include "BDSXtrackInterface.hh"
#include <cstring>
#include <cmath>
#include <BDSSamplerCustom.hh>
#include <iostream>
#include <vector>
#include <string>

#include "TSystem.h"
#include "TError.h"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include <fstream>
#include <unistd.h>
#include <fcntl.h>


// ------------------------------------------ //
// Logging functionality for Geant4 and BDSIM //
// ------------------------------------------ //

class XcollLogFileSession : public G4UIsession {
public:
  XcollLogFileSession()
  : logfile("geant4.out"), errfile("geant4.err") {}

  virtual ~XcollLogFileSession() { logfile.close(); errfile.close();}

  G4int ReceiveG4cout(const G4String& msg) override {
    logfile << msg;
    return 0;
  }
  G4int ReceiveG4cerr(const G4String& msg) override {
    errfile << msg;
    return 0;
  }

private:
  std::ofstream logfile;
  std::ofstream errfile;
};

void RedirectGeant4() {
  auto ui = G4UImanager::GetUIpointer();
  static XcollLogFileSession session;
  ui->SetCoutDestination(&session);
}

static std::ofstream rootInfoLog("root.out");
static std::ofstream rootErrLog("root.err");

static void MyRootErrorHandler(int level, Bool_t abort, const char *location, const char *msg)
{
    if (level >= kError) {
        rootErrLog << location << ": " << msg << '\n';
    } else {
        rootInfoLog << location << ": " << msg << '\n';
    }
    // mimic ROOT's abort behaviour if requested
    if (abort) ::abort();
}

class FDRedirect {
    int saved_out{-1}, saved_err{-1}, out_fd{-1}, err_fd{-1};
public:
    FDRedirect(const char* out_path, const char* err_path=nullptr) {
        saved_out = dup(STDOUT_FILENO);
        saved_err = dup(STDERR_FILENO);
        out_fd = ::open(out_path, O_WRONLY|O_CREAT|O_APPEND, 0644);
        if (!err_path) err_path = out_path;
        err_fd = ::open(err_path, O_WRONLY|O_CREAT|O_APPEND, 0644);
        dup2(out_fd, STDOUT_FILENO);
        dup2(err_fd, STDERR_FILENO);
    }
    ~FDRedirect() {
        fsync(STDOUT_FILENO);
        fsync(STDERR_FILENO);
        dup2(saved_out, STDOUT_FILENO); close(saved_out); close(out_fd);
        dup2(saved_err, STDERR_FILENO); close(saved_err); close(err_fd);
    }
};


// --------------------------- //
// Xcoll BDSIM implementations //
// --------------------------- //

BDSParticleDefinition* PrepareBDSParticleDefition(long long int pdgIDIn, double momentumIn, 
                                                  double kineticEnergyIn, double ionChargeIn)
{
    G4int pdgID = (G4int) pdgIDIn;

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefGeant;

    // Wrap in our class that calculates momentum and kinetic energy.
    // Requires that one of E, Ek, P be non-zero (only one).
    BDSParticleDefinition* particleDefinition = nullptr;
    BDSIonDefinition* ionDef = nullptr;

    // PDG for ions = 10LZZZAAAI
    if (pdgID > 1000000000) // is an ion
    {   

        G4IonTable* ionTable = particleTable->GetIonTable();
        particleDefGeant = ionTable->GetIon(pdgID);
        if (!particleDefGeant)
        {throw BDSException("BDSXtrackInterface> Ion \"" + std::to_string(pdgID) + "\" not found");}

        G4int ionCharge = ionChargeIn==0 ? particleDefGeant->GetAtomicNumber() : (G4int) ionChargeIn; 
        ionDef = new BDSIonDefinition(particleDefGeant->GetAtomicMass(),
                                      particleDefGeant->GetAtomicNumber(),
                                      ionCharge);

        // correct the particle mass in the case of partially-stripped ions
        G4double mass   = ionTable->GetIonMass(ionDef->Z(), ionDef->A());
        G4double charge = ionDef->Charge(); // correct even if overridden
        mass += ionDef->NElectrons()*G4Electron::Definition()->GetPDGMass();
        // The constructor with a custom name requires the bdsim name of the ion
        std::string bdsimPartName = "ion " + std::to_string(static_cast<int>(ionDef->A())) 
                                  + " " + std::to_string(static_cast<int>(ionDef->Z())) 
                                  + " " + std::to_string(static_cast<int>(charge));

        particleDefinition = new BDSParticleDefinition(bdsimPartName, mass, charge, 0, 
                                                       kineticEnergyIn, momentumIn, 1, ionDef, pdgID);
    }
    else
    {
        particleDefGeant = particleTable->FindParticle(pdgID);
        if (!particleDefGeant)
        {throw BDSException("BDSXtrackInterface> Particle \"" + std::to_string(pdgID) + "\" not found");}
        particleDefinition = new BDSParticleDefinition(particleDefGeant, 0, kineticEnergyIn, momentumIn, 1, nullptr);
    }

    return particleDefinition;
}

XtrackInterface::XtrackInterface(const  std::string& bdsimConfigFile,
                                 long long int       referencePdgIdIn,
                                 double              referenceEkIn,
                                 double              relativeEnergyCutIn,
                                 int                 seedIn,
                                 int                 referenceIonChargeIn,
                                 bool                batchMode):
        relativeEnergyCut(relativeEnergyCutIn),
        seed(seedIn)
{
    SetErrorHandler(MyRootErrorHandler);
    RedirectGeant4();
    fdredir = std::make_unique<FDRedirect>("engine.out","engine.err");

    stp = new BDSBunchSixTrackLink();
    bds = new BDSIMLink(stp);

	for (char* arg : argv){
		free(arg);  // Free dynamically allocated memory
	}
	argv.clear();  // Clear the vector

    std::string seedStr = std::to_string(seed);
    std::vector<std::string> arguments = {"--verbose",
                                          "--file=" + bdsimConfigFile,
                                          //"--file=" + bdsimConfigFile,
                                          //"--vis_debug",
                                          "--output=none",
                                          "--seed=" + seedStr,
                                          "--outfile=output_" + seedStr};

    for(auto & argument : arguments){
        argv.push_back(strdup(argument.c_str()));
    }
    if (batchMode){
        std::string batch_flag = "--batch";
        argv.push_back(strdup(batch_flag.c_str()));
        argv.push_back(strdup(batch_flag.c_str()));
        argv.push_back(strdup(batch_flag.c_str()));
        argv.push_back(strdup(batch_flag.c_str()));
        argv.push_back(strdup(batch_flag.c_str()));
    }
    // argv.push_back(nullptr);

    double referenceEk = referenceEkIn * CLHEP::GeV;

    double relEKCut = relativeEnergyCut;
    if (relEKCut < 1e-6) // defaults to 0 which means 0eV cut which is bad
    { relEKCut = 1.0; }

    double minimumEK = relEKCut * (referenceEk);

    G4cout << "Minimum kinetic energy " << minimumEK << " MeV" << G4endl;
    auto data = argv.data();

	// Print arguments
	std::cout << "Initialise called with arguments: " << std::endl;
	std::cout << "argc: " << argv.size() - 1 << std::endl;
	for (size_t i = 0; i < argv.size(); i++) {
		std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
	}
	G4cout << "minimumEK / CLHEP::GeV: " << minimumEK / CLHEP::GeV << G4endl;
	std::cout.flush();

    try
    { bds->Initialise(argv.size(), &argv[0], true, minimumEK / CLHEP::GeV, false); } // minimumEk in GeV
    catch (const std::exception &e)
    {
        std::cout << e.what() << std::endl;
        exit(1);
    }

    G4double ionCharge = (G4double) referenceIonChargeIn;
    refParticleDefinition = PrepareBDSParticleDefition(referencePdgIdIn, 0, referenceEk, ionCharge);
}


XtrackInterface::~XtrackInterface()
{

	// Clean up dynamically allocated memory in argv
	for (char* arg : argv) {
		free(arg);
	}
	argv.clear();  // Clear the vector
    //std::system("reset"); // spawn new c++ process
    //std::cout << "reset" << std::endl;
	// Clean up other resources
    delete bds;
    delete stp;
    delete refParticleDefinition;
    fdredir.reset();
}


void XtrackInterface::addCollimator(const std::string&   name,
                                    const std::string&   material,
                                    double lengthIn,
                                    double apertureLeftIn,
                                    double apertureRightIn,
                                    double rotationIn,
                                    double xOffsetIn,
                                    double yOffsetIn,
                                    double jawTiltLeft,
                                    double jawTiltRight,
                                    int    side,
                                    bool isACrystal)
    {

        bool buildLeft  = side == 0 || side == 1;
        bool buildRight = side == 0 || side == 2;

        bds->AddLinkCollimatorJaw(name,
                                  material,
                                  lengthIn * CLHEP::m,
                                  apertureLeftIn * CLHEP::m,
                                  apertureRightIn * CLHEP::m,
                                  rotationIn * CLHEP::rad,
                                  xOffsetIn * CLHEP::m,
                                  yOffsetIn * CLHEP::m,
                                  jawTiltLeft * CLHEP::rad,
                                  jawTiltRight * CLHEP::rad,
                                  buildLeft,
                                  buildRight,
                                  isACrystal,
                                  0);
    }


void XtrackInterface::addParticle(double xIn,
                                  double yIn,
                                  double pxIn,
                                  double pyIn,
                                  double ctIn,
                                  double deltapIn,
                                  double chiIn,
                                  double chargeRatioIn,
                                  double sIn,
                                  int64_t pdgIDIn,
                                  int64_t trackidIn
                                  )
{

    auto x  = (G4double) xIn;
    auto y  = (G4double) yIn;
    auto px  = (G4double) pxIn;
    auto py  = (G4double) pyIn;
    auto ct  = (G4double) ctIn;
    auto deltap = (G4double) deltapIn;
    auto chi = (G4double) chiIn;
    auto charge_ratio = (G4double) chargeRatioIn;
    auto s = (G4double) sIn;
    auto trackid = (G4int) trackidIn;
    auto pdgid = (long long int) pdgIDIn;

    //G4double mass_ratio = charge_ratio / chi;
    G4double q = charge_ratio * refParticleDefinition->Charge();
    G4double mass_ratio = charge_ratio / chi;
    G4double p = refParticleDefinition->Momentum() * (deltap + 1) * mass_ratio;

    // If the pdg id for the particle is 0 (default), take the reference pdg id
    long long int pdgidPart = (pdgid==0) ? refParticleDefinition->PDGID() : pdgid; 
    auto partDef = PrepareBDSParticleDefition(pdgidPart, p, 0, q);

    G4double t = - ct * CLHEP::m / (refParticleDefinition->Beta() * CLHEP::c_light); // this is time difference in ns

    G4double oneplusdelta = (1 + deltap);
    // G4double pz = std::sqrt(oneplusdelta*oneplusdelta - px*px - py*py);
    G4double xp = px / oneplusdelta;
    G4double yp = py / oneplusdelta;

    // Zp0 is 1 as here we assume no back-scatterd particles, e.g p>0
    G4double zp = BDSBunch::CalculateZp(xp, yp, 1);

    BDSParticleCoordsFull coords = BDSParticleCoordsFull(x * CLHEP::m,
                                                         y * CLHEP::m,
                                                         0,
                                                         xp,
                                                         yp,
                                                         zp,
                                                         t,
                                                         0,
                                                         partDef->TotalEnergy(),
                                                         1);

    stp->AddParticle(partDef, coords, trackid, trackid);
}


void XtrackInterface::addParticles(const py::list& coordinates)
{
    //TODO get the charge and mass ratios
    // Obtain the arrays from the list and cast them to the correct array type
    py::array_t<double> x = py::cast<py::array>(coordinates[0]);
    py::array_t<double> y = py::cast<py::array>(coordinates[1]);
    py::array_t<double> px = py::cast<py::array>(coordinates[2]);
    py::array_t<double> py = py::cast<py::array>(coordinates[3]);
    py::array_t<double> zeta = py::cast<py::array>(coordinates[4]);
    py::array_t<double> delta = py::cast<py::array>(coordinates[5]);
    py::array_t<double> chi = py::cast<py::array>(coordinates[6]);
    py::array_t<double> charge_ratio = py::cast<py::array>(coordinates[7]);
    py::array_t<double> s = py::cast<py::array>(coordinates[8]);
    py::array_t<int64_t> pdgid = py::cast<py::array>(coordinates[9]);
    py::array_t<int64_t> trackid = py::cast<py::array>(coordinates[10]);
    py::array_t<int64_t> state = py::cast<py::array>(coordinates[11]);
    py::array_t<int64_t> at_element = py::cast<py::array>(coordinates[12]);
    py::array_t<int64_t> at_turn = py::cast<py::array>(coordinates[13]);


    // Obtain the buffers
    py::buffer_info x_buff = x.request();
    py::buffer_info y_buff = y.request();
    py::buffer_info px_buff = px.request();
    py::buffer_info py_buff = py.request();
    py::buffer_info zeta_buff = zeta.request();
    py::buffer_info delta_buff = delta.request();
    py::buffer_info chi_buff = chi.request();
    py::buffer_info charge_ratio_buff = charge_ratio.request();
    py::buffer_info s_buff = s.request();
    py::buffer_info pdgid_buff = pdgid.request();
    py::buffer_info id_buff = trackid.request();
    py::buffer_info state_buff = state.request();
    py::buffer_info at_element_buff = at_element.request();
    py::buffer_info at_turn_buff = at_turn.request();

    // Get the pointers for iteration
    auto x_ptr = static_cast<double *>(x_buff.ptr);
    auto y_ptr = static_cast<double *>(y_buff.ptr);
    auto px_ptr = static_cast<double *>(px_buff.ptr);
    auto py_ptr = static_cast<double *>(py_buff.ptr);
    auto zeta_ptr = static_cast<double *>(zeta_buff.ptr);
    auto delta_ptr = static_cast<double *>(delta_buff.ptr);
    auto chi_ptr = static_cast<double *>(chi_buff.ptr);
    auto charge_ratio_ptr = static_cast<double *>(charge_ratio_buff.ptr);
    auto s_ptr = static_cast<double *>(s_buff.ptr);
    auto pdgid_ptr = static_cast<int64_t *>(pdgid_buff.ptr);
    auto trackid_ptr = static_cast<int64_t *>(id_buff.ptr);
    auto state_ptr = static_cast<int64_t *>(state_buff.ptr);
    auto at_element_ptr = static_cast<int64_t *>(at_element_buff.ptr);
    auto at_turn_ptr = static_cast<int64_t *>(at_turn_buff.ptr);

    long n = 1; // Get the number of elements in the array, assume all arrays have the same number of elements
    for (auto r: x_buff.shape) {
        n *= r;
    }

    for(int i=0; i < n; i++)
    {
        auto x_part = x_ptr[i];
        auto y_part = y_ptr[i];
        auto px_part = px_ptr[i];
        auto py_part = py_ptr[i];
        auto zeta_part = zeta_ptr[i];
        auto delta_part = delta_ptr[i];
        auto chi_part = chi_ptr[i];
        auto charge_ratio_part = charge_ratio_ptr[i];
        auto s_part = s_ptr[i];
        auto pdgid_part = pdgid_ptr[i];
        auto trackid_part = trackid_ptr[i];
        auto state_part = state_ptr[i];
        auto at_element_part = at_element_ptr[i];
        auto at_turn_part = at_turn_ptr[i];

        // The particles with state=0 are inactive and should not be processed at all, but need to keep track of them
        // The internal processing in BDSIM does not feature an active state check, so this must be done at a higher
        // level here

        if (state_part == 1) // State == 1 means that the particle is active
        {
            particleActiveState.push_back(true);
            addParticle(x_part, y_part, px_part, py_part,
                        zeta_part, delta_part, chi_part,
                        charge_ratio_part, s_part, pdgid_part, trackid_part);

            maxParticleID = std::max(maxParticleID, trackid_part);
            bds->SetCurrentMaximumExternalParticleID(maxParticleID);
        }
        else
        {
            if (state_part == -999999999) // This is reserve memory for new particles, no need to process
            {
                break;
            }
            particleActiveState.push_back(false);
        }

        auto particle_coords = new XtrackCoordinates{x_part, y_part, px_part, py_part, zeta_part,
                                                     delta_part, chi_part, charge_ratio_part, s_part,
                                                     pdgid_part, trackid_part, state_part,
                                                     at_element_part, at_turn_part};
        xtrackParticles.push_back(particle_coords);
    }
}


void XtrackInterface::collimate()
{
    bds->BeamOn((G4int)stp->Size());
}


void XtrackInterface::selectCollimator(const std::string& collimatorName)
{
    currentCollimatorName = collimatorName;
    // This doesn't throw an error if the element doesn't exist
    bds->SelectLinkElement(collimatorName);

    // Check if the element exists by querying the index: -1 means it doesn't exist
    if (bds->GetLinkIndex(collimatorName) == -1)
        {throw std::runtime_error("Element not found " + collimatorName);}
}


void XtrackInterface::clearData()
{
    bds->ClearSamplerHits();
    // A malloc error about freeing a pointer not allocated is thrown if the run is terminated after
    // the bunch is manually cleared. Consider using an alternative to vector.clear() in the BDSIM bunch class:
    // vector<T>().swap(x);   // clear x reallocating  (https://www.cplusplus.com/reference/vector/vector/clear/)
    // to ensure the pointer is still allocated after clearing
    stp->ClearParticles();

    for (auto part : xtrackParticles)
    {
       delete part;
    }

    std::vector<XtrackCoordinates*>().swap(xtrackParticles);
    std::vector<bool>().swap(particleActiveState);
    currentCollimatorName.clear();

    maxParticleID = 0;

}


py::dict XtrackInterface::collimateReturn(const py::list& coordinates)
{
    // Prepare the buffers for modifying the primary particle coordinates in place

    // Obtain the arrays from the list and cast them to the correct array type
    py::array_t<double> x = py::cast<py::array>(coordinates[0]);
    py::array_t<double> y = py::cast<py::array>(coordinates[1]);
    py::array_t<double> px = py::cast<py::array>(coordinates[2]);
    py::array_t<double> py = py::cast<py::array>(coordinates[3]);
    py::array_t<double> zeta = py::cast<py::array>(coordinates[4]);
    py::array_t<double> delta = py::cast<py::array>(coordinates[5]);
    py::array_t<double> chi = py::cast<py::array>(coordinates[6]);
    py::array_t<double> charge_ratio = py::cast<py::array>(coordinates[7]);
    py::array_t<double> s = py::cast<py::array>(coordinates[8]);
    py::array_t<int64_t> pdgid = py::cast<py::array>(coordinates[9]);
    py::array_t<int64_t> trackid = py::cast<py::array>(coordinates[10]);
    py::array_t<int64_t> state = py::cast<py::array>(coordinates[11]);
    py::array_t<int64_t> at_element = py::cast<py::array>(coordinates[12]);
    py::array_t<int64_t> at_turn = py::cast<py::array>(coordinates[13]);

    // Obtain the buffers
    py::buffer_info x_buff = x.request();
    py::buffer_info y_buff = y.request();
    py::buffer_info px_buff = px.request();
    py::buffer_info py_buff = py.request();
    py::buffer_info zeta_buff = zeta.request();
    py::buffer_info delta_buff = delta.request();
    py::buffer_info chi_buff = chi.request();
    py::buffer_info charge_ratio_buff = charge_ratio.request();
    py::buffer_info s_buff = s.request();
    py::buffer_info pdgid_buff = pdgid.request();
    py::buffer_info id_buff = trackid.request();
    py::buffer_info state_buff = state.request();
    py::buffer_info at_element_buff = at_element.request();
    py::buffer_info at_turn_buff = at_turn.request();

    // Get the pointers for iteration
    auto x_ptr = static_cast<double *>(x_buff.ptr);
    auto y_ptr = static_cast<double *>(y_buff.ptr);
    auto px_ptr = static_cast<double *>(px_buff.ptr);
    auto py_ptr = static_cast<double *>(py_buff.ptr);
    auto zeta_ptr = static_cast<double *>(zeta_buff.ptr);
    auto delta_ptr = static_cast<double *>(delta_buff.ptr);
    auto chi_ptr = static_cast<double *>(chi_buff.ptr);
    auto charge_ratio_ptr = static_cast<double *>(charge_ratio_buff.ptr);
    auto s_ptr = static_cast<double *>(s_buff.ptr);
    auto pdgid_ptr = static_cast<int64_t *>(pdgid_buff.ptr);
    auto trackid_ptr = static_cast<int64_t *>(id_buff.ptr);
    auto state_ptr = static_cast<int64_t *>(state_buff.ptr);
    auto at_element_ptr = static_cast<int64_t *>(at_element_buff.ptr);
    auto at_turn_ptr = static_cast<int64_t *>(at_turn_buff.ptr);

    // N.B book keeping needed
    // Note: Load and and transfer for surviving primaries
    // Secondaries inherit the element and turn of parent particle
    //auto at_element = static_cast<int64_t *>(at_element.ptr);
    //auto turn = static_cast<int64_t *>(turn.ptr);


    // Prepare the buffers for writing out the products
    // New numpy arrays are allocated for writing the products

    // Access the sampler hits - particles reaching the planes for transport back
    const BDSHitsCollectionSamplerLink* hits = bds->SamplerHits();

    //size_t hitsCount = hits ? hits->GetSize() : 0;

    size_t hitsCount = 0;
    if (hits)
    {
        hitsCount = hits->GetSize();
    }
    else
    {
        // There were no hits - check if there were any active particles at all coming in
        if (!stp->Size())
        {
            // A particle needs to be added to the bunch and the GetNextParticleLocal method
            // must be called to initialise the variables and ensure a safe deletion of the
            // underlying BDSBunch object. This will not be needed in the next release of BDSIM
            addParticle(0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0);
            // The dummy particle added is never used, as the collimateReturn method terminates
            // the processing
            stp->GetNextParticleLocal();
        }
    }


    // Count the number of secondary particles
    size_t secondaryCount = 0;
    for (size_t i = 0; i < hitsCount; i++)
    {
        auto hit = (*hits)[i];
        if (hit->externalParticleID != hit->externalParentID) { secondaryCount++; }
    }

    // The output arrays has slots for all primary particles, regardless if lost or not, and for secondary particles
    // size_t output_size = secondaryCount;
    size_t output_size = x.size();

    // TODO: this is mostly default buffers, so there are probably simpler constructors to use
    // Prepare the numpy array that will be returned
    auto x_out = py::array(py::buffer_info(
            nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sizeof(double),     /* Size of one item */
            py::format_descriptor<double>::value, /* Buffer format */
            1,          /* How many dimensions? */
            { output_size },  /* Number of elements for each dimension */
            { sizeof(double) }  /* Strides for each dimension */
    ));

    auto y_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto px_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto py_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto zeta_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto delta_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto s_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto pdgid_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto trackid_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto state_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto at_element_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto at_turn_out = py::array(py::buffer_info(
            nullptr,
            sizeof(int64_t),
            py::format_descriptor<int64_t>::value,
            1,
            { output_size },
            { sizeof(int64_t) }
    ));

    auto massratio_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto chargeratio_out = py::array(py::buffer_info(
            nullptr,
            sizeof(double),
            py::format_descriptor<double>::value,
            1,
            { output_size },
            { sizeof(double) }
    ));

    auto x_prod_buf = x_out.request();
    auto y_prod_buf = y_out.request();
    auto px_prod_buf = px_out.request();
    auto py_prod_buf = py_out.request();
    auto zeta_prod_buf = zeta_out.request();
    auto delta_prod_buf = delta_out.request();
    auto s_prod_buf = s_out.request();
    auto pdgid_prod_buf = pdgid_out.request();
    auto trackid_prod_buf = trackid_out.request();
    auto state_prod_buf = state_out.request();
    auto at_element_prod_buf = at_element_out.request();
    auto at_turn_prod_buf = at_turn_out.request();
    auto massratio_prod_buf = massratio_out.request();
    auto chargeratio_prod_buf = chargeratio_out.request();

    auto *x_prod_ptr = (double *) x_prod_buf.ptr;
    auto *y_prod_ptr = (double *) y_prod_buf.ptr;
    auto *px_prod_ptr = (double *) px_prod_buf.ptr;
    auto *py_prod_ptr = (double *) py_prod_buf.ptr;
    auto *zeta_prod_ptr = (double *) zeta_prod_buf.ptr;
    auto *delta_prod_ptr = (double *) delta_prod_buf.ptr;
    auto *s_prod_ptr = (double *) s_prod_buf.ptr;
    auto *pdgid_prod_ptr = (int64_t *) pdgid_prod_buf.ptr;
    auto *trackid_prod_ptr = (int64_t *) trackid_prod_buf.ptr;
    auto *state_prod_ptr = (int64_t *) state_prod_buf.ptr;
    auto *at_element_prod_ptr = (int64_t *) at_element_prod_buf.ptr;
    auto *at_turn_prod_ptr = (int64_t *) at_turn_prod_buf.ptr;
    auto *massratio_prod_ptr = (double *) massratio_prod_buf.ptr;
    auto *chargeratio_prod_ptr = (double *) chargeratio_prod_buf.ptr;

    long n = 1;
    for (auto r: x_buff.shape) {
        n *= r;
    }

    for (int i=0; i<n; i++) {
        x_prod_ptr[i] = x_ptr[i];
        y_prod_ptr[i] = y_ptr[i];
        px_prod_ptr[i] = px_ptr[i];
        py_prod_ptr[i] = py_ptr[i];
        zeta_prod_ptr[i] = zeta_ptr[i];
        delta_prod_ptr[i] = delta_ptr[i];
        s_prod_ptr[i] = s_ptr[i];
        pdgid_prod_ptr[i] = pdgid_ptr[i];
        trackid_prod_ptr[i] = trackid_ptr[i];
        state_prod_ptr[i] = state_ptr[i];

        at_element_prod_ptr[i] = at_element_ptr[i];
        at_turn_prod_ptr[i] = at_turn_ptr[i];
        massratio_prod_ptr[i] = charge_ratio_ptr[i]/chi_ptr[i];
        chargeratio_prod_ptr[i] = charge_ratio_ptr[i];
    }

    // Loop through the particles in the *original* bunch - the primaries
    size_t hits_index = 0;
    bool prim_survied = false;
    // double sum_deltaplusone_sec = 0.0;
    double sum_secondary_energy = 0.0;

    size_t prod_write_index = particleActiveState.size();
    for (size_t i=0; i < particleActiveState.size(); i++){

        if (!particleActiveState.at(i)){
            continue; // This was an inactive particle that hasn't been processed, do not change it
        }

        auto xtrack_part = xtrackParticles.at(i); // Get the cached coordinates of the original xtrack particle

        auto part = stp->GetNextParticle(); // Advance through the bunch
        auto prim_part_id = stp->CurrentExternalParticleID(); // Get the ID of the primary particle

        // Now start looping over the hits - the particles to be returned to the tracker
        // These can be primary or secondary particles. Each primary can produce 0, 1, or 2+ products
        // The products need to be sorted to keep the array order - surviving primary particles are all
        // filled in first. If a primary didn't survive, keep the original coordinates and make it inactive.
        // The hits are ordered by primary event, so just need one loop.
        while (hits_index < hitsCount)
        {
            BDSHitSamplerLink* hit = (*hits)[hits_index];
            if (hit->externalParentID != prim_part_id) { // The hits corresponding to the current primary are exhausted
                break;
            }

            const BDSParticleCoordsFull &coords = hit->coords;

            double qratio = hit->charge / refParticleDefinition->Charge();
            double mratio = hit->mass / refParticleDefinition->Mass();

            double dp = (hit->momentum / mratio - refParticleDefinition->Momentum()) / refParticleDefinition->Momentum();

            double collLength = bds->GetArcLengthOfLinkElement(currentCollimatorName);
            /// Need to compensate for the geometry construction in BDSIM
            /// There is a safety margin that is added to the collimator legnth
            double collMargin = 2.5 * BDSSamplerCustom::ChordLength();
            double zt = refParticleDefinition->Beta() * CLHEP::c_light * ((collLength + collMargin) / (CLHEP::c_light * refParticleDefinition->Beta()) - coords.T);

            double oneplusdelta = (1 + dp) * mratio;

            auto track_id = hit->externalParticleID;
            auto parent_id = hit->externalParentID;
            auto pdg_id = hit->pdgID;

            if (track_id == parent_id){
                // This is a primary particle as its parent is itself
                prim_survied = true;

                x_prod_ptr[i] = coords.x / CLHEP::m;
                y_prod_ptr[i] = coords.y / CLHEP::m;
                px_prod_ptr[i] = coords.xp * oneplusdelta; // convert back to px proper
                py_prod_ptr[i] = coords.yp * oneplusdelta;
                zeta_prod_ptr[i] = zt / CLHEP::m;
                delta_prod_ptr[i] = dp;
                s_prod_ptr[i] = xtrack_part->s + collLength / CLHEP::m;
                pdgid_prod_ptr[i] = pdg_id;
                //trackid_ptr[i] = track_id; // Don't touch the primary particle id
                trackid_prod_ptr[i] = track_id;
                state_prod_ptr[i] = 1; // active

                at_element_prod_ptr[i] = xtrack_part->at_element; // active
                at_turn_prod_ptr[i] = xtrack_part->at_turn; // active
                massratio_prod_ptr[i] = mratio;
                chargeratio_prod_ptr[i] = qratio;
            }
            else
            {
				//if (dp > 1.01 ) {
				//std::cout << "in secondary particles" << std::endl;
				//std::cout << hits_index << std::endl;
				//std::cout << hitsCount << std::endl;
				//std::cout << particleActiveState.size() << std::endl;
				//std::cout << particleActiveState.at(hits_index) << std::endl;
				//std::cout << prod_write_index << std::endl;
//std::cout << "coords.x = " << coords.x / CLHEP::m << std::endl;
//std::cout << "coords.y = " << coords.y / CLHEP::m << std::endl;
//std::cout << "coords.xp = " << coords.xp << std::endl;
//std::cout << "coords.yp = " << coords.yp << std::endl;
//std::cout << "oneplusdelta = " << oneplusdelta << std::endl;
//std::cout << "zt = " << zt / CLHEP::m << std::endl;
//std::cout << "dp = " << dp << std::endl;
//std::cout << "xtrack_part->s = " << xtrack_part->s << std::endl;
//std::cout << "collLength = " << collLength / CLHEP::m << std::endl;
//std::cout << "pdg_id = " << pdg_id << std::endl;
//std::cout << "parent_id = " << parent_id << std::endl;
//std::cout << "at_element = " << xtrack_part->at_element << std::endl;
//std::cout << "at_turn = " << xtrack_part->at_turn << std::endl;
//std::cout << "massratio = " << mratio << std::endl;
//std::cout << "chargeratio = " << qratio << std::endl;
//};
                // Secondary particles are populated in newly allocated arrays
                x_prod_ptr[prod_write_index] = coords.x / CLHEP::m;
                y_prod_ptr[prod_write_index] = coords.y / CLHEP::m;
                px_prod_ptr[prod_write_index] = coords.xp * oneplusdelta; // convert back to px proper
                py_prod_ptr[prod_write_index] = coords.yp * oneplusdelta; // convert back to py proper;
                zeta_prod_ptr[prod_write_index] = zt / CLHEP::m;
                delta_prod_ptr[prod_write_index] = dp;
                s_prod_ptr[prod_write_index] = xtrack_part->s + collLength / CLHEP::m;
                pdgid_prod_ptr[prod_write_index] = pdg_id;
                trackid_prod_ptr[prod_write_index] = parent_id;
                state_prod_ptr[prod_write_index] = 1; // active
                at_element_prod_ptr[prod_write_index] = xtrack_part->at_element; // active
                at_turn_prod_ptr[prod_write_index] = xtrack_part->at_turn; // active
                massratio_prod_ptr[prod_write_index] = mratio;
                chargeratio_prod_ptr[prod_write_index] = qratio;

                sum_secondary_energy += std::sqrt(std::pow(hit->momentum,2) + std::pow(hit->mass,2));
                prod_write_index++;
            }

            hits_index++;
        }

        if (!prim_survied) // Primary didn't survive - set inactive
        {
            state_prod_ptr[i] = -333; // inactive
            // Correct the energy of the lost primary particle to account for the production of secondaries
            // The effective delta is such that the lost particle has the effective delta
            // which corresponds to the energy in - energy out for this primary

            // reconstruct the incoming primary particle energy
            G4double qprim = charge_ratio_ptr[i] * refParticleDefinition->Charge();
            G4double mass_ratio_prim = charge_ratio_ptr[i] / chi_ptr[i];
            G4double p_prim = refParticleDefinition->Momentum() * (delta_ptr[i] + 1) * mass_ratio_prim;
            G4double mass_prim = mass_ratio_prim * refParticleDefinition->Mass();
            G4double energy_prim = std::sqrt(std::pow(p_prim, 2) + std::pow(mass_prim, 2));

            // compute the effective delta
            G4double energy_diff = energy_prim - sum_secondary_energy;
            G4double p_eff;
            G4double squared_energy_mass_diff = std::pow(energy_diff, 2) - std::pow(mass_prim, 2);
            if (squared_energy_mass_diff < 0){
                // This means that the total energy escaping includes part of the rest mass of 
                // the primary. Tolerate the error for now, as otherwise need to adjust also the
                // mass and PDG id of the lost primary particle
                p_eff = 0;
            }
            else
            {
                p_eff = std::sqrt(squared_energy_mass_diff);
            }
            G4double delta_eff = (p_eff / mass_ratio_prim - refParticleDefinition->Momentum()) / refParticleDefinition->Momentum();
            delta_prod_ptr[i] = delta_eff;
        }
        prim_survied = false; // reset for next particle
        sum_secondary_energy = 0.0;
    }

    auto result = py::dict();

    result["s"] = s_out;
    result["x"] = x_out;
    result["px"] = px_out;
    result["y"] = y_out;
    result["py"] = py_out;
    result["zeta"] = zeta_out;
    result["delta"] = delta_out;
    result["pdg_id"] = pdgid_out;
    result["at_element"] = at_element_out;
    result["at_turn"] = at_turn_out;
    result["mass_ratio"] = massratio_out;
    result["charge_ratio"] = chargeratio_out;
    result["parent_particle_id"] = trackid_out;
    result["state"] = state_out;

    return result;
}
