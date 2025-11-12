#include "BDSXtrackInterface.hh"
#include "xcoll_logging.h"

#include <cstring>
#include <cmath>
#include <BDSSamplerCustom.hh>
#include <iostream>
#include <vector>
#include <string>


BDSParticleDefinition* PrepareBDSParticleDefition(long long int pdgIDIn, double momentumIn,
                                                  double kineticEnergyIn, double ionChargeIn){
    G4int pdgID = (G4int) pdgIDIn;

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDefGeant;

    // Wrap in our class that calculates momentum and kinetic energy.
    // Requires that one of E, Ek, P be non-zero (only one).
    BDSParticleDefinition* particleDefinition = nullptr;
    BDSIonDefinition* ionDef = nullptr;

    // PDG for ions = 10LZZZAAAI
    if (pdgID > 1000000000){ // is an ion

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
    } else {
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
        seed(seedIn){
    // Redirect ROOT messages to file
    rootlog::init("root.out", "root.err", /*append=*/true);
    // Redirect Geant4 messages to file
    RedirectGeant4();
    // Redirect stdout and stderr to file (everything else like BDSIM)
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

    double referenceEk = referenceEkIn * CLHEP::eV;
    double relEKCut = relativeEnergyCut;
    // if (relEKCut < 1e-6)
    // { relEKCut = 1.0; }
    double minimumEK = relEKCut * referenceEk;
    G4cout << "Reference kinetic energy " << referenceEk << " MeV" << G4endl;
    G4cout << "Minimum kinetic energy " << minimumEK << " MeV" << G4endl;

	// Print arguments
    auto data = argv.data();
	std::cout << "Initialise called with arguments: " << std::endl;
	std::cout << "argc: " << argv.size() - 1 << std::endl;
	for (size_t i = 0; i < argv.size(); i++) {
		std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
	}
	G4cout << "minimumEK / CLHEP::GeV: " << minimumEK / CLHEP::GeV << G4endl;
	std::cout.flush();

    try
    { bds->Initialise(argv.size(), &argv[0], true, minimumEK / CLHEP::GeV, false); } // Initialise takes minimumEk in GeV
    catch (const std::exception &e)
    {
        std::cout << e.what() << std::endl;
        exit(1);
    }

    G4double ionCharge = (G4double) referenceIonChargeIn;
    refParticleDefinition = PrepareBDSParticleDefition(referencePdgIdIn, 0, referenceEk, ionCharge);
}


XtrackInterface::~XtrackInterface(){
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
                                    bool isACrystal){

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
                                  double xpIn,
                                  double yIn,
                                  double ypIn,
                                  double zetaIn,
                                  double pIn,
                                  double qIn,
                                  double weightIn,
                                  int64_t pdgIDIn,
                                  int64_t trackidIn
                                  ){

    auto x  = (G4double) xIn;
    auto xp = (G4double) xpIn;
    auto y  = (G4double) yIn;
    auto yp = (G4double) ypIn;
    auto zeta  = (G4double) zetaIn;
    auto p = (G4double) pIn;
    auto q = (G4double) qIn;
    auto weight = (G4double) weightIn;
    auto trackid = (G4int) trackidIn;
    auto pdgid = (long long int) pdgIDIn;

    // If the pdg id for the particle is 0 (default), take the reference pdg id
    long long int pdgidPart = (pdgid==0) ? refParticleDefinition->PDGID() : pdgid; 
    auto partDef = PrepareBDSParticleDefition(pdgidPart, p * CLHEP::eV, 0, q);

    G4double t = - zeta * CLHEP::m / (refParticleDefinition->Beta() * CLHEP::c_light); // this is time difference in ns

    // Zp0 is 1 as here we assume no back-scatterd particles, e.g p>0
    G4double zp = BDSBunch::CalculateZp(xp * CLHEP::rad, yp * CLHEP::rad , 1);

    BDSParticleCoordsFull coords = BDSParticleCoordsFull(x * CLHEP::m,
                                                         y * CLHEP::m,
                                                         0,
                                                         xp * CLHEP::rad,
                                                         yp * CLHEP::rad,
                                                         zp,
                                                         t,  // already in ns
                                                         0,
                                                         partDef->TotalEnergy(),
                                                         weight);

    stp->AddParticle(partDef, coords, trackid, trackid);
}


void XtrackInterface::addParticles(const py::list& coordinates){
    // Obtain the arrays from the list and cast them to the correct array type
    py::array_t<double> x = py::cast<py::array>(coordinates[0]);
    py::array_t<double> xp = py::cast<py::array>(coordinates[1]);
    py::array_t<double> y = py::cast<py::array>(coordinates[2]);
    py::array_t<double> yp = py::cast<py::array>(coordinates[3]);
    py::array_t<double> zeta = py::cast<py::array>(coordinates[4]);
    py::array_t<double> p = py::cast<py::array>(coordinates[5]);
    py::array_t<double> q = py::cast<py::array>(coordinates[6]);
    py::array_t<double> weight = py::cast<py::array>(coordinates[7]);
    py::array_t<int64_t> pdgid = py::cast<py::array>(coordinates[8]);
    py::array_t<int64_t> trackid = py::cast<py::array>(coordinates[9]);

    // Obtain the buffers
    py::buffer_info x_buff = x.request();
    py::buffer_info xp_buff = xp.request();
    py::buffer_info y_buff = y.request();
    py::buffer_info yp_buff = yp.request();
    py::buffer_info zeta_buff = zeta.request();
    py::buffer_info p_buff = p.request();
    py::buffer_info q_buff = q.request();
    py::buffer_info weight_buff = weight.request();
    py::buffer_info pdgid_buff = pdgid.request();
    py::buffer_info id_buff = trackid.request();

    // Get the pointers for iteration
    auto x_ptr = static_cast<double *>(x_buff.ptr);
    auto xp_ptr = static_cast<double *>(xp_buff.ptr);
    auto y_ptr = static_cast<double *>(y_buff.ptr);
    auto yp_ptr = static_cast<double *>(yp_buff.ptr);
    auto zeta_ptr = static_cast<double *>(zeta_buff.ptr);
    auto p_ptr = static_cast<double *>(p_buff.ptr);
    auto q_ptr = static_cast<double *>(q_buff.ptr);
    auto weight_ptr = static_cast<double *>(weight_buff.ptr);
    auto pdgid_ptr = static_cast<int64_t *>(pdgid_buff.ptr);
    auto trackid_ptr = static_cast<int64_t *>(id_buff.ptr);

    long n = 1; // Get the number of elements in the array, assume all arrays have the same number of elements
    for (auto r: x_buff.shape){
        n *= r;
    }

    for(int i=0; i < n; i++){
        auto x_part = x_ptr[i];
        auto xp_part = xp_ptr[i];
        auto y_part = y_ptr[i];
        auto yp_part = yp_ptr[i];
        auto zeta_part = zeta_ptr[i];
        auto p_part = p_ptr[i];
        auto q_part = q_ptr[i];
        auto weight_part = weight_ptr[i];
        auto pdgid_part = pdgid_ptr[i];
        auto trackid_part = trackid_ptr[i];

        addParticle(x_part, xp_part, y_part, yp_part,
                    zeta_part, p_part, q_part, weight_part,
                    pdgid_part, trackid_part);

        maxParticleID = std::max(maxParticleID, trackid_part);
        bds->SetCurrentMaximumExternalParticleID(maxParticleID);
    }
}


void XtrackInterface::collimate(){
    bds->BeamOn((G4int)stp->Size());
}


void XtrackInterface::selectCollimator(const std::string& collimatorName){
    currentCollimatorName = collimatorName;
    // This doesn't throw an error if the element doesn't exist
    bds->SelectLinkElement(collimatorName);

    // Check if the element exists by querying the index: -1 means it doesn't exist
    if (bds->GetLinkIndex(collimatorName) == -1)
        {throw std::runtime_error("Element not found " + collimatorName);}
}


void XtrackInterface::clearData(){
    bds->ClearSamplerHits();
    // A malloc error about freeing a pointer not allocated is thrown if the run is terminated after
    // the bunch is manually cleared. Consider using an alternative to vector.clear() in the BDSIM bunch class:
    // vector<T>().swap(x);   // clear x reallocating  (https://www.cplusplus.com/reference/vector/vector/clear/)
    // to ensure the pointer is still allocated after clearing
    stp->ClearParticles();
    currentCollimatorName.clear();
    maxParticleID = 0;
}


py::dict XtrackInterface::collimateReturn(size_t num_sent, size_t output_size){
    // Access the sampler hits - particles reaching the planes for transport back
    const BDSHitsCollectionSamplerLink* hits = bds->SamplerHits();

    size_t hitsCount = hits ? hits->GetSize() : 0;

    // Count the number of secondary particles
    size_t secondaryCount = 0;
    for (size_t i = 0; i < hitsCount; i++){
        auto hit = (*hits)[i];
        if (hit->externalParticleID != hit->externalParentID) {secondaryCount++;}
    }
    if (secondaryCount + num_sent > output_size){
        throw BDSException("BDSXtrackInterface> Too many secondaries!\nProduced " + std::to_string(secondaryCount) \
                         + ", but only space for " + std::to_string(output_size - num_sent) + "!");
    }

    auto x_out    = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto xp_out   = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto y_out    = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto yp_out   = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto zeta_out = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto p_out    = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto m_out    = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto q_out    = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {output_size}, {sizeof(double)}));
    auto weight_out  = py::array(py::buffer_info(nullptr, sizeof(double),  py::format_descriptor<double>::value, 1, {output_size},  {sizeof(double)}));
    auto pdgid_out   = py::array(py::buffer_info(nullptr, sizeof(int64_t), py::format_descriptor<int64_t>::value, 1, {output_size}, {sizeof(int64_t)}));
    auto trackid_out = py::array(py::buffer_info(nullptr, sizeof(int64_t), py::format_descriptor<int64_t>::value, 1, {output_size}, {sizeof(int64_t)}));
    auto state_out   = py::array(py::buffer_info(nullptr, sizeof(int64_t), py::format_descriptor<int64_t>::value, 1, {output_size}, {sizeof(int64_t)}));

    auto *x_prod_ptr = (double *) x_out.request().ptr;
    auto *xp_prod_ptr = (double *) xp_out.request().ptr;
    auto *y_prod_ptr = (double *) y_out.request().ptr;
    auto *yp_prod_ptr = (double *) yp_out.request().ptr;
    auto *zeta_prod_ptr = (double *) zeta_out.request().ptr;
    auto *p_prod_ptr = (double *) p_out.request().ptr;
    auto *m_prod_ptr = (double *) m_out.request().ptr;
    auto *q_prod_ptr = (double *) q_out.request().ptr;
    auto *weight_prod_ptr = (double *) weight_out.request().ptr;
    auto *pdgid_prod_ptr = (int64_t *) pdgid_out.request().ptr;
    auto *trackid_prod_ptr = (int64_t *) trackid_out.request().ptr;
    auto *state_prod_ptr = (int64_t *) state_out.request().ptr;

    // Loop through the particles in the *original* bunch - the primaries
    size_t hits_index = 0;
    bool primary_survived = false;
    double sum_secondary_energy = 0.0;

    double collLength = bds->GetArcLengthOfLinkElement(currentCollimatorName);
    /// Need to compensate for the geometry construction in BDSIM
    /// There is a safety margin that is added to the collimator legnth
    collLength += 2.5 * BDSSamplerCustom::ChordLength();
    double velocity = refParticleDefinition->Beta() * CLHEP::c_light;

    size_t prod_write_index = num_sent;
    for (size_t i=0; i < num_sent; i++){
        auto part = stp->GetNextParticle(); // Advance through the bunch
        auto primary_part_id = stp->CurrentExternalParticleID(); // Get the ID of the primary particle

        // Now start looping over the hits - the particles to be returned to the tracker
        // These can be primary or secondary particles. Each primary can produce 0, 1, or 2+ products
        // The products need to be sorted to keep the array order - surviving primary particles are all
        // filled in first. If a primary didn't survive, keep the original coordinates and make it inactive.
        // The hits are ordered by primary event, so just need one loop.
        while (hits_index < hitsCount){
            BDSHitSamplerLink* hit = (*hits)[hits_index];
            if (hit->externalParentID != primary_part_id) {
                // The hits corresponding to the current primary are exhausted
                break;
            }

            const BDSParticleCoordsFull &coords = hit->coords;
            double zeta = collLength - coords.T * velocity;
            auto track_id = hit->externalParticleID;
            auto parent_id = hit->externalParentID;

            if (track_id == parent_id){
                // This is a primary particle as its parent is itself
                primary_survived = true;
                x_prod_ptr[i] = coords.x / CLHEP::m;
                xp_prod_ptr[i] = coords.xp / CLHEP::rad;
                y_prod_ptr[i] = coords.y / CLHEP::m;
                yp_prod_ptr[i] = coords.yp / CLHEP::rad;
                zeta_prod_ptr[i] = zeta / CLHEP::m;
                p_prod_ptr[i] = hit->momentum / CLHEP::eV;
                m_prod_ptr[i] = hit->mass / CLHEP::eV;
                q_prod_ptr[i] = hit->charge;
                weight_prod_ptr[i] = coords.weight;
                pdgid_prod_ptr[i] = hit->pdgID;
                trackid_prod_ptr[i] = track_id;
                state_prod_ptr[i] = 1; // active

            } else {
                // Secondary particles are populated in newly allocated arrays
                x_prod_ptr[prod_write_index] = coords.x / CLHEP::m;
                xp_prod_ptr[prod_write_index] = coords.xp / CLHEP::rad;
                y_prod_ptr[prod_write_index] = coords.y / CLHEP::m;
                yp_prod_ptr[prod_write_index] = coords.yp / CLHEP::rad;
                zeta_prod_ptr[prod_write_index] = zeta / CLHEP::m;
                p_prod_ptr[prod_write_index] = hit->momentum / CLHEP::eV;
                m_prod_ptr[prod_write_index] = hit->mass / CLHEP::eV;
                q_prod_ptr[prod_write_index] = hit->charge;
                weight_prod_ptr[prod_write_index] = coords.weight;
                pdgid_prod_ptr[prod_write_index] = hit->pdgID;
                trackid_prod_ptr[prod_write_index] = parent_id;
                state_prod_ptr[prod_write_index] = 1; // active;

                prod_write_index++;
            }

            hits_index++;
        }

        if (!primary_survived){ // Primary didn't survive - set inactive
            state_prod_ptr[i] = -300; // inactive (LOST_WITHOUT_SPEC)
        }
        primary_survived = false; // reset for next particle
    }

    auto result = py::dict();

    result["x"] = x_out;
    result["xp"] = xp_out;
    result["y"] = y_out;
    result["yp"] = yp_out;
    result["zeta"] = zeta_out;
    result["p"] = p_out;
    result["m"] = m_out;
    result["q"] = q_out;
    result["weight"] = weight_out;
    result["pdg_id"] = pdgid_out;
    result["parent_particle_id"] = trackid_out;
    result["state"] = state_out;

    return result;
}
