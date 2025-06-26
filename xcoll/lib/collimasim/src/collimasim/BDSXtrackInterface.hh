#include "BDSBunchSixTrackLink.hh"
#include "BDSException.hh"
#include "BDSIMLink.hh"
#include "BDSIonDefinition.hh"
#include "BDSParticleCoordsFull.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"

#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Types.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// The struct is only used for inactive particle coodrinates for now
struct XtrackCoordinates{
    double x;
    double y;
    double px;
    double py;
    double zeta;
    double delta;
    double chi;
    double charge_ratio;
    double s;
    int64_t pdgid;
    int64_t trackid;
    int64_t state;
    int64_t at_element;
    int64_t at_turn;
};

BDSParticleDefinition* PrepareBDSParticleDefition(long long int pdgIDIn, double momentumIn, 
                                                  double kineticEnergyIn, double ionChargeIn);

/// TODO: make a base class for the interface classes as there is a lot of shared functionality
class XtrackInterface
{
public:
    XtrackInterface() = delete;  // No default constructor

    XtrackInterface(const std::string&  bdsimConfigFile,
                    long long int       referencePdgIdIn,
                    double              referenceEkIn,
                    double              relativeEnergyCutIn,
                    int                 seedIn,
                    int                 referenceIonChargeIn=0,
                    bool                batchMode=true);

    virtual ~XtrackInterface();

    void addCollimator(const  std::string&   name,
                       const  std::string&   material,
                       double lengthIn,
                       double apertureLeftIn,
                       double apertureRightIn,
                       double rotationIn,
                       double xOffsetIn,
                       double yOffsetIn,
                       double jawTiltLeft,
                       double jawTiltRight,
                       int    side);

    void addParticles(const py::list& coordinates);

    // This is not exposed to python
    void addParticle(double xIn,
                     double yIn,
                     double pxIn,
                     double pyIn,
                     double ctIn,
                     double deltapIn,
                     double chiIn,
                     double chargeRatioIn,
                     double sIn,
                     int64_t pdgIDIn,
                     int64_t trackidIn);

    void collimate();
    void clearData();
    void selectCollimator(const std::string& name);

    double getReferenceMass() { return refParticleDefinition->Mass() / CLHEP::GeV; }

    py::dict collimateReturn(const py::list& coordinates);

private:
    BDSIMLink* bds = nullptr;
    BDSBunchSixTrackLink* stp = nullptr;
    std::vector<char *> argv;

    std::vector<bool> particleActiveState;
    std::vector<XtrackCoordinates*> xtrackParticles;

    BDSParticleDefinition* refParticleDefinition = nullptr;

    std::string currentCollimatorName;

    int64_t maxParticleID=0;
    double relativeEnergyCut=0.0;
    int seed=0;
};