#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>

#include "BDSXtrackInterface.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(g4interface, m) {
  m.doc() = R"pbdoc(Python interface to BDSIM (Geant4). The main purpose of the interface is to enable collimation studies, including particle-matter interaction in collimators, for pure particle tracking codes.)pbdoc";

  py::class_<XtrackInterface>(m, "XtrackInterface")
            .def(py::init<const std::string&, int, double, double, int, int, bool>(),
                 "bdsimConfigFile"_a, "referencePdgId"_a, "referenceEk"_a,
                 "relativeEnergyCut"_a, "seed"_a, "referenceIonCharge"_a=0,"batchMode"_a=true,
                 py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
            .def("addCollimator", &XtrackInterface::addCollimator,
                 "name"_a, "material"_a, "tipMaterial"_a, "tipThickness"_a, "length"_a, "apertureLeft"_a,
                 "apertureRight"_a, "rotation"_a, "xOffset"_a, "yOffset"_a, "jawTiltLeft"_a,
                 "jawTiltRight"_a, "side"_a, "isACrystal"_a,
                 py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
            /// The cast here to disambiguate overloaded methods. The other versions of
            /// addParticles were removed, but keep the cast for reference for now.
            .def("addParticles", static_cast<void (XtrackInterface::*)
            (const py::list& coordinates)>(&XtrackInterface::addParticles), "coordinates"_a)
            /// The C++ ostream and estream are redirected to the Python streams
            /// so all io can be handled on the Python side
            .def("collimate", &XtrackInterface::collimate,
                 py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
            .def("selectCollimator", &XtrackInterface::selectCollimator, "collimatorName"_a)
            .def("collimateReturn", &XtrackInterface::collimateReturn,
                 "num_sent"_a,
                 py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
            .def("clearData", &XtrackInterface::clearData)
            .def("getReferenceMass", &XtrackInterface::getReferenceMass);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
