#include <pybind11/pybind11.h>

#include "utilities.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pycali, m)
{
  py::class_<Config>(m, "Config")
    .def(py::init([]() {return new Config(); }))
    .def(py::init([](const std::string& fname) {return new Config(fname);}));

  py::class_<Cali>(m, "Cali")
    .def(py::init<>())
    .def(py::init([](Config& cfg) {return new Cali(cfg);}))
    .def("mcmc", &Cali::mcmc)
    .def("get_best_params", &Cali::get_best_params)
    .def("align_with_error", &Cali::align_with_error)
    .def("output", &Cali::output)
    .def("recon", &Cali::recon);
}