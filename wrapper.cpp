#include <pybind11/pybind11.h>

#include "utilities.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pycali, m)
{
  py::class_<Config>(m, "Config")
    .def(py::init([]() {return new Config(); }))
    .def(py::init([](const std::string& fname) {return new Config(fname);}))
    .def("get_param_filename", &Config::get_param_filename)
    .def("setup", &Config::setup, py::arg("fcont"), py::arg("fline")="",
                  py::arg("nmcmc")=2000, py::arg("pdiff")=0.7, 
                  py::arg("scale_low")=0.5, py::arg("scale_up")=1.5,
                  py::arg("shift_low")=-1.0, py::arg("shift_up")=1.0,
                  py::arg("sigma_low")=1.0e-4, py::arg("sigma_up")=1.0,
                  py::arg("tau_low")=1.0, py::arg("tau_up")=1.0e4
                  )
    .def("print_cfg", &Config::print_cfg)
    .def_readwrite("fcont", &Config::fcont)
    .def_readwrite("fline", &Config::fline)
    .def_readwrite("nmcmc", &Config::nmcmc)
    .def_readwrite("scale_low", &Config::scale_low)
    .def_readwrite("scale_up", &Config::scale_up)
    .def_readwrite("shift_low", &Config::shift_low)
    .def_readwrite("shift_up", &Config::shift_up)
    .def_readwrite("sigma_low", &Config::sigma_low)
    .def_readwrite("sigma_up", &Config::sigma_up)
    .def_readwrite("tau_low", &Config::tau_low)
    .def_readwrite("tau_up", &Config::tau_up);

  py::class_<Cali>(m, "Cali")
    .def(py::init<>())
    .def(py::init([](Config& cfg) {return new Cali(cfg);}))
    .def("mcmc", &Cali::mcmc)
    .def("get_best_params", &Cali::get_best_params)
    .def("align_with_error", &Cali::align_with_error)
    .def("output", &Cali::output)
    .def("recon", &Cali::recon);
}