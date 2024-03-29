#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "utilities.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pycali, m)
{
  py::class_<Config>(m, "Config")
    .def(py::init([]() {return new Config(); }))
    .def(py::init([](const std::string& fname) {return new Config(fname);}))
    .def("get_param_filename", &Config::get_param_filename)
    .def("setup", &Config::setup, py::arg("fcont"), py::arg("fline")=list<string>({}),
                  py::arg("nmcmc")=2000, py::arg("ptol")=0.7, 
                  py::arg("scale_range_low")=0.5, py::arg("scale_range_up")=1.5,
                  py::arg("shift_range_low")=-1.0, py::arg("shift_range_up")=1.0,
                  py::arg("syserr_range_low")=0.0, py::arg("syserr_range_up")=0.1,
                  py::arg("errscale_range_low")=0.1, py::arg("errscale_range_up")=2.0,
                  py::arg("sigma_range_low")=1.0e-4, py::arg("sigma_range_up")=1.0,
                  py::arg("tau_range_low")=1.0, py::arg("tau_range_up")=1.0e4,
                  py::arg("fixed_scale")=false, py::arg("fixed_shift")=false,
                  py::arg("fixed_syserr")=true, py::arg("fixed_error_scale")=true,
                  py::arg("fixed_codes")=vector<int>({}),
                  py::arg("fixed_scalecodes")=vector<int>({}),
                  py::arg("flag_norm")=true
                  )
    .def("print_cfg", &Config::print_cfg)
    .def("test", &Config::test, py::arg("range")=vector<double>({0.5, 1.5}))
    .def_readwrite("fcont", &Config::fcont)
    .def_readwrite("fline", &Config::fline)
    .def_readwrite("nmcmc", &Config::nmcmc)
    .def_readwrite("scale_range_low", &Config::scale_range_low)
    .def_readwrite("scale_range_up", &Config::scale_range_up)
    .def_readwrite("shift_range_low", &Config::shift_range_low)
    .def_readwrite("shift_range_up", &Config::shift_range_up)
    .def_readwrite("syserr_range_low", &Config::syserr_range_low)
    .def_readwrite("syserr_range_up", &Config::syserr_range_up)
    .def_readwrite("errscale_range_low", &Config::errscale_range_low)
    .def_readwrite("errscale_range_up", &Config::errscale_range_up)
    .def_readwrite("sigma_range_low", &Config::sigma_range_low)
    .def_readwrite("sigma_range_up", &Config::sigma_range_up)
    .def_readwrite("tau_range_low", &Config::tau_range_low)
    .def_readwrite("tau_range_up", &Config::tau_range_up)
    .def_readwrite("fixed_scale", &Config::fixed_scale)
    .def_readwrite("fixed_shift", &Config::fixed_shift)
    .def_readwrite("fixed_syserr", &Config::fixed_syserr)
    .def_readwrite("fixed_error_scale", &Config::fixed_error_scale)
    .def_readwrite("fixed_codes", &Config::fixed_codes)
    .def_readwrite("fixed_scalecodes", &Config::fixed_scalecodes)
    .def_readwrite("flag_norm", &Config::flag_norm);

  py::class_<Cali>(m, "Cali")
    .def(py::init<>())
    .def(py::init([](Config& cfg) {return new Cali(cfg);}))
    .def("mcmc", &Cali::mcmc)
    .def("get_best_params", &Cali::get_best_params)
    .def("align_with_error", &Cali::align_with_error)
    .def("output", &Cali::output)
    .def("recon", &Cali::recon)
    .def("get_norm_cont", &Cali::get_norm_cont)
    .def("get_norm_line", &Cali::get_norm_line)
    .def_readwrite("ncode", &Cali::ncode)
    .def_readwrite("num_params", &Cali::num_params);
}