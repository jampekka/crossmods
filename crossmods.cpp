#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "crossmods.hpp"
namespace py = pybind11;
PYBIND11_MODULE(crossmods, m){
using namespace pybind11::literals;
py::class_<::crossmods::Grid1d>(m, "Grid1d")
.def_readwrite("_low", &::crossmods::Grid1d::_low)
.def_readwrite("dx", &::crossmods::Grid1d::dx)
.def_readwrite("N", &::crossmods::Grid1d::N)
.def(py::init<double, double, size_t>(), "low"_a, "high"_a, "N"_a)
.def("bin", (::size_t ( ::crossmods::Grid1d::* )( double )const) &::crossmods::Grid1d::bin, "x"_a)
.def("low", (double ( ::crossmods::Grid1d::* )(  )) &::crossmods::Grid1d::low)
.def("high", (double ( ::crossmods::Grid1d::* )(  )const) &::crossmods::Grid1d::high)
.def("__getitem__", (double ( ::crossmods::Grid1d::* )( ::size_t )const) &::crossmods::Grid1d::operator[], "i"_a)
.def(py::init<crossmods::Grid1d const &>(), "arg0"_a);
py::class_<::crossmods::CrossingPdf>(m, "CrossingPdf")
.def_readwrite("grid", &::crossmods::CrossingPdf::grid)
.def_readwrite("ps", &::crossmods::CrossingPdf::ps)
.def_readwrite("uncrossed", &::crossmods::CrossingPdf::uncrossed)
.def(py::init<size_t, double>(), "dur"_a, "dt"_a)
.def("__call__", (double ( ::crossmods::CrossingPdf::* )( double )) &::crossmods::CrossingPdf::operator(), "ct"_a)
.def("loglikelihood", (double ( ::crossmods::CrossingPdf::* )( ::crossmods::invec,double )) &::crossmods::CrossingPdf::loglikelihood, "cts"_a, "slack"_a=0.)
.def(py::init<crossmods::CrossingPdf const &>(), "arg0"_a);
py::class_<::crossmods::NormalDistribution>(m, "NormalDistribution")
.def_readwrite("mean", &::crossmods::NormalDistribution::mean)
.def_readwrite("std", &::crossmods::NormalDistribution::std)
.def(py::init<double, double>(), "mean"_a=0., "std"_a=1.)
.def("cdf", (double ( ::crossmods::NormalDistribution::* )( double )) &::crossmods::NormalDistribution::cdf, "x"_a)
.def(py::init<crossmods::NormalDistribution const &>(), "arg0"_a);
py::class_<::crossmods::LognormalDistribution>(m, "LognormalDistribution")
.def_readwrite("mean", &::crossmods::LognormalDistribution::mean)
.def_readwrite("std", &::crossmods::LognormalDistribution::std)
.def(py::init<double, double>(), "mean"_a=0., "std"_a=1.)
.def("cdf", (double ( ::crossmods::LognormalDistribution::* )( double )) &::crossmods::LognormalDistribution::cdf, "x"_a)
.def(py::init<crossmods::LognormalDistribution const &>(), "arg0"_a);
py::class_<::crossmods::LognormalTdm>(m, "LognormalTdm")
.def(py::init<double, double, double, double, double>(), "thm"_a, "ths"_a, "lagm"_a, "lags"_a, "pass_th"_a)
.def("decisions", (::crossmods::CrossingPdf & ( ::crossmods::LognormalTdm::* )( ::crossmods::CrossingPdf &,::crossmods::invec,::size_t )) &::crossmods::LognormalTdm::decisions, "out"_a, "taus"_a, "start"_a=0)
.def("blocker_decisions", (::crossmods::CrossingPdf & ( ::crossmods::LognormalTdm::* )( ::crossmods::CrossingPdf &,::crossmods::invec,::crossmods::invec )) &::crossmods::LognormalTdm::blocker_decisions, "out"_a, "taus"_a, "taus_b"_a)
.def("blocker_decisions", (::std::unique_ptr<crossmods::CrossingPdf, std::default_delete<crossmods::CrossingPdf> > ( ::crossmods::LognormalTdm::* )( ::crossmods::invec,::crossmods::invec,double )) &::crossmods::LognormalTdm::blocker_decisions, "taus"_a, "taus_b"_a, "dt"_a)
.def("decisions", (::std::unique_ptr<crossmods::CrossingPdf, std::default_delete<crossmods::CrossingPdf> > ( ::crossmods::LognormalTdm::* )( ::crossmods::invec,double )) &::crossmods::LognormalTdm::decisions, "taus"_a, "dt"_a)
.def(py::init<crossmods::LognormalTdm const &>(), "arg0"_a);
m.def("stdnormcdf", (double (*)( double )) &::crossmods::stdnormcdf, "x"_a);
m.def("normcdf", (double (*)( double,double,double )) &::crossmods::normcdf, "x"_a, "m"_a, "v"_a);
py::class_<::crossmods::Vddm>(m, "Vddm")
.def_readwrite("dt", &::crossmods::Vddm::dt)
.def_readwrite("std", &::crossmods::Vddm::std)
.def_readwrite("damping", &::crossmods::Vddm::damping)
.def_readwrite("tau_threshold", &::crossmods::Vddm::tau_threshold)
.def_readwrite("pass_threshold", &::crossmods::Vddm::pass_threshold)
.def_readwrite("scale", &::crossmods::Vddm::scale)
.def_readwrite("act_threshold", &::crossmods::Vddm::act_threshold)
.def(py::init<double, double, double, double, double, double>(), "dt"_a, "std"_a, "damping"_a, "tau_threshold"_a, "scale"_a=1., "act_threshold"_a=1.)
.def("step", (double ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,double,double const *,double *,double )const) &::crossmods::Vddm::step, "acts"_a, "tau"_a, "prev_weights"_a, "new_weights"_a, "decision_prob"_a)
.def("step", (::std::tuple<std::vector<double, std::allocator<double> >, double> ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,double,::crossmods::invec,double )const) &::crossmods::Vddm::step, "acts"_a, "tau"_a, "prev_weights"_a, "decision_prob"_a)
.def("decisions", (double ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,double const *,double *,::size_t )const) &::crossmods::Vddm::decisions, "acts"_a, "taus"_a, "decidedpdf"_a, "dur"_a)
.def("blocker_decisions", (double ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,double const *,double const *,double *,::size_t )const) &::crossmods::Vddm::blocker_decisions, "acts"_a, "taus"_a, "taus_b"_a, "decidedpdf"_a, "dur"_a)
.def("decisions", (::crossmods::CrossingPdf & ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,::crossmods::CrossingPdf &,::crossmods::invec )const) &::crossmods::Vddm::decisions, "grid"_a, "out"_a, "taus"_a)
.def("decisions", (::std::unique_ptr<crossmods::CrossingPdf, std::default_delete<crossmods::CrossingPdf> > ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,::crossmods::invec )const) &::crossmods::Vddm::decisions, "grid"_a, "taus"_a)
.def("blocker_decisions", (::crossmods::CrossingPdf & ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,::crossmods::CrossingPdf &,::crossmods::invec,::crossmods::invec )const) &::crossmods::Vddm::blocker_decisions, "grid"_a, "out"_a, "taus"_a, "taus_b"_a)
.def("blocker_decisions", (::std::unique_ptr<crossmods::CrossingPdf, std::default_delete<crossmods::CrossingPdf> > ( ::crossmods::Vddm::* )( ::crossmods::Grid1d const &,::crossmods::invec,::crossmods::invec )const) &::crossmods::Vddm::blocker_decisions, "grid"_a, "taus"_a, "taus_b"_a)
.def(py::init<crossmods::Vddm const &>(), "arg0"_a);
}
