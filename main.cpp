#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include <cassert>
#define assertm(exp, msg) assert(((void)msg, exp))

#include "tridiag.hpp"
#include "binomial.hpp"
#include "binomial1d.hpp"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(optionpricing, mod)
{
    mod.doc() = "";

    using BinomialMeshDouble = BinomialMesh<double>;

    py::class_<BinomialMeshDouble>(mod, "BinomialMesh")
        .def(py::init<int>())
        .def(py::init<int, double>())
        .def("set_initial_condition", &BinomialMeshDouble::set_initial_condition)
        .def("calculate_serial", &BinomialMeshDouble::calculate_serial)
        .def("calculate_parallel", &BinomialMeshDouble::calculate_parallel)
        .def("get", &BinomialMeshDouble::get)
        .def("print", &BinomialMeshDouble::print);

    using BinomialMeshDouble1D = BinomialMesh1D<double>;

    py::class_<BinomialMeshDouble1D>(mod, "BinomialMesh1D")
        .def(py::init<int>())
        .def(py::init<int, double>())
        .def("set_initial_condition", &BinomialMeshDouble1D::set_initial_condition)
        .def("calculate_serial", &BinomialMeshDouble1D::calculate_serial)
        .def("calculate_parallel", &BinomialMeshDouble1D::calculate_parallel);

    py::class_<BinomialPricing>(mod, "BinomialPricing")
        .def(py::init<int>())
        .def("sum_serial", &BinomialPricing::sum_serial)
        .def("sum_parallel", &BinomialPricing::sum_parallel);

#ifdef VERSION_INFO
    mod.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    mod.attr("__version__") = "dev";
#endif
}