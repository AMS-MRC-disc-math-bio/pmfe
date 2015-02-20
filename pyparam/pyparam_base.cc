// Copyright (c) 2014 Andrew Gainer-Dewar.

#include "mfe.h"
#include "pmfe_types.h"
#include <gmpxx.h>

#include <boost/python.hpp>
#include <boost/python/long.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/args.hpp>
#include <boost/python/operators.hpp>
#include <boost/operators.hpp>

namespace py = boost::python;

struct mpz_to_python_int
{
    static PyObject* convert(mpz_class const& val) {
        return py::incref(py::object(int(val.get_si())).ptr());
    }
};

BOOST_PYTHON_MODULE(pyparam_base)
{
    // TODO: Document
    py::class_<pmfe::ScoreVector>("ScoreVector", py::init<>())
        .def(py::init<py::tuple, py::tuple, py::tuple, py::tuple, py::tuple>())
        .def_readonly("_multiloops", &pmfe::ScoreVector::multiloops)
        .def_readonly("_branches", &pmfe::ScoreVector::branches)
        .def_readonly("_unpaired", &pmfe::ScoreVector::unpaired)
        .def_readonly("_w", &pmfe::ScoreVector::w)
        .def_readonly("_energy", &pmfe::ScoreVector::energy)
        .def("_as_pairs", &pmfe::ScoreVector::as_pairs)
        .def(py::self == py::self).def(py::self != py::self) // Inherit equality and inequality
        ;

    // TODO: Document
    py::class_<pmfe::ParameterVector>("ParameterVector", py::init<>())
        .def(py::init<py::tuple, py::tuple, py::tuple, py::tuple>())
        .def_readonly("_multiloop_penalty", &pmfe::ParameterVector::multiloop_penalty)
        .def_readonly("_unpaired_penalty", &pmfe::ParameterVector::unpaired_penalty)
        .def_readonly("_branch_penalty", &pmfe::ParameterVector::branch_penalty)
        .def_readonly("_dummy_scaling", &pmfe::ParameterVector::dummy_scaling)
        .def("_as_pairs", &pmfe::ParameterVector::as_pairs)
        .def(py::self == py::self).def(py::self != py::self) // Inherit equality and inequality
        ;

    // TODO: Document
    // TODO: Support dangle_model without param_dir
    py::def("get_mfe_scores", pmfe::mfe_pywrap, "Get scores for MFE structure for a given sequence and set of parameters",
            (
             py::arg("seq_file"),
             py::arg("output_file"),
             py::arg("params"),
             py::arg("param_dir") = "/usr/local/share/pmfe/Turner99/pmfe",
             py::arg("dangle_model") = 1
             )
            );
}
