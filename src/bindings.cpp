#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "viability_ik.h"

namespace py = pybind11;

PYBIND11_MODULE(pyviability_ik, m)
{
     py::class_<VIABILITY_IK>(m, "VIABILITY_IK")
         .def(py::init<>())
         .def(py::init<MatrixXd const &, VectorXd const &, VectorXd const &, VectorXd const &, double>(),
              py::arg("A"), py::arg("t_q_lim"), py::arg("v_lim"), py::arg("a_lim"), py::arg("dt"))
         .def("setProblem", &VIABILITY_IK::setProblem,
              py::arg("A"), py::arg("t_q_lim"), py::arg("v_lim"), py::arg("a_lim"), py::arg("dt"))
         .def("solve", &VIABILITY_IK::solve,
              py::arg("J"), py::arg("b"), py::arg("q0"), py::arg("dq0"))
         .def("result", &VIABILITY_IK::result)
         .def("get_t_a_lim_u", &VIABILITY_IK::get_t_a_lim_u)
         .def("get_t_a_lim", &VIABILITY_IK::get_t_a_lim);
}