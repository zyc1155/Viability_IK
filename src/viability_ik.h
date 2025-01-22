#pragma once

#include <Eigen/Eigen>
#include "casadi/casadi.hpp"
#include "MILPSolver.h"

using namespace Eigen;

class VIABILITY_IK
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VIABILITY_IK();
    VIABILITY_IK(MatrixXd const &A, VectorXd const &t_q_lim, VectorXd const &v_lim, VectorXd const &a_lim, double dt);

    // Algorithm 1
    void setProblem(MatrixXd const &A, VectorXd const &t_q_lim, VectorXd const &v_lim, VectorXd const &a_lim, double dt);

    // IK, used in every loop of online computation stage
    bool solve(MatrixXd const &J, VectorXd const &b, VectorXd const &q0, VectorXd const &dq0);

    // Get \f$ \tilde{\bm{a}}^{\mathrm{lim}}_u \f$
    VectorXd get_t_a_lim_u();

    // Get \f$ \tilde{\bm{a}}^{\mathrm{lim}} \f$
    VectorXd get_t_a_lim();

    // Get result
    VectorXd result();

private:
    /**
     * SCIP is utilized to solve the following mixed-integer linear programming (MILP) problem:
     * \f{align}{
     *   \underset{\bm{x}\in \mathbb{R}^{m},\bm{y}\in \{0,1\}^{n}}{\text{minimize}} &\ bm{c}_1^T \bm{x} + \bm{c}_2^T \bm{y} \nonumber \\
     *   \text{subject to} & \ \bm{A}_1  \bm{x} +  \bm{A}_2  \bm{y} \leq \bm{b} \\
     *   &\  \bm{x}_{lb} \leq \bm{x} \leq \bm{x}_{ub}.
     * \f}
     *
     */
    double milp(MatrixXd const &A_i, VectorXd const &v_lim, VectorXd const &a_lim);

    /**
     * The low-level Quadratic programming (QP) interface (casadi::conic) of CasADi is utilized, which solves the follows:
     * \f{align}{
     *   \underset{\bm{x}}{\text{minimize}} & \ \frac{1}{2}{\bm{x}^T\bm{H}_{qp}\bm{x}} + {\bm{g}_{qp}^T\bm{x}} \nonumber \\
     *   \text{subject to} & \ \bm{x}_{lb} \leq \bm{x} \leq \bm{x}_{ub} \\
     *   & \ \bm{a}_{lb} \leq \bm{A}_{qp}\bm{x} \leq \bm{a}_{ub}
     * \f}
     *
     * OSQP is utilized as the QP solver.
     */
    VectorXd qp(casadi::DM const &A, casadi::DM const &a_lim, casadi::DM const &t_a_lim_u, std::vector<SparseVector<int>> const &S);

    casadi::DM EigenMatrixToCasadiDM(MatrixXd const &A);

    bool initialized;
    int m;          // \f$ m \f$
    int n;          // \f$ n \f$
    int cart_S;     // Cardinality of \f$ \mathcal{S} \f$
    int non_zero_S; // Total number of non-zeros of all elements belonging to \f$ \mathcal{S} \f$
    double dt;      // Sampling time \f$ T \f$

    MatrixXd A;                       // \f$ \bm{A} \f$
    VectorXd t_q_lim;                 // \f$ \tilde{\bm{q}}^{\mathrm{lim}} \f$
    VectorXd v_lim;                   // \f$ \bm{v}}^{\mathrm{lim}} \f$
    VectorXd a_lim;                   // \f$ \bm{a}}^{\mathrm{lim}} \f$
    VectorXd t_a_lim_u;               // \f$ \tilde{\bm{a}}^{\mathrm{lim}}_u \f$
    VectorXd t_a_lim;                 // \f$ \tilde{\bm{a}}^{\mathrm{lim}} \f$
    VectorXd dq_max_A;                // \f$ \overline{\bm{q}}^{k-1}_A \f$
    VectorXd dq_max, dq_min;          // \f$ \overline{\bm{q}}^{k-1}, \underline{\bm{q}}^{k-1} \f$
    std::vector<SparseVector<int>> S; // \f$ \mathcal{S} \f$

    VectorXd _result;
    casadi::DM ca_A; // \f$ \bm{A} \f$ stored in casadi form
    casadi::Function qpsolver;
    MILPSolver milpsolver;
};