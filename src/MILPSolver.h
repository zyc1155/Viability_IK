#pragma once

#include <vector>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;

/**
 * SCIP is utilized to solve the following mixed-integer linear programming (MILP) problem:
 * \f{align}{
 *   \underset{\bm{x}\in \mathbb{R}^{m},\bm{y}\in \{0,1\}^{n}}{\text{minimize}} &\ bm{c}_1^T \bm{x} + \bm{c}_2^T \bm{y} \nonumber \\
 *   \text{subject to} & \ \bm{A}_1  \bm{x} +  \bm{A}_2  \bm{y} \leq \bm{b} \\
 *   &\  \bm{x}_{lb} \leq \bm{x} \leq \bm{x}_{ub}.
 * \f}
 *
 * Here, $\bm{A}_1 \in \mathbb{R}^{l\times m}$, $\bm{A}_2 \in \mathbb{R}^{l\times n}$, $\bm{b} \in \mathbb{R}^{l}$, $\bm{c}_1 \in \mathbb{R}^{m}$, $\bm{c}_2 \in \mathbb{R}^{n}$, and $\bm{x}_{lb}, \bm{x}_{ub} \in \mathbb{R}^{m}$.
 */
class MILPSolver
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MILPSolver();
    MILPSolver(int m, int n);
    ~MILPSolver();
    void problem(int m, int n);
    void solve(SparseMatrix<double> &A1, SparseMatrix<double> &A2, VectorXd &b, VectorXd &x_lb, VectorXd &x_ub, VectorXd &c1, VectorXd &c2);
    VectorXd get_bestsol();
    double get_bestobj();
    SCIP *get_scip_ptr();

private:
    SCIP *_scip;
    size_t _n;
    size_t _m;
    size_t _l;
    std::vector<SCIP_VAR *> _vars;
    std::vector<SCIP_CONS *> _cons;
    VectorXd d_opt;
    double obj;
};
