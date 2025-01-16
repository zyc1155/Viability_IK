#include "viability_ik.h"
#include "scip_exception.hpp"
#include "eigen-cdd/Polyhedron.h"

VIABILITY_IK::VIABILITY_IK()
{
    initialized = false;
}

VIABILITY_IK::VIABILITY_IK(MatrixXd const &A_, VectorXd const &t_q_lim_, VectorXd const &v_lim_, VectorXd const &a_lim_, double dt_)
{
    setProblem(A_, t_q_lim_, v_lim_, a_lim_, dt_);
}

void VIABILITY_IK::setProblem(MatrixXd const &A_, VectorXd const &t_q_lim_, VectorXd const &v_lim_, VectorXd const &a_lim_, double dt_)
{
    if (A_.cols() != v_lim_.size() || A_.cols() != a_lim_.size())
        throw std::invalid_argument("A_.cols()!=v_lim_.size() || A_.cols()!=a_lim_.size()");
    if (A_.rows() != t_q_lim_.size())
        throw std::invalid_argument("A_.rows() != t_q_lim_.size()");

    m = A_.rows();
    n = A_.cols();
    non_zero_h_S = 0;

    A = A_;
    t_q_lim = t_q_lim_;
    v_lim = v_lim_;
    a_lim = a_lim_;
    dt = dt_;
    _result.resize(n);

    /********************************************************/
    // Obtaining \f$ \tilde{\bm{a}}^{\mathrm{max}} \f$.
    t_a_lim_u = VectorXd::Zero(m);
    milpsolver.problem(2 * n, n);
    for (int i = 0; i < m; i++)
        t_a_lim_u(i) = milp(A.row(i), v_lim, a_lim);

    /********************************************************/
    // Obtaining \f$ \hat{\mathcal{S}} \f$.
    VectorXd psi = VectorXd::Constant(m, 1e-5);
    Polyhedron poly;
    auto success = poly.setHrep(A, t_q_lim);
    if (!success)
        throw std::runtime_error("Vertex calculation Failed!");

    auto vrep = poly.vrep();
    cart_h_S = vrep.first.rows();

    {
        VectorXd err;
        VectorXi s;
        for (int i = 0; i < cart_h_S; ++i)
        {
            err = A * vrep.first.row(i).transpose() - t_q_lim + psi;
            s = (err.array() >= .0).cast<int>();
            h_S.push_back(s.sparseView());
            non_zero_h_S += h_S.back().nonZeros();
        }
    }
    /********************************************************/
    // Obtaining \f$ \tilde{\bm{a}}^{\mathrm{lim}} \f$.
    ca_A = EigenMatrixToCasadiDM(A);
    casadi::DM ca_a_lim = casadi::DM::zeros(n), ca_t_a_lim_u = casadi::DM::zeros(m);
    Eigen::Map<Eigen::VectorXd>(ca_a_lim.ptr(), n, 1) = a_lim;
    Eigen::Map<Eigen::VectorXd>(ca_t_a_lim_u.ptr(), m, 1) = t_a_lim_u;

    t_a_lim = qp(ca_A, ca_a_lim, ca_t_a_lim_u, h_S);
    /********************************************************/

    casadi::Dict opts;
    opts["osqp.verbose"] = false;
    opts["osqp.eps_abs"] = 1e-4;
    opts["osqp.eps_rel"] = 1e-4;
    opts["osqp.eps_prim_inf"] = 1e-5;
    opts["osqp.eps_dual_inf"] = 1e-5;

    casadi::SpDict qp({{"h", casadi::Sparsity::dense(n, n)},
                       {"a", ca_A.sparsity()}});
    qpsolver = casadi::conic("qp_solver", "osqp", qp, opts);

    initialized = true;
}

bool VIABILITY_IK::solve(MatrixXd const &J, VectorXd const &b, VectorXd const &q0, VectorXd const &dq0)
{
    if (initialized)
    {
        if (J.cols() != n || q0.size() != n || dq0.size() != n)
            throw std::invalid_argument("J.cols() != n || q0.size() != n || dq0.size() != n");

        if (J.rows() != b.size())
            throw std::invalid_argument("J.rows() != b.size()");

        MatrixXd H = J.transpose() * J;
        H += VectorXd::Constant(n, 0.01).asDiagonal();
        VectorXd g = -b.transpose() * J;

        VectorXd dq_max, dq_min;

        // VectorXd::Constant(n, 1e-6) is added to improve the numerical stability.
        dq_max = v_lim.cwiseMin(dt * a_lim + dq0) + VectorXd::Constant(n, 1e-6);
        dq_min = (-v_lim).cwiseMax(-dt * a_lim + dq0) - VectorXd::Constant(n, 1e-6);

        VectorXd t_q = t_q_lim - A * q0 - dt / 2 * A * dq0;
        dq_max_A = (-dt / 2) * t_a_lim + (2 * t_a_lim.cwiseProduct((t_a_lim * dt * dt / 8).cwiseMax(t_q))).cwiseSqrt();

        casadi::DM ca_H = casadi::DM::zeros(n, n), ca_g = casadi::DM::zeros(n), lbx = casadi::DM::zeros(n), ubx = casadi::DM::zeros(n), uba = casadi::DM::zeros(m);
        Eigen::Map<Eigen::MatrixXd>(ca_H.ptr(), n, n) = H;
        Eigen::Map<Eigen::VectorXd>(ca_g.ptr(), n, 1) = g;
        Eigen::Map<Eigen::VectorXd>(lbx.ptr(), n, 1) = dq_min;
        Eigen::Map<Eigen::VectorXd>(ubx.ptr(), n, 1) = dq_max;
        Eigen::Map<Eigen::VectorXd>(uba.ptr(), m, 1) = dq_max_A;

        casadi::DMDict args;
        args["h"] = ca_H;
        args["g"] = ca_g;
        args["lbx"] = lbx;
        args["ubx"] = ubx;
        args["uba"] = uba;
        args["a"] = ca_A;

        casadi::DMDict res = qpsolver(args);
        for (int i = 0; i < n; i++)
            _result(i) = (double)res["x"](i);

        return qpsolver.stats()["success"];
    }
    else
        return false;
}

VectorXd VIABILITY_IK::get_t_a_lim_u()
{
    return t_a_lim_u;
}

VectorXd VIABILITY_IK::get_t_a_lim()
{
    return t_a_lim;
}

VectorXd VIABILITY_IK::result()
{
    return _result;
}

double VIABILITY_IK::milp(MatrixXd const &A_i, VectorXd const &v_lim, VectorXd const &a_lim)
{
    SparseMatrix<double> A1_milp(4 * n + 1, 2 * n), A2_milp(4 * n + 1, n);
    VectorXd b_milp(4 * n + 1), c1_milp(2 * n), c2_milp(n), x_lb(2 * n), x_ub(2 * n);
    VectorXd A_iT_abs = A_i.transpose().cwiseAbs();

    VectorXd vector_M = A_iT_abs.cwiseProduct(a_lim.cwiseMax(2 * v_lim / dt - a_lim));
    double M = vector_M.maxCoeff();

    VectorXd intermediate1 = A_iT_abs.cwiseProduct(a_lim), intermediate2 = A_iT_abs.cwiseProduct(v_lim) / dt;
    for (int i = 0; i < n; i++)
    {
        // A1 and A2
        A1_milp.insert(i, i) = 1.0;
        A1_milp.insert(n + i, i) = -1.0;
        A1_milp.insert(2 * n + i, i) = 1.0;
        A1_milp.insert(3 * n + i, i) = -1.0;

        A1_milp.insert(2 * n + i, n + i) = -A_i(i) / dt;
        A1_milp.insert(3 * n + i, n + i) = A_i(i) / dt;
        A1_milp.insert(4 * n, n + i) = -A_i(i);

        A2_milp.insert(n + i, i) = M;
        A2_milp.insert(3 * n + i, i) = -M;
    }

    // b
    b_milp.head(n) = intermediate1;
    b_milp.segment(n, n) = -intermediate1 + VectorXd::Constant(n, M);
    b_milp.segment(2 * n, n) = intermediate2;
    b_milp.segment(3 * n, n) = -intermediate2;
    b_milp(4 * n) = .0;

    // c1 and c2
    c1_milp.head(n) = VectorXd::Ones(n);
    c1_milp.tail(n) = VectorXd::Zero(n);
    c2_milp = VectorXd::Zero(n);

    // x_lb and x_ub
    x_lb.head(n) = VectorXd::Constant(n, -SCIPinfinity(milpsolver.get_scip_ptr()));
    x_lb.tail(n) = -v_lim;

    x_ub.head(n) = VectorXd::Constant(n, SCIPinfinity(milpsolver.get_scip_ptr()));
    x_ub.tail(n) = v_lim;

    try
    {
        milpsolver.solve(A1_milp, A2_milp, b_milp, x_lb, x_ub, c1_milp, c2_milp);
    }
    catch (SCIPException const &exp) // catch SCIP errors
    {
        std::cerr << "SCIP Error: " << exp.what() << std::endl;
        abort();
    }
    catch (std::logic_error &exp) // catch other errors
    {
        std::cerr << "Error: " << exp.what() << std::endl;
        abort();
    }
    catch (...) // catch other errors
    {
        // do nothing, but abort in debug mode
        abort();
    }

    return milpsolver.get_bestobj();
}

VectorXd VIABILITY_IK::qp(casadi::DM const &A_, casadi::DM const &a_lim_, casadi::DM const &t_a_lim_u_, std::vector<SparseVector<int>> const &h_S_)
{
    // The original objective function is in least square form.
    // Thus, the following convention is utilized:
    // \f$ \frac{1}{2} \| \bm{R}_{qp}\bm{x}-\bm{s}_{qp} \|^2_{\bm{W}_{qp}} \propto \frac{1}{2}{\bm{x}^T\bm{H}_{qp}\bm{x}} + {\bm{g}_{qp}^T\bm{x}} \f$,
    // where \f$ \bm{H}_{qp}= \bm{R}_{qp}^T\bm{W}_{qp}bm{R}_{qp} \f$ and \f$ \bm{g}_{qp}= -\bm{R}_{qp}^T\bm{W}_{qp}\bm{s}_{qp} \f$.
    casadi::DM R_qp(m + 1 + cart_h_S * n, m + 1 + cart_h_S * n), s_qp(m + 1 + cart_h_S * n, 1), W_qp(m + 1 + cart_h_S * n, m + 1 + cart_h_S * n);
    for (int i = 0; i < m; i++)
    {
        R_qp(i, i) = 1;
        W_qp(i, i) = 1 / casadi::DM::dot(A_(i, casadi::Slice()), A_(i, casadi::Slice()));
    }
    R_qp(m, m) = 1;
    W_qp(m, m) = 1000.0;

    double c0 = 0.5;
    s_qp(casadi::Slice(0, m)) = t_a_lim_u_;
    s_qp(m) = c0;

    casadi::DM para = casadi::DM::mtimes(R_qp.T(), W_qp);
    casadi::DM H_qp = casadi::DM::mtimes(para, R_qp), g_qp = -casadi::DM::mtimes(para, s_qp);
    casadi::DM lbx = casadi::DM::vertcat({-casadi::DM::inf(m), 0, casadi::DM::repmat(-a_lim_, cart_h_S, 1)});
    casadi::DM ubx = casadi::DM::vertcat({t_a_lim_u_, c0, casadi::DM::repmat(a_lim_, cart_h_S, 1)});
    casadi::DM A_qp(m + non_zero_h_S, m + 1 + cart_h_S * n), uba = casadi::DM::zeros(m + non_zero_h_S);

    // \f$ \tilde{c} \tilde{\bm{a}}^{\mathrm{max}} \leq \tilde{\bm{a}}^{\mathrm{lim}} \f$
    A_qp(casadi::Slice(0, m), casadi::Slice(0, m)) = -1 * casadi::DM::eye(m);
    A_qp(casadi::Slice(0, m), m) = t_a_lim_u_;

    int index = m;
    // \f$ \forall \bm{s} \in \hat{\mathcal{S}} \f$
    for (int i = 0; i < cart_h_S; i++)
    {
        // Only construct the constraints corresponding to \f$ s_j=1 \f$, where \f$ j \f$ corresponds to it.index()
        for (SparseVector<int>::InnerIterator it(h_S_[i]); it; ++it)
        {
            A_qp(index, it.index()) = 1.0;                                                                    // Coefficients corresponding to \f$ \tilde{\bm{a}}^{\mathrm{lim}} \f$
            A_qp(index, casadi::Slice(m + 1 + i * n, m + 1 + (i + 1) * n)) = A_(it.index(), casadi::Slice()); // Coefficients corresponding to \f$ \hat{\ddot{\bm{q}}}_i \f$
            index++;
        }
    }

    casadi::Dict opts;
    opts["osqp.verbose"] = false;

    casadi::SpDict qp({{"h", H_qp.sparsity()},
                       {"a", A_qp.sparsity()}});

    casadi::Function solver = casadi::conic("qp_solver", "osqp", qp, opts);
    casadi::DMDict args;
    args["h"] = H_qp;
    args["g"] = g_qp;
    args["lbx"] = lbx;
    args["ubx"] = ubx;
    args["a"] = A_qp;
    args["uba"] = uba;
    args["x0"] = casadi::DM::vertcat({c0 * t_a_lim_u_, c0, casadi::DM::zeros(cart_h_S * n)});

    casadi::DMDict res = solver(args);

    if ((double)res["x"](m) < 1e-6)
        throw std::runtime_error("t_a_lim negative!");

    VectorXd t_a_lim_(m);
    for (int i = 0; i < m; i++)
        t_a_lim_(i) = (double)res["x"](i);
    return t_a_lim_;
}

casadi::DM VIABILITY_IK::EigenMatrixToCasadiDM(MatrixXd const &A_)
{
    SparseMatrix<double> sp_A = A_.sparseView(1.0, 1e-5);
    std::vector<casadi_int> row_indices;
    std::vector<casadi_int> col_indices;
    std::vector<double> values;
    for (int k = 0; k < sp_A.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(sp_A, k); it; ++it)
        {
            row_indices.push_back(it.row());
            col_indices.push_back(it.col());
            values.push_back(it.value());
        }
    }
    casadi::Sparsity SpA = casadi::Sparsity::triplet(
        sp_A.rows(),
        sp_A.cols(),
        row_indices,
        col_indices);

    casadi::DM ca_A(SpA);
    for (int k = 0; k < values.size(); k++)
    {
        ca_A(row_indices[k], col_indices[k]) = values[k];
    }

    return ca_A;
}