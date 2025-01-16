#include "MILPSolver.hpp"
#include "scip_exception.hpp"

#include <sstream>
#include <iostream>

MILPSolver::MILPSolver()
{
    _m = 0;
    _n = 0;
    _scip = 0;
}
MILPSolver::MILPSolver(int m, int n)
{
    _m = 0;
    _n = 0;
    _scip = 0;
    problem(m, n);
}

MILPSolver::~MILPSolver()
{
    try
    {
        for (SCIP_VAR *var : _vars)
            SCIP_CALL_EXC(SCIPreleaseVar(_scip, &var));
        _vars.clear();

        for (SCIP_CONS *cons : _cons)
            SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
        _cons.clear();

        // after releasing all vars and cons we can free the scip problem
        // remember this has allways to be the last call to scip
        SCIP_CALL_EXC(SCIPfree(&_scip));
    }
    catch (SCIPException const &exp) // catch SCIP errors
    {
        std::cerr << "SCIP Error: " << exp.what() << std::endl;
        abort();
    }
}
void MILPSolver::problem(int m, int n)
{
    if (_m != 0 || _n != 0)
    {
        try
        {
            for (SCIP_VAR *var : _vars)
                SCIP_CALL_EXC(SCIPreleaseVar(_scip, &var));
            _vars.clear();

            for (SCIP_CONS *cons : _cons)
                SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
            _cons.clear();

            // after releasing all vars and cons we can free the scip problem
            // remember this has allways to be the last call to scip
            SCIP_CALL_EXC(SCIPfree(&_scip));
        }
        catch (SCIPException const &exp) // catch SCIP errors
        {
            std::cerr << "SCIP Error: " << exp.what() << std::endl;
            abort();
        }
    }

    _m = m;
    _n = n;
    _vars.resize(m + n);
    d_opt.resize(m + n);

    try
    { // initialize scip
        SCIP_CALL_EXC(SCIPcreate(&_scip));

        // load default plugins linke separators, heuristics, etc.
        SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));

        // disable scip output to stdout
        SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(_scip), TRUE);

        // create an empty problem
        SCIP_CALL_EXC(SCIPcreateProb(_scip, "MILPSolver", NULL, NULL, NULL, NULL, NULL, NULL, NULL));

        std::ostringstream namebuf;
        // create continous variables
        for (int i = 0; i < _m; i++)
        {
            SCIP_VAR *var = NULL;
            namebuf.str("");
            namebuf << "x#" << i + 1;

            // create the SCIP_VAR object
            SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), -SCIPinfinity(_scip), SCIPinfinity(_scip), 1.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));

            // add the SCIP_VAR object to the scip problem
            SCIP_CALL_EXC(SCIPaddVar(_scip, var));

            // storing the SCIP_VAR pointer for later access
            _vars[i] = var;
        }

        // create binary variables
        for (int i = 0; i < _n; i++)
        {
            SCIP_VAR *var = NULL;
            namebuf.str("");
            namebuf << "y#" << i + 1;

            // create the SCIP_VAR object
            SCIP_CALL_EXC(SCIPcreateVar(_scip, &var, namebuf.str().c_str(), 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));

            // add the SCIP_VAR object to the scip problem
            SCIP_CALL_EXC(SCIPaddVar(_scip, var));

            // storing the SCIP_VAR pointer for later access
            _vars[_m + i] = var;
        }
    }
    catch (SCIPException const &exp) // catch SCIP errors
    {
        std::cerr << "SCIP Error: " << exp.what() << std::endl;
        abort();
    }
}
void MILPSolver::solve(SparseMatrix<double> &A1, SparseMatrix<double> &A2, VectorXd &b, VectorXd &x_lb, VectorXd &x_ub, VectorXd &c1, VectorXd &c2)
{
    if (A1.cols() != _m || c1.size() != _m || x_lb.size() != _m || x_ub.size() != _m)
        throw std::invalid_argument("A1.cols() != m || c1.size() != m || x_lb.size() != m || x_ub.size() != m");

    if (A2.cols() != _n || c2.size() != _n)
        throw std::invalid_argument("A2.cols() != n || c2.size() != n");

    if (A1.rows() != b.size() || A2.rows() != b.size())
        throw std::invalid_argument("A1.rows() != b.size() || A2.rows() != b.size()");

    // delete old constraints
    if (!_cons.empty())
    {
        // reset SCIP's state
        SCIP_CALL_EXC(SCIPfreeTransform(_scip));
        for (SCIP_CONS *cons : _cons)
        {
            SCIP_CALL_EXC(SCIPdelCons(_scip, cons));
            // SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
        }
        _cons.clear();
    }

    // reset the variable's objective value according to $\bm{c}_1$ and $\bm{c}_2$
    for (size_t i = 0; i < _m; i++)
        SCIP_CALL_EXC(SCIPchgVarObj(_scip, _vars[i], c1(i)));
    for (size_t i = 0; i < _n; i++)
        SCIP_CALL_EXC(SCIPchgVarObj(_scip, _vars[_m + i], c2(i)));

    // reset the continous variable's bound according to $\bm{x}_{lb}$ and $bm{x}_{ub}$
    for (size_t i = 0; i < _m; i++)
    {
        SCIP_CALL_EXC(SCIPchgVarLb(_scip, _vars[i], x_lb(i)));
        SCIP_CALL_EXC(SCIPchgVarUb(_scip, _vars[i], x_ub(i)));
    }

    // add new constraints
    _l = b.size();
    std::ostringstream namebuf;
    for (size_t i = 0; i < _l; i++)
    {
        SCIP_CONS *cons = NULL;
        namebuf.str("");
        namebuf << "cons#" << i + 1;

        // create SCIP_CONS object
        SCIP_CALL_EXC(SCIPcreateConsLinear(_scip, &cons, namebuf.str().c_str(), 0, NULL, NULL, -SCIPinfinity(_scip), b(i),
                                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

        // store the constraint for later on
        _cons.push_back(cons);
    }
    // Iterate over the non-zero elements of the sparse matrices
    for (int i = 0; i < A1.outerSize(); ++i)
    {
        for (SparseMatrix<double>::InnerIterator it(A1, i); it; ++it)
            SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, _cons[it.row()], _vars[it.col()], it.value()));
    }

    for (int i = 0; i < A2.outerSize(); ++i)
    {
        for (SparseMatrix<double>::InnerIterator it(A2, i); it; ++it)
            SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, _cons[it.row()], _vars[_m + it.col()], it.value()));
    }

    // add constraint to the problem
    for (size_t i = 0; i < _l; i++)
        SCIP_CALL_EXC(SCIPaddCons(_scip, _cons[i]));

    // this tells scip to start the solution process
    SCIP_CALL_EXC(SCIPsolve(_scip));

    // get the best found solution from scip
    SCIP_SOL *sol = SCIPgetBestSol(_scip);
    if (sol != NULL)
    {
        obj = SCIPgetPrimalbound(_scip);
        for (size_t i = 0; i < _m + _n; i++)
        {
            d_opt(i) = SCIPgetSolVal(_scip, sol, _vars[i]);
        }
    }
    else
        throw std::domain_error("No solution found!");
}
VectorXd MILPSolver::get_bestsol()
{
    return d_opt;
}

double MILPSolver::get_bestobj()
{
    return obj;
}

SCIP *MILPSolver::get_scip_ptr()
{
    return _scip;
}
