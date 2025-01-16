#include <iostream>
#include <iomanip>
#include <Eigen/Core>

#include "viability_ik.h"
using namespace Eigen;

MatrixXd Jacobian_matrix(VectorXd const &q)
{
    return MatrixXd::Identity(q.size(), q.size());
}

VectorXd Forward_kinematics(VectorXd const &q)
{
    return q;
}
void print_1_2_row(VectorXd const &q0, VectorXd const &dq0, VectorXd const &ddq0);
void print_i_row(double time, VectorXd const &q1, VectorXd const &dq1, VectorXd const &ddq1);

int main()
{

    int m = 8, n = 3;
    double dt = 0.05, simulation_time = .0;
    VIABILITY_IK vik;
    MatrixXd A(m, n);
    VectorXd t_q_lim(m);
    VectorXd v_lim(n), a_lim(n);

    bool success;
    double tau = 0.2;
    VectorXd p_d(n), q0(n), dq0(n);
    VectorXd q1(n), dq1(n), ddq1(n);

    A << 0, 0, -1,
        -1, 4.5, 1.65,
        1, 5.2, 3.2,
        -1.625, -2.25, 1,
        -6, -1, 1.5,
        -1, -34.6667, 6.66667,
        2.82143, 1, 1.57143,
        1, -1, -6.66667;
    t_q_lim << 0,
        5.6,
        9.92,
        0,
        -0,
        0,
        6.11429,
        0;
    v_lim << 1, 1.5, 1.2;
    a_lim << 10.0, 15.0, 12.0;

    p_d << -0.4, -0.6, 1.3;
    q0 << .0, .0, 0.1;
    dq0 = VectorXd::Zero(n);

    vik.setProblem(A, t_q_lim, v_lim, a_lim, dt);

    VectorXd t_a_lim_u = vik.get_t_a_lim_u();
    VectorXd t_a_lim = vik.get_t_a_lim();

    std::cout << "Given constraint parameters:" << std::endl;
    std::cout << "A:" << "\n"
              << A << std::endl;
    std::cout << "t_q_lim:" << "\n"
              << t_q_lim.transpose() << std::endl;
    std::cout << "v_lim:" << "\n"
              << v_lim.transpose() << std::endl;
    std::cout << "a_lim:" << "\n"
              << a_lim.transpose() << std::endl;
    std::cout << std::endl;

    std::cout << "parameters obtained from Algorithm 1:" << std::endl;
    std::cout << "t_a_lim_u:" << "\n"
              << t_a_lim_u.transpose() << std::endl;
    std::cout << "t_a_lim:" << "\n"
              << t_a_lim.transpose() << std::endl;
    std::cout << std::endl;

    print_1_2_row(q0, dq0, VectorXd::Zero(3));

    do
    {
        success = vik.solve(Jacobian_matrix(q0), (p_d - Forward_kinematics(q0)) / tau, q0, dq0);
        if (success)
        {
            dq1 = vik.result();
            q1 = q0 + dt * (dq0 + dq1) / 2.0;
            ddq1 = (dq1 - dq0) / dt;

            simulation_time += dt;
            if (simulation_time > 100)
            {
                std::cout << "Time up!" << std::endl;
                break;
            }
            print_i_row(simulation_time, q1, dq1, ddq1);

            q0 = q1;
            dq0 = dq1;
        }
    } while (success && dq1.norm() > 1e-3 && ddq1.norm() > 1e-3);

    if (success)
        std::cout << "Finished!" << std::endl;
    else
        std::cout << "Failed!" << std::endl;
    return 0;
}

// Function to output a number with color
void printColored(double value, bool highlight, int precision)
{
    if (highlight)
    {
        // ANSI escape code for red text
        std::cout << "\033[1;31m"; // Bold Red
    }
    std::cout << std::fixed << std::setprecision(precision) << value;
    if (highlight)
    {
        // Reset color
        std::cout << "\033[0m";
    }
}

void print_1_2_row(VectorXd const &q0, VectorXd const &dq0, VectorXd const &ddq0)
{
    std::cout << "IK start:" << std::endl;
    std::cout << "time";
    for (int i = 0; i < 3; i++)
        std::cout << "\t" << "q_" + std::to_string(i);
    for (int i = 0; i < 3; i++)
        std::cout << "\t" << "dq_" + std::to_string(i);
    for (int i = 0; i < 3; i++)
        std::cout << "\t" << "ddq_" + std::to_string(i);
    std::cout << std::endl;

    std::cout << std::fixed << std::setprecision(2) << .0;
    for (int i = 0; i < 3; i++)
    {
        std::cout << "\t";
        printColored(q0(i), false, 3);
    }
    for (int i = 0; i < 3; i++)
    {
        std::cout << "\t";
        printColored(dq0(i), false, 3);
    }
    for (int i = 0; i < 3; i++)
    {
        std::cout << "\t";
        printColored(ddq0(i), false, 3);
    }
    std::cout << std::endl;
}
void print_i_row(double time, VectorXd const &q1, VectorXd const &dq1, VectorXd const &ddq1)
{
    std::cout << std::fixed << std::setprecision(2) << time;
    for (int i = 0; i < 3; i++)
    {
        std::cout << "\t";
        printColored(q1(i), false, 3);
    }
    for (int i = 0; i < 3; i++)
    {
        std::cout << "\t";
        printColored(dq1(i), false, 3);
    }
    for (int i = 0; i < 3; i++)
    {
        std::cout << "\t";
        printColored(ddq1(i), false, 3);
    }
    std::cout << std::endl;
}