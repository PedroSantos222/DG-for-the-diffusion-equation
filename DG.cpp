#include <iostream>
#define _USE_MATH_DEFINES // Must be before #include <cmath>
#include <Eigen/Dense>
#include<cmath>
//DG method for diffusion equation, we are going to be approximating with 4 polinomials

constexpr double pi = 3.14159265358979323846;

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd AK(4, 4);
MatrixXd AL(4, 4);
MatrixXd AR(4, 4);
MatrixXd AB(4, 4);
MatrixXd AT(4, 4);


const int Nx{ 20 };
const int Ny{ 20 };
double h = 1.0 / Nx;
int N{ Nx * Ny * 4 };

int cellIndex(int i, int j, int Ny) {
    return 4 * (i * Ny + j);
}

double phi1(double x, double y,
    double xa, double xb,
    double ya, double yb)
{
    return ((xb - x) * (yb - y))
        / ((xb - xa) * (yb - ya));
}

double phi2(double x, double y,
    double xa, double xb,
    double ya, double yb)
{
    return ((x - xa) * (yb - y))
        / ((xb - xa) * (yb - ya));
}

double phi3(double x, double y,
    double xa, double xb,
    double ya, double yb)
{
    return ((x - xa) * (y - ya))
        / ((xb - xa) * (yb - ya));
}

double phi4(double x, double y,
    double xa, double xb,
    double ya, double yb)
{
    return ((xb - x) * (y - ya))
        / ((xb - xa) * (yb - ya));
}

Eigen::Vector4d getCoefficients(const Eigen::VectorXd& solution,int i,int j,int Ny){
    int index = 4 * (i * Ny + j);
    return solution.segment<4>(index);
}

double UhCell(int i, int j, double x, double y, const Eigen::VectorXd& solution, int Ny, double h)
{// Get the four coefficients for this cell
    int index = 4 * (i * Ny + j);
    double c1 = solution(index);
    double c2 = solution(index + 1);
    double c3 = solution(index + 2);
    double c4 = solution(index + 3);

    // Coordinates of the cell
    double xa = i * h;
    double xb = (i + 1) * h;

    double ya = j * h;
    double yb = (j + 1) * h;

    // DG approximation
    return c1 * phi1(x, y, xa, xb, ya, yb)+ c2 * phi2(x, y, xa, xb, ya, yb)+ c3 * phi3(x, y, xa, xb, ya, yb)+ c4 * phi4(x, y, xa, xb, ya, yb);
}

double UhSolution(double x,double y,const Eigen::VectorXd& solution,int Nx,int Ny,double h)
{
    int i = static_cast<int>(std::floor(x / h));
    int j = static_cast<int>(std::floor(y / h));

    // Optional safety check
    if (i < 0 || i >= Nx || j < 0 || j >= Ny)
    {
        throw std::out_of_range("Point is outside the domain.");
    }

    return UhCell(i, j, x, y, solution, Ny, h);
}

double u0(double x, double y) //Initial condition for when we consider the time variable
{
    return std::sin(pi * x) * std::sin(pi * y);
}

Eigen::Vector4d initialCoefficients(int i,int j,double h){
    double xa = i * h;
    double xb = (i + 1) * h;

    double ya = j * h;
    double yb = (j + 1) * h;

    const double g = 1.0 / std::sqrt(3.0);

    double points[2] = { -g, g };

    Eigen::Vector4d result = Eigen::Vector4d::Zero();

    for (int a = 0; a < 2; ++a)
    {
        for (int b = 0; b < 2; ++b)
        {
            double xi = points[a];
            double eta = points[b];

            double x = 0.5 * (xa + xb) + 0.5 * h * xi;
            double y = 0.5 * (ya + yb) + 0.5 * h * eta;

            double weight = (h * h) / 4.0;

            result(0) += weight * u0(x, y)
                * phi1(x, y, xa, xb, ya, yb);

            result(1) += weight * u0(x, y)
                * phi2(x, y, xa, xb, ya, yb);

            result(2) += weight * u0(x, y)
                * phi3(x, y, xa, xb, ya, yb);

            result(3) += weight * u0(x, y)
                * phi4(x, y, xa, xb, ya, yb);
        }
    }

    return result;
}

int main()
{
    AK << 84, 19, -2, 19,
        19, 84, 19, -2,
        -2, 19, 84, 19,
        19, -2, 19, 84;

    AK /= 12.0;

    AL << -2, -38, -19, -1,
        0, 0, 0, 0,
        0, 0, 0, 0,
        -1, -19, -38, -2;

    AL /= 12.0;

    AR << 0, 0, 0, 0,
        -38, -2, -1, -19,
        -19, -1, -2, -38,
        0, 0, 0, 0;

    AR /= 12.0;

    AB << -2, -1, -19, -38,
        -1, -2, -38, -19,
        0, 0, 0, 0,
        0, 0, 0, 0;

    AB /= 12.0;

    AT << 0, 0, 0, 0,
        0, 0, 0, 0,
        -19, -38, -2, -1,
        -38, -19, -1, -2;

    AT /= 12.0;

    int N = 4 * Nx * Ny;

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);

    Eigen::Vector4d leftBC = Eigen::Vector4d::Constant(1.0);
    Eigen::Vector4d rightBC = Eigen::Vector4d::Constant(2.0);
    Eigen::Vector4d bottomBC = Eigen::Vector4d::Constant(3.0);
    Eigen::Vector4d topBC = Eigen::Vector4d::Constant(3.0);

    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            int row = cellIndex(i, j, Ny);

            // Current cell
            A.block<4, 4>(row, row) += AK;

            // Left
            if (i > 0)
            {
                int col = cellIndex(i - 1, j, Ny);
                A.block<4, 4>(row, col) += AL;
            }
            else
            {
                b.segment<4>(row) -= AL * leftBC;
            }

            // Right
            if (i < Nx - 1)
            {
                int col = cellIndex(i + 1, j, Ny);
                A.block<4, 4>(row, col) += AR;
            }
            else
            {
                b.segment<4>(row) -= AR * rightBC;
            }

            // Bottom
            if (j > 0)
            {
                int col = cellIndex(i, j - 1, Ny);
                A.block<4, 4>(row, col) += AB;
            }
            else
            {
                b.segment<4>(row) -= AB * bottomBC;
            }

            // Top
            if (j < Ny - 1)
            {
                int col = cellIndex(i, j + 1, Ny);
                A.block<4, 4>(row, col) += AT;
            }
            else
            {
                b.segment<4>(row) -= AT * topBC;
            }
        }
    }
    Eigen::VectorXd solution = A.fullPivLu().solve(-b);

    double value0 = UhSolution(0.25,0.37,solution,Nx,Ny,h);
    std::cout << "The value for the time-independent solution at (0.25,0.37) is " << value0 << std::endl;

    // Now we are interested in adding the time variable 

    Eigen::Matrix4d Mlocal;
    Mlocal << 4, 2, 1, 2,
        2, 4, 2, 1,
        1, 2, 4, 2,
        2, 1, 2, 4;

    Mlocal /= 324.0;

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);

    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            int index = cellIndex(i, j, Ny);

            M.block<4, 4>(index, index) = Mlocal;
        }
    }

    Eigen::MatrixXd B = -M.fullPivLu().solve(A);

    Eigen::VectorXd U0 = Eigen::VectorXd::Zero(N);

    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            int index = cellIndex(i, j, Ny);

            // m_k = integral of u0 * phi_k
            Eigen::Vector4d m =
                initialCoefficients(i, j, h);

            // Mlocal * c = m
            // Therefore c = Mlocal^(-1) * m
            Eigen::Vector4d coefficients =
                Mlocal.fullPivLu().solve(m);

            U0.segment<4>(index) = coefficients;
        }
    }


    int Nt = 2000;
    double dt = 1.0 / Nt;

    std::vector<Eigen::VectorXd> sol;

    sol.push_back(U0);

    for (int n = 0; n < Nt; ++n)
    {
        Eigen::VectorXd Unew =
            sol.back() + dt * B * sol.back();

        sol.push_back(Unew);
    }

    double value1 = UhSolution(0.25,0.37,sol[6],Nx, Ny, h);
    std::cout << "The value of the time-dependent solution at (t=0.5,0.25,0.37) is " << value1 << std::endl;
}
