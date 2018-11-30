#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Optimization
{
private:
    int m_sampling;
    double m_sigma_sqrd;
    int m_M;
    int m_N;
public:
    Optimization() {}
    int init(int sampling, double sigma_sqrd, int M, int N);
    void Gradient_a(const VectorXd &a, VectorXd &da);
    void Gradient_b(const VectorXd &e, VectorXd &db);
    void Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW);
    void Total_Gradient_a(const double EL_avg, const double MC, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &DA);
    void Total_Gradient_b(const double EL_avg, const double MC, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &DB);
    void Total_Gradient_W(const double EL_avg, const double MC, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &DW);
    void GD_a(const double eta, const double EL_avg, const double MC, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a);
    void GD_b(const double eta, const double EL_avg, const double MC, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b);
    void GD_W(const double eta, const double EL_avg, const double MC, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W);
};
