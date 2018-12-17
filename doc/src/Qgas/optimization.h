#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Optimization
{
public:
    Optimization() {}
    void Gradient_a(const VectorXd &a, VectorXd &da);
    void Gradient_b(const VectorXd &e, VectorXd &db);
    void Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW);
    void Total_Gradient_a(const double EL_avg, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &DA);
    void Total_Gradient_b(const double EL_avg, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &DB);
    void Total_Gradient_W(const double EL_avg, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &DW);
    void GD_a(const double EL_avg, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a);
    void GD_b(const double EL_avg, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b);
    void GD_W(const double EL_avg, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W);
    void ADAM_a(int i, VectorXd &m, VectorXd &v, const double b1, const double b2, const double EL_avg, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a);
    void ADAM_b(int i, VectorXd &m, VectorXd &v, const double b1, const double b2, const double EL_avg, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b);
    void ADAM_W(int i, MatrixXd &m, MatrixXd &v, const double b1, const double b2, const double EL_avg, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W);
};
