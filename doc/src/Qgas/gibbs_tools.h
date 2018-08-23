#pragma once

double x_sampling(const VectorXd &a, const VectorXd &h, const MatrixXd &W, double sigma_sqrd, int i);

double h_sampling(const VectorXd &b, int i);
