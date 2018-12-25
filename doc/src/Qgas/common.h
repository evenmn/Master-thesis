#pragma once
#include <eigen3/Eigen/Dense>

// Parameters defined in main
extern int     P;
extern int     D;
extern int     N;
extern int     MC;
extern int     O;
extern int     iterations;
extern int     sampling;
extern int     optimization;
extern bool    interaction;
extern bool    one_body;
extern double  sigma;
extern double  omega;
extern double  dx;
extern double  dt;
extern double  eta;
extern double  Diff;

// Parameters defined in vmc
extern int     M;
extern int     P_half;
extern int     M_half;
extern double  sigma_sqrd;

extern Eigen::MatrixXd W;
extern Eigen::VectorXd X;
extern Eigen::VectorXd v;
extern Eigen::VectorXd e;
extern Eigen::VectorXd e_p;
extern Eigen::VectorXd a;
extern Eigen::VectorXd b;
extern Eigen::VectorXd Xa;
extern Eigen::VectorXd X_new;
extern Eigen::VectorXd X_newa;
extern Eigen::VectorXd v_new;
extern Eigen::VectorXd diff;
extern Eigen::MatrixXd Dist;
extern Eigen::MatrixXd Dist_inv;
extern Eigen::MatrixXd A_up_inv;
extern Eigen::MatrixXd A_dn_inv;
extern Eigen::MatrixXd dA_up;
extern Eigen::MatrixXd dA_dn;
