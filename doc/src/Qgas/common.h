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
extern double  steplength;
extern double  timestep;
extern double  eta;
extern double  Diff;

// Parameters defined in vmc
extern int     M;
extern int     P_half;
extern double  sigma_sqrd;

extern Eigen::MatrixXd W;
extern Eigen::VectorXd X;
extern Eigen::VectorXd v;
extern Eigen::VectorXd e;
extern Eigen::VectorXd a;
extern Eigen::VectorXd b;
extern Eigen::VectorXd Xa;
