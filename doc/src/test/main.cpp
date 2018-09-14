#include <iostream>
#include <ctime>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/LU"

using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd a = MatrixXd::Zero(2,2);
    MatrixXd b = MatrixXd::Zero(2,2);

    a << 1,2,
         3,4;

    b << 5,6,
         7,8;


    VectorXd e = VectorXd::Zero(2,1);

    e(0) = 1;
    e(1) = 2;

    //cout << e.cwiseAbs2() << endl;

    MatrixXd c = MatrixXd::Random(28,28);

    clock_t start_time = clock();
    cout << c.determinant() << endl;
    clock_t end_time = clock();

    double CPU_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    cout << "CPU-time: " << CPU_time << endl;

    clock_t start_time2 = clock();
    cout << c.inverse() << endl;
    clock_t end_time2 = clock();

    double CPU_time2 = (double)(end_time2 - start_time2)/CLOCKS_PER_SEC;
    cout << "CPU-time: " << CPU_time2 << endl;

    return 0;
}
