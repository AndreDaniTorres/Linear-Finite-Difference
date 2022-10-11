#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include "diff2.h"
#include <iostream>
using namespace std;
using namespace Eigen;


// y'' = P(x) y' + Q(x) y + R(x) 
// Funciones ingresada por el usuario.
double p(double x) {
    return -2 / x;
}

double q(double x) {
    return 2 / pow(x, 2);
}

double r(double x) {
    return sin(log(x)) / pow(x, 2);
}


int main()
{
    diff num(1, 2, 9, 1, 2, &p, &q, &r);

    Matrix<double, 11, 1> a = num.sol();
    num.saveSol();
    cout << a;

    return 0;
}