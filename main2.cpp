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
    return -1;
}

double q(double x) {
    return -100;
}

double r(double x) {
    return 0;
}

        
int main()
{
    diff num(0, 15, 100, 0, 3, &p, &q, &r);

    Matrix<double, Dynamic, Dynamic> a = num.sol();
    num.saveSol();
   // num.print();
    //cout << a;

    return 0;
}