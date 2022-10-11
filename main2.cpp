#include "diff.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include "diff.h"
using namespace std;
using namespace Eigen;


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

    diff num(1, 2, 1, 2, &p, &q, &r);

    Matrix<double, 11, 1> a = num.sol();
    num.saveSol();
    cout << a;

    return 0;



}