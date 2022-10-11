#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include "diff.h"
using namespace std;
using namespace Eigen;

static const int N = 9;

diff::diff(double a, double b, double fa, double fb, double (*func1)(double), double (*func2)(double), double (*func3)(double)) {

    this->a = a;
    this->b = b;
    this->fa = fa;
    this->fb = fb;
    p.resize(1, N);
    q.resize(1, N);
    r.resize(1, N);

    h = (b - a) / (N + 1);

    for (int i = 0; i < N; i++) {

        p(0, i) = func1((i + 1) * h + a);
        q(0, i) = func2((i + 1) * h + a);
        r(0, i) = func3((i + 1) * h + a);

    }

}

void diff::setA(double a) {
    this->a = a;

}
void diff::setB(double b) {
    this->b = b;
}
void diff::setFa(double fa) {
    this->fa = fa;

}
void diff::setFb(double fb) {
    this->fb = fb;

}

void diff::getP() {
    cout << p;

}
void diff::getQ() {
    cout << q;

}
void diff::getR() {
    cout << r;

}


Eigen::Matrix<double,N,N>  diff::trigonal() {

    Matrix<double, N, N> m(N, N);
    for (int i = 0; i < (N); i++) {
        for (int j = 0; j < (N); j++) {
            if (i == j) {
                m(i, i) = 2 + h * h * q(0, j);

            }
            else if (j == i + 1 && i != N - 1) {
                m(i, j) = -1 + h / 2 * p(0, i);



            }
            else if (i == j + 1 && j != N - 1) {
                m(i, j) = -1 - h / 2 * p(0, i);

            }
            else {
                m(i, j) = 0;
            }


        }

    }
    return m;
}

Eigen::Matrix<double, 1, N> diff::b_eq() {

    Matrix<double, 1, N> b1(1, N);
    b1(0, 0) = -h * h * r(0, 0) + fa * (1 + h / 2 * p(0, 0));
    b1(0, N - 1) = -h * h * r(0, N - 1) + fb * (1 - h / 2 * p(0, N - 1));

    for (int i = 1; i < (N - 1); i++) {

        b1(0, i) = -h * h * r(0, i);

    }

    return b1;
}

Eigen::Matrix<double, N + 2, 1> diff::sol() {

    Matrix<double, N, 1> b_ = b_eq().transpose();
    Matrix<double, N, N> A_ = trigonal().inverse();
    Matrix<double, N, 1> soli = A_ * b_;
    Matrix<double, N + 2, 1> Sol;
    Sol(0, 0) = fa;
    Sol(N + 1, 0) = fb;
    for (int i = 1; i < N + 1; i++) {
        Sol(i, 0) = soli(i - 1, 0);
    }
    return Sol;
}

void diff::saveSol() {
    ofstream MyFile("filename.txt");
    Matrix<double, N + 2, 1> k = sol();
    for (int i = 0; i < N + 2; i++) {
        MyFile << (h * i + a) << setw(8) << k(i, 0) << endl;

    }
}



