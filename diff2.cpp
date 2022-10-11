// linear boundary-value problem
// Sea una ecuación del tipo: y'' = P(x) y' + Q(x) y + R(x) 

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include <iostream> 
#include "diff2.h"
using namespace std;
using namespace Eigen;

// Pasos


// --------------------------------
// --------- CONSTRUCTOR ----------
// --------------------------------

diff::diff(double a, double b,int N, double fa, double fb, double (*func1)(double), double (*func2)(double), double (*func3)(double)) {

    this->a = a;
    this->b = b;
    this->N = N;
    this->fa = fa;
    this->fb = fb;
    this->setP(func1, a, b, N);
    this->setQ(func2, a, b, N);
    this->setR(func3, a, b, N);

    

}

// --------------------------------
// --------- GET ----------
// --------------------------------
void diff::getP() {
    cout << p<<endl;

}
void diff::getQ() {
    cout << q << endl;

}
void diff::getR() {
    cout << r << endl;

}

// --------------------------------
// --------- SET ----------
// --------------------------------

void diff::setA(double a) {
    this->a = a;
}
void diff::setB(double b) {
    this->b = b;
}
void diff::setFa(double fa)
{
    this->fa = fa;
}
void diff::setFb(double fb) {
    this->fb = fb;
}

void diff::setP(double (*func1)(double), double a, double b,int N) {

    p.resize(1, N);
    h = (b - a) / (N + 1);

    for (int i = 0; i < N; i++) {

        p(0, i) = func1((i + 1) * h + a);

    }
}
void diff::setQ(double (*func1)(double), double a, double b, int N) {
    q.resize(1, N);
    h = (b - a) / (N + 1);
 
    for (int i = 0; i < N; i++) {

        q(0, i) = func1((i + 1) * h + a);

    }
}

void diff::setR(double (*func1)(double), double a, double b, int N) {
    r.resize(1, N);
    h = (b - a) / (N + 1);
   

    for (int i = 0; i < N; i++) {

        r(0, i) = func1((i + 1) * h + a);
    }
}


// --------------------------------
// --------- Matriz tridiagonal ----------
// --------------------------------
Eigen::Matrix<double, Dynamic, Dynamic>  diff::trigonal() {

    MatrixXd  m(N, N);
    for (int i = 0; i < (N); i++) {
        for (int j = 0; j < (N); j++) {
            if (i == j) { //Diagonal 
                m(i, i) = 2 + h * h * q(0, j);
            }
            else if (j == i + 1 && i != N - 1) { //Diagonal superior
                m(i, j) = -1 + h / 2 * p(0, i);
            }
            else if (i == j + 1 && j != N - 1) { //Diagonal inferior
                m(i, j) = -1 - h / 2 * p(0, i);
            }
            else {
                m(i, j) = 0;
            }
        }
    }
    return m;
}

// --------------------------------
// --------- b ----------
// --------------------------------
Eigen::Matrix<double, Dynamic, Dynamic> diff::b_eq() {
    MatrixXd  b1(1, N);
    b1(0, 0) = -h * h * r(0, 0) + fa * (1 + h / 2 * p(0, 0));
    b1(0, N - 1) = -h * h * r(0, N - 1) + fb * (1 - h / 2 * p(0, N - 1));

    for (int i = 1; i < (N - 1); i++) {
        b1(0, i) = -h * h * r(0, i);
    }

    return b1;
}

// --------------------------------
// --------- Soluciones ----------
// --------------------------------


Eigen::Matrix<double, Dynamic, Dynamic> diff::sol() {

    MatrixXd  b_ = b_eq().transpose();
    MatrixXd A_ = trigonal().inverse();
    MatrixXd soli = A_ * b_;
    Matrix<double, Dynamic, Dynamic> Sol;
    Sol.resize(N + 2, 1);
    Sol(0, 0) = fa;
    Sol(N + 1, 0) = fb;
    for (int i = 1; i < N + 1; i++) {
        Sol(i, 0) = soli(i - 1, 0);
    }
    return Sol;
}

// --------------------------------
// --------- Save coordinate and solutions ----------
// --------------------------------

void diff::saveSol() {
    //Create file

    ofstream MyFile("filename.txt");
    MatrixXd k = sol();
    for (int i = 0; i < N + 2; i++) {
        MyFile << (h * i + a) << setw(8) << k(i, 0) << endl;

    }
}



