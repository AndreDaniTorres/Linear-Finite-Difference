#ifndef DIFF_H
#define DIFF_H

#include <iostream>
#include <Eigen/Dense>



class diff {

public:

    static const int N = 9;

    diff(double a, double b, double fa, double fb, double (*func1)(double), double (*func2)(double), double (*func3)(double));

    void setA(double a);
    void setB(double b);
    void setFa(double fa);
    void setFb(double fb);

    void getP();
    void getQ();
    void getR();

    Eigen::Matrix<double, N, N> trigonal();

    Eigen::Matrix<double, 1, N> b_eq();

    Eigen::Matrix<double, N + 2, 1> sol();

    void saveSol();



private:
    Eigen::MatrixXd p;
    Eigen::MatrixXd q;
    Eigen::MatrixXd r;
    double a;
    double b;
    double fa;
    double fb;
    double h;

};

#endif