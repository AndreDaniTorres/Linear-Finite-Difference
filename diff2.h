#ifndef DIFF_H
#define DIFF_H

#include <Eigen/Dense>


//class

class diff {

public:

    
    //Constructor
    diff(double , double , int, double , double , double (*)(double), double (*)(double), double (*)(double)); //(*func2)(double) //a,b,N,fa,fb,P,Q,R

    // Se fijan los puntos finales y las condiciones de frontera.
    void setA(double );
    void setB(double );
    void setFa(double );
    void setFb(double );

    void setP(double (*)(double),double , double,int);
    void setQ(double (*)(double), double, double,int);
    void setR(double (*)(double), double, double,int);


    //Muestra las funciones 
    void getP();  
    void getQ();
    void getR();

    // Solucion a Ax=b

    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> trigonal(); //A
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_eq(); //b
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sol(); //x

    void saveSol();



private:
    //funciones
    Eigen::MatrixXd p;
    Eigen::MatrixXd q;
    Eigen::MatrixXd r;

    //Puntos finales, condiciones iniciales y paso.
    double a;
    double b;
    double fa;
    double fb;
    double h;
    int N;

};

#endif