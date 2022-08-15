// linear boundary-value problem

// Sea una ecuación del tipo: y'' = P(x) y' + Q(x) y + R(x) 

#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>

#include "finite-differences.h"

// using namespace std;



FiniteDifferences::FiniteDifferences(){


} 


void FiniteDifferences::setN( int value ){
    N = value; 
} 

void FiniteDifferences::setA( float value ){
    a = value; 
} 

void FiniteDifferences::setB( float value ){
    b = value; 
} 

void FiniteDifferences::setAlpha( float value ){
    alpha = value; 
} 

void FiniteDifferences::setBeta( float value ){
    beta = value; 
} 





float FiniteDifferences::p(float x){
    return (-1*2)/x;
}

float FiniteDifferences::q(float x){
    return 2/pow(x,2);
}

float FiniteDifferences::r(float x){
    return sin(log(x))/pow(x,2);
}



void FiniteDifferences::centeredDifference(std::vector<double> &vectA, std::vector<double> &vectB, std::vector<double> &vectC, std::vector<double> &vectD){
    
    float h = (b - a) / (N + 1);
    float x = a + h;

    std::cout << "x=" << x << ", h="<< h <<std::endl;
    
    // 1
    vectA.push_back( 2 + pow(h, 2)*q(x) );
    vectB.push_back( -1 + (h/2)*p(x) );
    vectC.push_back( 0 );
    vectD.push_back( -1*pow(h, 2)*r(x) + ( 1+(h/2)*p(x) )*alpha);



    
    
    // From 2 to N-1
    for (size_t i = 2; i <= N - 1; i++)
    {
        x = a + i * h;

        vectA.push_back( 2 + pow(h, 2)*q(x) );
        vectB.push_back( -1 + (h/2)*p(x) );
        vectC.push_back( -1 - (h/2)*p(x) );
        vectD.push_back( -1*pow(h, 2)*r(x) );

    }

    // N
    x = b - h;
    
    vectA.push_back( 2 + pow(h, 2)*q(x) );
    vectB.push_back( 0 );
    vectC.push_back( -1 - (h/2)*p(x)  );
    vectD.push_back( -1*pow(h, 2)*r(x) + ( 1-(h/2)*p(x) )*beta );

}

std::vector<double> FiniteDifferences::tridiagonaLinearSystem(){


    std::vector<double> vectA;
    std::vector<double> vectB;
    std::vector<double> vectC;
    std::vector<double> vectD;

    centeredDifference(vectA, vectB, vectC, vectD);


    std::vector<double> vectL;
    std::vector<double> vectU;
    std::vector<double> vectZ; 


    // 1
    vectL.push_back( vectA.at(0) );
    vectU.push_back( vectB.at(0)/vectL.at(0) );
    vectZ.push_back( vectD.at(0)/vectL.at(0) );
    

    // From 2 to N-1
    for (size_t i = 1; i < N - 1; i++)
    {
        vectL.push_back( vectA.at(i) - vectC.at(i)*vectU.at(i-1));
        vectU.push_back( vectB.at(i)/vectL.at(i));
        vectZ.push_back( (vectD.at(i) - vectC.at(i)*vectZ.at(i-1) )/vectL.at(i) );
        
    }


    // N
    vectL.push_back( vectA.at(N - 1) - vectC.at(N - 1)*vectU.at(N - 2)); // A(8), C(8), U(7)
    vectU.push_back( 0 );
    vectZ.push_back( (vectD.at(N - 1) - vectC.at(N - 1)*vectZ.at(N - 2) )/vectL.at(N - 1) );
    

    std::cout << "Definir  W" <<std::endl;

    std::vector<double> vectW(N + 2); // se agrega W(O), W(N+1)
    
    vectW.at(0) = alpha;
    vectW.at(N + 1) = beta; // W(10) 
    vectW.at(N) = vectZ.at(N - 1); // W(9) , Z(8)


    for (size_t i = N; i >= 1; i--) // va de 8 A 1
    {
        // N para Z y para U es Z(9), U(9), pero lo que busco es Z(8) 

        vectW.at(i) = vectZ.at(i - 1) - vectU.at(i - 1)*vectW.at(i + 1);  
    
    }

    
    
    return vectW;
}
