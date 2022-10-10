#include <iostream>
#include <vector>


    
// Class
class FiniteDifferences
{
    public:
    FiniteDifferences(); 

    // Set Boundary conditions and endpoints
    void setN(/*int const&*/);
    void setA();
    void setB();
    void setAlpha();
    void setBeta();

    // Functions boundary-value problem
    float p(float x);
    float q(float x);
    float r(float x);

    // Linear Finite-Difference solution
    std::vector<double> tridiagonaLinearSystem();
    
    private:

    // Centered-difference formula
    void centeredDifference(std::vector<double> &, std::vector<double> &, 
                            std::vector<double> &, std::vector<double> &);

    double verificarDatos();
    
    // Boundary conditions and endpoints
    int N;
    float a;
    float b;
    float alpha;
    float beta;
  
};