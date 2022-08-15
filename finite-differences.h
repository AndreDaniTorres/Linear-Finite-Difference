#include <iostream>
#include <vector>

// using namespace std;


    
// clase 
class FiniteDifferences
{
    public:
    FiniteDifferences(); 

    void centeredDifference(std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
    std::vector<double> tridiagonaLinearSystem();
    void setN(int);
    void setA(float);
    void setB(float);
    void setAlpha(float);
    void setBeta(float);
    float p(float x);
    float q(float x);
    float r(float x);
    

    
    private:
    int N;
    float a;
    float b;
    float alpha;
    float beta;
  
};