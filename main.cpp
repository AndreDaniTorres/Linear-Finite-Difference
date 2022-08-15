#include <iostream>
#include <vector>


#include "finite-differences.h"

using namespace std;

int main()
{
    FiniteDifferences finiteDifferences; 
    finiteDifferences.setN(9);
    finiteDifferences.setA(1.0);
    finiteDifferences.setB(2.0);
    finiteDifferences.setAlpha(1.0);
    finiteDifferences.setBeta(2.0);

    vector< double > k = finiteDifferences.tridiagonaLinearSystem();

    

    for (size_t i = 0; i < k.size(); i++)
    {
        cout << k.at(i) <<endl;
    }
    
    return 0; 
}

