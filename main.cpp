#include <iostream>
#include <vector>


#include "finite-difference.h"

using namespace std;

int main()
{
    FiniteDifferences finiteDifferences; 

    // Solve system
    vector< double > k = finiteDifferences.tridiagonaLinearSystem();

    // Get X 
    vector< double > x = finiteDifferences.getX();

    // Print and Save solution file
    finiteDifferences.saveSolutions(k, x);

    return 0; 
}

