#include <iostream>
#include <vector>


#include "finite-difference.h"

using namespace std;

int main()
{
    FiniteDifferences finiteDifferences; 


    vector< double > k = finiteDifferences.tridiagonaLinearSystem();

    

    for (size_t i = 0; i < k.size(); i++)
    {
        cout << k.at(i) <<endl;
    }
    
    return 0; 
}

