#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
using namespace std;
/**********************************************************************************/
template<typename T, int n>
vect<T,n> get_F_lambda(vect<T,n> u_cropped, double lambda)
{

    vect<T,n> F_lambda; //dF/d (lambda), this is N-2x1 because of boundary conditions. 
    int i;
    
    for (i=0;i<n;i++) 
    {
        F_lambda(i) = exp(u_cropped(i));
    }

//cout<<"F="<<'\n'<<F_u<<endl;
    
    return F_lambda;
}
