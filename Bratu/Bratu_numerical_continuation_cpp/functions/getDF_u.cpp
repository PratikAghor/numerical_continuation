#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
using namespace std;
/**********************************************************************************/
template<typename T, int n>
matrix2d<double,n,n> getDF_u(double h, vect<T,n> u_cropped, double lambda ) //n = N-2
{
    //u_cropped is u with boundaries cropped off. Hence, u(i) = u_cropped(i-1);

    double oneByhSquare = (double) 1.0/ (double)(h*h);
    matrix2d<double, n, n> F_u; //dF/dU, this is n-1xn-1 because of boundary conditions. 
    int i;
    
    for (i=1;i<n-1;i++) 
    {
        F_u(i,i-1)  = oneByhSquare;
        F_u(i,i)    = -2*oneByhSquare + lambda*exp(u_cropped(i));  
        F_u(i,i+1)  = oneByhSquare;         
    }
    //write the first and the last row of F_u
         F_u(0,0) = -2*oneByhSquare + lambda*exp(u_cropped(0));
         F_u(0,1) = oneByhSquare; 
            
         F_u(n-1, n-2) = oneByhSquare;    
         F_u(n-1, n-1) = -2*oneByhSquare + lambda*exp(u_cropped(n-1));

//cout<<"F_u="<<'\n'<<F_u<<endl;
    
    return F_u;
}
