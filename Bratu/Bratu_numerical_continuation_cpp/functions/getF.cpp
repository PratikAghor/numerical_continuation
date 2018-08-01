#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
using namespace std;
/**********************************************************************************/
template<typename T, int n>
vect<T,n> getF(double h, vect<T,n> u_cropped, double lambda ) // n = N-2
{
    //u_cropped is u with boundaries cropped off. Hence, u(i) = u_cropped(i-1);
    double oneByhSquare = (double) 1.0/ (double)(h*h);
    vect<T,n> F; //F, this is N-2x1 (i.e., nx1) because of boundary conditions. 
    int i;
    
    for (i=1;i<n-1;i++) 
    {
        F(i)  = oneByhSquare*(u_cropped(i+1)-2.0*u_cropped(i)+u_cropped(i-1))+lambda*exp(u_cropped(i));
    }
        //write the first and the last entry of the F vector
        F(0) = oneByhSquare*(u_cropped(1)-2.0*u_cropped(0))+lambda*exp(u_cropped(0)); // this is actually evaluated at the node 1, hence u(1). u(0) = 0 because of the boundary condition.
        F(n-1) = oneByhSquare*(-2.0*u_cropped(n-1)+u_cropped(n-2))+lambda*exp(u_cropped(n-1)); // this is actually evaluated at the node N-2, hence u(n). u(N-1)= u(n+1) = 0 because of the boundary condition.
//cout<<"F="<<'\n'<<F_u<<endl;
    
    return F;
}
