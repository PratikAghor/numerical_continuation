#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"


//#include "getF_u.cpp"
//#include "getF.cpp"
using namespace std;

//Newton method to find a solution of DF_u*Delta_u = -f_guess
/**********************************************************************************/
template<typename T, int n>
vect<T,n> Newton(double h, vect<T,n> u_guess, double lambda, double tolerance )
{
/*******************************************************/
    vect<T,n> Delta_u, u_new;
    
    matrix2d<T,n,n> DF_u; //Jacobian matrix for Newton iterations, evaluated at u_guess
    vect<T,n>       f_guess;//f evaluated at u_guess    
    
    double Newton_iterations	=	20; //Newton iterations
    for(int i=0; i<Newton_iterations; i++)
    {
        f_guess		=  getF(h, u_guess, lambda);

        DF_u  		=  getDF_u(h, u_guess, lambda);
        
   
        //solve for Delta_u
        Delta_u = LU(DF_u, f_guess);
        
        u_new = u_guess - Delta_u; //I have already overloaded the operator - in vector.h
        
        //update the guess if it is not within tolerance
		if (fabs(get_max_value(Delta_u))< tolerance) //check for convergence       
		 {
			cout<<"converged in"<<'\t'<<i<<'\t'<<"iterations"<<endl;
			break;
         }
        else
        {
			//update the guess
			u_guess = u_new;

        }

        if(i == (Newton_iterations-1) && fabs(get_max_value(Delta_u))> tolerance)
        {
			cout<<"Did not converge, try a better guess."<<endl;
        }
    }
/*******************************************************/ 
	return u_new;   
}
