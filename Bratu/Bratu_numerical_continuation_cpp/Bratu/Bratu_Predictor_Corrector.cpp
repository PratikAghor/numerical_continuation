#include<iostream>
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
#include<math.h>


#include "../functions/lu.cpp"
#include "../functions/gauss.cpp"

#include "../functions/get_F_lambda.cpp"
#include "../functions/getDF_u.cpp"
#include "../functions/getF.cpp"

#include "../functions/Newton.cpp"
#include "../functions/get_max_value.cpp"
int main()
{

//Continuation for Bratu's problem in 1d. u"(x)+lambda*exp(u)=0 for x belonging to [0,1] and u(0)=0, u(1)=0;
/*************************************************/
    const int N = 8; //no. of grid points
    vect< double, N-2> F, F_lambda, negative_F_lambda, u0_cropped, u1_cropped, u0_dot_cropped, u_guess_cropped, u_increment_cropped;
    vect<double, N> u0, u1;
    matrix2d<double, N-2, N-2> DF_u; //dF/dU, this is N-2xN-2 because of boundary conditions. 

    double L = 1.0; //length of the domain
    double h = (double) L/(double)(N-1);
    double oneByhSquare = (double) 1.0/ (double)(h*h);
   
    int i, j;
    double lambda0, lambda1, lambda_max; //continuation parameter
    double Delta_lambda = 0.01; //step size
    double tolerance = 1e-8;

    double max_u1;
   
   lambda0=0.0;
//   lambda_max=8;
   
   int nmax = 350;

    for (int n_counter=0; n_counter < nmax; n_counter++)
   {

    //get u0_cropped from u0
        for(int i = 0; i < N-2; i++)
        {
            u0_cropped(i) = u0(i+1);
        }

        cout<<"lambda0 = "<<'\t'<<lambda0<<endl;
//calculate dF/d(Lambda); 
        F_lambda = get_F_lambda(u0_cropped, lambda0);
        cout<<"F_lambda="<<'\n'<<F_lambda<<endl;

//create the marix dF/du = F_u
        DF_u = getDF_u(h, u0_cropped, lambda0);
        cout<<"DF_u="<<'\n'<<DF_u<<endl;

//calculate u0_dot by solving F_u(u0,lambda0)*u0_dot= -F_lambda(u0,lambda0);
        negative_F_lambda = -1.0*F_lambda;
        u0_dot_cropped = LU(DF_u,negative_F_lambda);

        cout<<"u0_dot_cropped="<<'\n'<<u0_dot_cropped<<endl;

//update lambda and u_guess for Newton iterations
        lambda1 = lambda0 + Delta_lambda;
        
        u_increment_cropped = Delta_lambda*u0_dot_cropped;
        u_guess_cropped = u0_cropped + u_increment_cropped;

//do Newton iterations to get the next value of u        
        u1_cropped = Newton(h, u_guess_cropped, lambda1, tolerance);

//now get u1 from u1_cropped by padding up boundary conditions
        for(int i = 1; i < N-1; i++)
        {
            u1(i) = u1_cropped(i-1);
        }
        u1(0) = 0.0; u1(N-1) = 0.0;


        max_u1 = get_max_value(u1);

//write in a txt file        

 	std::ofstream fOut;
	fOut.open("max_u_vs_lambda/max_u_vs_lambda.txt",ios::out | ios::app);
            fOut<< lambda1 << '\t' <<max_u1<<endl;		
	fOut.close();
//make u0=u1, lambda0=lambda1 and repeat
        
        u0=u1;
        lambda0=lambda1;
     }   
        
       
/*************************************************/
return 0;
}
