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
    vect< double, N-2> F, F_guess, F1, F_lambda, F_lambda_guess, F_lambda1, negative_F_lambda, u0_cropped, u1_cropped, u0_dot_cropped, u1_dot_cropped, u1_guess_cropped, u1_guess_cropped_new, temp, temp2; // F_lambda = dF/d(\lambda) at (u0, lambda0); F_lambda_guess = dF/d(\lambda) at (u1_guess_cropped, lambda1_guess)
    
    vect<double, N> u0, u1;
    
    matrix2d<double, N-2, N-2> F_u, F_u1_guess, F_u1; //dF/dU, this is N-2xN-2 because of boundary conditions, F_u1_guess is dF/dU at (u1_guess_cropped, lambda1_guess)

    matrix2d<double, (N-2)+1, (N-2)+1> A, Newton_extended_matrix; // `A' is the extended matrix to find the tangent vector and `Newton_extended_matrix' is the one to solve Newton iterations
    
    vect<double, (N-2)+1> tangent_vector, rhs_for_tangent_vector, extended_rhs_vector, increment_vector; //the tangent vector and the rhs for tangent vector, increment_vector is (Delta_u_cropped, Delta_lambda)

    double L = 1.0; //length of the domain
    double h = (double) L/(double)(N-1);
    double oneByhSquare = (double) 1.0/ (double)(h*h);
   
    int i, j;
    double lambda0, lambda1, lambda_max, lambda0_dot, lambda1_dot, lambda1_guess, lambda1_guess_new; //continuation parameter
//    double Delta_lambda = 0.01; //step size
    double tolerance = 1e-10;

    double max_u1;

    double temp_dot_temp, tangent_vector_magnitude, temp3; // temp3 = dot((u1_guess_cropped-u0), u0_dot_cropped);//used in constructing the extended rhs vector for Newton iterations
   	
   	double ds = 1e-1; //the arc length


   lambda0=0.0;
//   lambda_max=8;
   
   int nmax = 400;

    for (int n_counter=0; n_counter < nmax; n_counter++)
   {
   		cout<<"n_counter ="<<'\t'<<n_counter<<endl;

    //get u0_cropped from u0
        for(int i = 0; i < N-2; i++)
        {
            u0_cropped(i) = u0(i+1);
        }

        cout<<"lambda0 = "<<'\t'<<lambda0<<endl;
//calculate dF/d(Lambda); 
        F_lambda = get_F_lambda(u0_cropped, lambda0);
//        cout<<"F_lambda="<<'\n'<<F_lambda<<endl;

//create the marix dF/du = F_u
        F_u = getDF_u(h, u0_cropped, lambda0);
//        cout<<"F_u="<<'\n'<<F_u<<endl;
/******************************************************************/
        //first let's find u0_dot_cropped and lambda_0_dot of the extended system

        //first get lambda_0_dot = 1/sqrt((([F_u(u0_cropped, lambda_0)]^{-1}*[F_lambda(u0_cropped, lambda_0)])^{2}+1));
        temp = LU(F_u, F_lambda);
        temp = -1*temp;
 //       cout<<"temp = "<<'\n'<<temp<<endl;

        temp_dot_temp = dot(temp, temp);

        if(n_counter<350)
        {lambda0_dot = 1.0/((double)sqrt((temp_dot_temp+1.0)));}
	   	else
    	{lambda0_dot = -1.0/((double)sqrt((temp_dot_temp+1.0)));} //change the sign of lambda_dot in order to reduce lambda and follow the branch after the turning point.
       	

        cout<<"lambda0_dot = "<<'\t'<<lambda0_dot<<endl;

       	u0_dot_cropped = lambda0_dot*temp;
        // cout<<"u0_dot_cropped = "<<'\n'<<u0_dot_cropped<<endl;

       	//now we have lambda0_dot and u0_dot
/******************************************************************/
 /********************************/
 	//Now get the guess vector for Newton iterations
       	u1_guess_cropped = ds*u0_dot_cropped;
       	u1_guess_cropped = u1_guess_cropped + u0_cropped;

       	lambda1_guess = lambda0 + (ds*lambda0_dot); 

        // cout<<"u1_guess_cropped = "<<'\n'<<u1_guess_cropped<<endl;
        cout<<"lambda1_guess = "<<'\n'<<lambda1_guess<<endl;

/********************************/
    	double Newton_iterations	=	20; //Newton iterations

	    for(int i = 0; i < Newton_iterations; i++)
	    {

			//F_lambda_guess = dF/d(\lambda) at (u1_guess_cropped, lambda1_guess)
		    //F_u1_guess is dF/dU at (u1_guess_cropped, lambda1_guess) 	
		    
	            F_guess		   = getF(h, u1_guess_cropped, lambda1_guess);
		        F_u1_guess     = getDF_u(h, u1_guess_cropped, lambda1_guess);
	            F_lambda_guess = get_F_lambda(u1_guess_cropped, lambda1_guess);

	    //construct Newton_extended_matrix to be used in Newton-Raphson iterations

	       	for(int i = 0; i < N-2; i++)
	       	{
	       		for(int j = 0; j < N-2; j++)
	       		{
	       			Newton_extended_matrix(i,j) = F_u1_guess(i,j);
	       		}
	       		Newton_extended_matrix(i, N-2) = F_lambda_guess(i); 
	       	}

	       	for(int j = 0; j < N-2; j++)
	       	{
	       		Newton_extended_matrix(N-2, j) = u0_dot_cropped(j);
	       	}

	       	Newton_extended_matrix(N-2, N-2) = lambda0_dot;

	//       	cout<<"Newton_extended_matrix = "<<'\n'<<Newton_extended_matrix<<endl;

	/********************************/

	       for(int i = 0; i < N-2; i++)
	       {
	       	extended_rhs_vector(i) = F_guess(i);
	       }

	      temp2 = (u1_guess_cropped-u0_cropped); 
	      temp3 = dot(temp2, u0_dot_cropped);

	      extended_rhs_vector(N-2) = temp3 + (lambda1_guess - lambda0)*lambda0_dot - ds;
//	      cout<<"extended_rhs_vector = "<<'\n'<<extended_rhs_vector<<endl;

	/********************************/       	
	       	//increment_vector is (Delta_u_cropped, Delta_lambda)

			//solve for increment_vector
	      	increment_vector = LU(Newton_extended_matrix, extended_rhs_vector);
	      	
	      	//increment vector is an increment in u1_guess cropped and lambda1_guess, not the rhs

	      	for (int i = 0; i < N-2; ++i)
	      	{
	      		u1_guess_cropped_new(i) = u1_guess_cropped(i) - increment_vector(i);
	      	}
	      		lambda1_guess_new = lambda1_guess - increment_vector(N-2);

	        
	        //update the guess if it is not within tolerance
			if (fabs(get_max_value(increment_vector)) < tolerance) //check for convergence       
			 {
				cout<<"converged in"<<'\t'<<i<<'\t'<<"iterations"<<endl;
				break;
	         }
	        else
	        {
				//update the guess
				u1_guess_cropped = u1_guess_cropped_new;
				lambda1_guess 	 = lambda1_guess_new;

	        }

	        if(i == (Newton_iterations-1) && fabs(get_max_value(increment_vector))> tolerance)
	        {
				cout<<"Did not converge, try a better guess."<<endl;
	        }

	    }    
/********************************/
	    u1_cropped = u1_guess_cropped_new;
	    lambda1    = lambda1_guess_new;
/******************************************************************/ 

/******************************************************************/ 
//now get u1 from u1_cropped by padding up boundary conditions
        for(int i = 1; i < N-1; i++)
        {
            u1(i) = u1_cropped(i-1);
        }
        u1(0) = 0.0; u1(N-1) = 0.0;


        max_u1 = get_max_value(u1);

//write in a txt file        

 	std::ofstream fOut;
	fOut.open("max_u_vs_lambda/max_u_vs_lambda_arclength_continuation.txt",ios::out | ios::app);
            fOut<< lambda1 << '\t' <<max_u1<<endl;		
	fOut.close();
//make u0=u1, lambda0=lambda1 and repeat
        
        u0=u1;
        lambda0=lambda1;
     }   
        
       
/*************************************************/
return 0;
}
