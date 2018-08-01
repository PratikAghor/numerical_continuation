#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"


#include "getF_u.cpp"
#include "getF.cpp"
#include "Newton.cpp"
#include "get_max_value.cpp"
#include "lu.cpp"


using namespace std;

//Check the Newton method code to find a solution of F_u*Delta_u = -F
/**********************************************************************************/
int main()
{
	double x_guess,y_guess; //variables

	const int n = 2;
	
	vect<double, n> f_guess, delta_answer, answer, answer_guess;


	//build DF = [[df(0)/dx, df(0)/dy;
	//				[df(1)/dx, df(1)/dy]]

	matrix2d<double, n, n> DF;

	//define tolerance
	double tolerance = 1e-8;

	//define a guess for f_0
	answer_guess(0)	=	1.0;
	answer_guess(1)	=	-1.0;


	double Newton_iterations	=	10; //Newton iterations

	for (int i = 0; i < Newton_iterations; i++)
	{
			
			x_guess	= answer_guess(0);
			y_guess	= answer_guess(1);

			//build f
			f_guess(0)	=	4.0*x_guess*y_guess*y_guess + 2.0*x_guess*y_guess - 2.0;
			f_guess(1)	=	4.0*x_guess*x_guess*y_guess + x_guess*x_guess - 2.0*y_guess;

			

			DF(0,0)	=	4.0*y_guess*y_guess + 2.0*y_guess;	
			
			DF(0,1)	=	8.0*x_guess*y_guess + 2.0*x_guess;
			
			DF(1,0)	=	8.0*x_guess*y_guess + 2.0*x_guess;	
			
			DF(1,1)	=	4.0*x_guess*x_guess - 2.0;	

			cout<<"DF = "<<'\n'<<DF<<endl;

		//apply Newton's method to find the steady state solution
		delta_answer	=	LU(DF, f_guess);

	//	cout <<"iterations = " <<'\t'<<i<<'\t'<<"delta_answer = "<<delta_answer<<endl;

		answer 		 =	answer_guess - delta_answer;	//update

		if (fabs(get_max_value(delta_answer))< tolerance) //check for convergence
		{
			cout<<"converged in"<<'\t'<<i<<'\t'<<"iterations"<<endl;
			break;
		}

		else
		{
			//update the guess
			answer_guess = answer;
		}


	}
	cout<<"answer = " <<'\n'<<answer<<endl;
/**********************************/
	return 0;
}