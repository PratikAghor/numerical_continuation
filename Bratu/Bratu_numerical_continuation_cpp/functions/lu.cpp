#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
using namespace std;
/**********************************************************************************/
template<typename T, int n1, int n2>
vect<double,n1> LU(matrix2d<T,n1,n2>&a,vect<double,n1> &b )
{	

	int r1,c1;
	int i,j,k,l;
	

	matrix2d<double,n1,n2> L,U;
	vect<double,n1> y,x;

	for(i=0;i<n1;i++)
	{
		L(i,i)=1;
				
		for(j=0;j<i;j++)
		{	double sum1=0;
			for(k=0;k<j;k++)
			{
			sum1=sum1+L(i,k)*U(k,j);
			}
			L(i,j)= (a(i,j)-sum1)/U(j,j);
		}
		
		for(j=i;j<n2;j++)
		{
			double sum2=0;
			for(k=0;k<i;k++)
			{
			sum2=sum2+L(i,k)*U(k,j);
			}
			
			U(i,j)=a(i,j)-sum2;
		}
		
	}

/**********************************************************************/
/*cout<<"L="<<'\n'<<L<<endl;
cout<<"U="<<'\n'<<U<<endl;
/***** FINDING y; Ly=b*************************************************/
 
	for(i=0;i<n1;i++)
	    {                                        //forward subtitution method
		double sum=0;
		for(j=0;j<i;j++)
		{sum=sum+L(i,j)*y(j);}
		y(i)=(b(i)-sum)/L(i,i);
	    }
/********** FINDING x; Ux=y*********************************************/

    for(i=n1-1;i>=0;i--)
    {
       double sum=0;
        for(j=n2-1;j>i;j--)
            {sum=sum+U(i,j)*x(j);}
        x(i)=(y(i)-sum)/U(i,i);
    }
/***********************************************************************/
//disply solution
return x;
}

