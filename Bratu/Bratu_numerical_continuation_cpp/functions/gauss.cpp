#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
using namespace std;
/**********************************************************************************/
template<typename T1, int n3, int n4>
vect<double,n3> Gauss(matrix2d<T1,n3,n4>&a1,vect<double,n3> &y1 )
{
	
	int r1,c1,n;
	int i,j,k,l;
	vect<double,n3> x1;
	//we want to solve a*x1=y;
/************************************************************************************/
// -------------- FORWARD SWEEP --------------//
	for (j=0;j<n4-1;j++)  // column number
	{
		for (i=j;i<n3-1;i++) // Row number
		{
			double mIJ;
			mIJ=  a1(i+1,j)/ a1(j,j);
/*****************************************************************/			
			for (k=0;k<n4;k++)
			{
			a1(i+1,k) = a1(i+1,k) - mIJ * a1(j,k);
			}
			y1(i+1)= y1(i+1) - mIJ * y1(j);			
		}
	}

/***************************************************************************/
/*cout<<a<<endl;
cout<<y1<<endl;	
/*****************************************************************************/	
//back substitution to find x1 that solves ax=y1
//Note: a has become an upper triangular matrix2d1 now, also, y1 has been modified:both done in the forward sweep above

  for(i=n3-1;i>=0;i--)
    {
       double sum=0.0;
        for(j=n4-1;j>i;j--)
            	{sum=sum+a1(i,j)*x1(j);}
        	x1(i)=(y1(i)-sum)/a1(i,i);
    }
/****************************************************************************/
//give tolerance
double tol;
tol=1e-10;
	for (i=0;i<n3;i++)
	{
		for(j=0;j<n4;j++)
		{
			if(fabs(a1(i,j))<tol)
			{a1(i,j)=0.0;}
		}
			if(fabs(y1(i))<tol)
			{y1(i)=0.0;}
	}
/****************************************************************************/
//check
/*	cout<<a1<<endl;
	cout<<y1<<endl;
	cout<<"x1="<<'\t'<<x1<<endl;
/****************************************************************************/
				
return x1;
}
