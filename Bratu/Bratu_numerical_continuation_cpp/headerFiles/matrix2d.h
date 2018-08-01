#ifndef MATRIX2D_H //Include guard
#define MATRIX2D_H
#include<iostream>
#include<stdlib.h>
#include "vector.h"

using namespace std;

//matrix2d class
/******************************************************************************************************************/
/******************************************************************************************************************/
template <typename T, int row, int column>   //row by column
class matrix2d
{
	public:
  

	
	T operator()(const int i, const int j) const 
	{
		if(i>=0 && i<row && j>=0 && j<column)
			{return mat2[i*column + j];}

		else
			{cout<<'\t'<<"Index out of bounds"<<endl;}

	}//{return mat2[i*column+j];}
	


	T& operator()(const int i,const int j) 
	{
		if(i>=0 && i<row && j>=0 && j<column)
			{return mat2[i*column + j];}

		else
			{cout<<'\t'<<"Index out of bounds"<<endl;}

	}//{return mat2[i*column+j];}	

	void clear(); //clear function
	void copyOther(const matrix2d &other);//copy function	
		
	matrix2d(); //CONSTRUCTOR	
	~matrix2d(); //DESTRUCTOR	

	matrix2d(const matrix2d &other) // Copy constructor
	{
	copyOther(other);
//      	cout<<"matrix2d copy constructor invoked "<<endl;			
	}
	
	matrix2d & operator =(const matrix2d & other) // Assn. operator
	{	if(this != &other)
		{
			clear();
			copyOther(other);
//      			cout<<"matrix2d assignment operator invoked "<<endl;								
		}
	return (*this);
	}

	
private:
	T *mat2 ;


	

};


	template<typename T, int row, int column>  
	void matrix2d<T,row,column>::clear() //clear function
	{
        free(mat2);
//  	delete [] mat2;
	    mat2 = NULL;
	}
	
	template<typename T, int row, int column>  
	void matrix2d<T,row,column>::copyOther(const matrix2d &other)//copy function
	{
			mat2 = new T [row*column]; 
			for(int i=0;i<row*column;i++)
			{
				mat2[i]=other.mat2[i];
			}
	}
	
	template<typename T, int row, int column> 
	matrix2d<T,row,column>::matrix2d() //CONSTRUCTOR
	{
//        cout << "matrix2d Constructor invoked: " << endl;
		mat2 = new T[row*column];
		for(int i=0;i<row*column;i++)
		{mat2[i]= 0;}
	}

	template<typename T, int row, int column> 
	matrix2d<T,row,column>::~matrix2d()//DESTRUCTOR
	{
//		clean up allocated memory
		clear();
//		delete [] mat2;
//		mat2 = NULL;
//	cout << "matrix2d Destructor invoked: " << endl;
	}


	template<typename T, int row, int column> //overloading output operator 
	std::ostream & operator<<  (std::ostream& os, matrix2d<T,row,column> &mymatrix2d)
	{	
		int i,j;
		for(i=0;i<row;i++)
		{	
			for(j=0;j<column;j++)
			{
			//std::cout<<'\t';
			os<<mymatrix2d(i,j);
			std::cout<<'\t';
			}
			std::cout<<'\n';
		}
	
	return os;
	}

	template<typename T, int row, int column>  //overloading input operator 
	std::istream & operator>>  (std::istream& is, matrix2d<T,row,column> &mymatrix2d)
	{	
		int i,j;
		for(i=0;i<row;i++)
		{	
			
			for(j=0;j<column;j++)
			{	
			is>>mymatrix2d(i,j);
			}
		std::cout<<'\n';	
		}
	return is;
	}

	
	template<typename T, int row, int column> //overloading + operator 
	matrix2d<T,row,column> operator+  ( matrix2d<T,row,column> &a, matrix2d<T,row,column> &b)
	{
	int i,j;
	matrix2d<T,row,column> answer;
		for(i=0;i<row;i++)
		{	
			for(j=0;j<column;j++)
			{
			answer(i,j)=a(i,j)+b(i,j);
			}	
		
		}
	return answer;
	}

	
	template<typename T, int row, int column> //overloading - operator 
	matrix2d<T,row,column> operator-  ( matrix2d<T,row,column> &a, matrix2d<T,row,column> &b)
	{
	int i,j;
	matrix2d<T,row,column> answer;
		for(i=0;i<row;i++)
		{	
			for(j=0;j<column;j++)
			{
			answer(i,j)=a(i,j)-b(i,j);
			}	
		
		}
	return answer;
	}


	
	template<typename T, int row1, int column1, int row2, int column2> //overloading * operator 
	matrix2d<T,row1,column2> operator*  ( matrix2d<T,row1,column1> &a, matrix2d<T,row2,column2> &b)
	{
	int i,j,k,c1,c2,r1,r2;
	matrix2d<T,row1,column2> answer;

		if(column1!=row2)
        	{std::cout<<"matrix2d multiplication invalid. number of columns of A must be equal to number of rows of B"<<std::endl;}
	
		else
		{
			for(i=0;i<row1;i++)
			{
				for(k=0;k<column2;k++)
				{
					answer(i,k)=0;
					for(j=0;j<column1;j++)
					{
					answer(i,k)=answer(i,k)+a(i,j)*b(j,k);
					}
				}
			
			}
		}
	return answer;
	}

//matrix2d-VECTOR MULTIPLICATION

	template<typename T, int row1, int column1, int row2> //overloading * operator 
	vect<T,row2> operator*  ( matrix2d<T,row1,column1> &a, vect<T,row2> &b)
	{
	int i,j,k,c1,c2,r1,r2;
	vect<T,row2> answer;

		if(column1!=row2)
        	{std::cout<<"matrix2d-vector multiplication invalid. number of columns of A must be equal to number of rows of B"<<std::endl;}
	
		else
		{
			for(i=0;i<row1;i++)
			{
				
					answer[i]=0;
					for(j=0;j<column1;j++)
					{
					answer[i]=answer[i]+a(i,j)*b[j];
					}
				
			
			}
		}
	return answer;
	}


//SYMMETRIC matrix2d FLAG 
	template<typename T, int row, int column>
	bool symFlag  ( matrix2d<T,row,column> &a)
	{
		bool flag=1;	
		int i,j;
		for(i=0;i<row;i++)
		{	
			for(j=0;j<column;j++)
			{
				if(a(i,j)==a(j,i))
				{
				flag=flag*1;
				}	
				else
				{
				flag=flag*0;
				}
			}
		}
	return flag;
	}

	
//TRANSPOSE OF A MATRIX
	template<typename T, int row, int column>
	matrix2d<T,row,column> Transpose (matrix2d<T,row,column> &a)
	{
	int i,j;
	matrix2d<T,row,column> answer;
		for(i=0;i<row;i++)
		{	
			for(j=0;j<column;j++)
			{
			answer(i,j)=a(j,i);
			}	
		
		}
	return answer;
	}
/******************************************************************************************************************/
/******************************************************************************************************************/
#endif
