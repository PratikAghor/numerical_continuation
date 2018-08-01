#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"


#include "getDF_u.cpp"
#include "getF.cpp"
#include "Newton.cpp"
#include "get_max_value.cpp"
#include "lu.cpp"


using namespace std;

/**********************************************************************************/
int main()
{
	//check the dot product function

	const int n =3;

	vect<double, n> a, b;
	double c;

	for (int i = 0; i < n; i++)
	 {
	 	a(i) = i+1;
	 	b(i) = i+3;
	 } 

	 c = dot(a, b);

	 cout<<"a = "<<'\n'<<a<<endl;
	 cout<<"b = "<<'\n'<<b<<endl;
	 cout<<"c = <a.b> "<<'\t'<<c<<endl;
/**********************************/
	return 0;
}