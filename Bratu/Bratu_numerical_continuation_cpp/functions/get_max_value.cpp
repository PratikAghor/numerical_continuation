#include<iostream>
#include<fstream>
#include "math.h"
#include "../headerFiles/matrix2d.h"
#include "../headerFiles/vector.h"
using namespace std;
/**********************************************************************************/
template<typename T, int n>
T get_max_value(vect<T,n> u) // n = N-2
{
    T u_max_value = u(0);
    
    for (int i=0;i<n;i++) //search for the max value
		if(u(i) > u_max_value )
		{u_max_value = u(i);
//		cout<<"i="<<i<<'\t'<<"u_max_value="<<u_max_value<<endl;
		}
	
	
	
return u_max_value;
}
