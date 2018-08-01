import numpy as np  # Import numpy
import scipy
import scipy.linalg   # SciPy Linear Algebra Library

import math
import functions
##########################################

N = 8; #no. of grid points

#declare required arrays

#declare other parameters
tolerance = 1e-10;

lambda0 = 0.0;
u0 = np.zeros(N);
u1 = np.zeros(N);

u0_cropped = np.zeros(N-2);
u1_guess_cropped_new = np.zeros(N-2);
Newton_extended_matrix = np.zeros(((N-2+1), (N-2+1)), float); #Newton_extended_matrix for Newton iterations
extended_rhs_vector = np.zeros(N-2+1); #extended rhs vector used in Newton iterations

h = 1.0/(N-1);	#spacing
nmax = 500;	#maximum number of steps
ds = 1e-1; #the arc length
max_Newton_iterations = 20; #max no. of Newton iterations to converge.

for n_counter in range(0,nmax):
	print 'n_counter = ', n_counter 

	#get u0_cropped from u0 
	for i in range(0,N-2):
		u0_cropped[i] = u0[i+1];

	print 'lambda0 = ', lambda0 
	#************************************************#
	F_lambda = functions.getF_lambda(u0_cropped, lambda0)
	#	print 'F_lambda = ', F_lambda 

	F_u		 = functions.getF_u(h, u0_cropped, lambda0)
	#print 'F_u = \n', F_u

	#Solve F_u*temp = -F_lambda by LU
	# P0, L0, U0 = scipy.linalg.lu(F_u)	

	# temp = np.linalg.solve(L0,F_lambda)
	# temp = -np.linalg.solve(U0,temp);

	lu_and_piv_0 = scipy.linalg.lu_factor(F_u)
	temp = -1.0*(scipy.linalg.lu_solve(lu_and_piv_0, F_lambda))

		
	#Directly Solve F_u*temp = -F_lambda
	# temp = -1.0*(np.linalg.solve(F_u, F_lambda));
	#print 'temp = \n', temp
	    
	temp_dot_temp = np.dot(temp, temp);

	if (n_counter < 350):
		lambda0_dot = 1.0/(math.sqrt(temp_dot_temp+1.0));
	else:
		lambda0_dot = -1.0/(math.sqrt(temp_dot_temp+1.0));
		#change the sign of lambda_dot in order to reduce lambda and follow the branch after the turning point.

	print 'lambda0_dot = ', lambda0_dot

	u0_dot_cropped = lambda0_dot*temp;
	#print 'u0_dot_cropped = \n', u0_dot_cropped
	#now we have lambda0_dot and u0_dot
	#************************************************#
	#Now get the guess vector for Newton iterations
	u1_guess_cropped = u0_cropped + (ds*u0_dot_cropped);
	lambda1_guess 	 = lambda0 + (ds*lambda0_dot); 

	# print 'u1_guess_cropped = \n', u1_guess_cropped
	# print 'lambda1_guess = ', lambda1_guess
	#************************************************#
	#perform Newton iterations in order to converge to the next solution
	for m in range(0, max_Newton_iterations):
		F_guess		   = functions.getF(h, u1_guess_cropped, lambda1_guess);
		F_u1_guess     = functions.getF_u(h, u1_guess_cropped, lambda1_guess);
		F_lambda_guess = functions.getF_lambda(u1_guess_cropped, lambda1_guess);
		# print 'F_guess = \n', F_guess
		# print 'F_u1_guess = \n', F_u1_guess
		# print 'F_lambda_guess = \n', F_lambda_guess
		#******************#
		#construct Newton_extended_matrix to be used in Newton-Raphson iterations
		for i in range(0,N-2):
			for j in range(0,N-2):
				Newton_extended_matrix[i,j] = F_u1_guess[i,j];

			Newton_extended_matrix[i, N-2] = F_lambda_guess[i];

		for j in range(0,N-2):
			Newton_extended_matrix[N-2, j] = u0_dot_cropped[j];
		Newton_extended_matrix[N-2, N-2] = lambda0_dot;
		# print 'Newton_extended_matrix = \n', Newton_extended_matrix
		#******************#
		#construct extended_rhs_vector
		for i in range(0,N-2):
			extended_rhs_vector[i] = -1.0*F_guess[i];

		extended_rhs_vector[N-2] = -1.0*((np.dot((u1_guess_cropped-u0_cropped), u0_dot_cropped)) + (lambda1_guess - lambda0)*lambda0_dot - ds);
		# print 'extended_rhs_vector = \n', extended_rhs_vector
		#******************#
		#increment_vector is (Delta_u_cropped, Delta_lambda)
		#solve for increment_vector using LU: Note to use LU as it will be useful later when we are implimenting variable step or constructing a test function to see if a bifurcation point is nearby
		# P, L, U = scipy.linalg.lu(Newton_extended_matrix)

		# temp1 = np.linalg.solve(L, extended_rhs_vector)
		# increment_vector = np.linalg.solve(U, temp1);

		lu_and_piv = scipy.linalg.lu_factor(Newton_extended_matrix)
		increment_vector = scipy.linalg.lu_solve(lu_and_piv, extended_rhs_vector)
		# print 'increment_vector = \n', increment_vector

		for i in range(0, N-2):
			u1_guess_cropped_new[i] = u1_guess_cropped[i] + increment_vector[i];
			#here, I have already taken care of the minus sign while calculating the increment vector. Hence I am adding it. In the C++ code, I was calculating the increment vector without the minus sign and subtracting it. It's basically the same thing.

		lambda1_guess_new = lambda1_guess + increment_vector[N-2];

		# print 'u1_guess_cropped_new = \n', u1_guess_cropped_new

		#update the guess if it is not within tolerance
		#check for convergence
		if np.amax(abs(increment_vector))<=tolerance:
			print 'Converged in ', m, ' iterations.'
			break
		else:
			#update the guess
			# print ('updating the guess value, please wait.')
			u1_guess_cropped = u1_guess_cropped_new;
			lambda1_guess 	 = lambda1_guess_new;

		# print 'u1_guess_cropped = \n', u1_guess_cropped

		if (m == max_Newton_iterations-1) and (np.amax(abs(increment_vector))<tolerance):
			print 'Did not converge. Try a better guess.'
		#******************#
	u1_cropped = u1_guess_cropped_new;
	lambda1    = lambda1_guess_new;
	#******************#
	#now get u1 from u1_cropped by padding up boundary conditions
	for i in range(0,N-2):
		u1[i+1] = u1_cropped[i];

	u1[0] = 0.0; 
	u1[N-1] = 0.0;

	max_u1 = np.amax(u1);

	#Use file module (more general, enables custom formatting):
	f = open("data/Bratu_max_u1_vs_lambda_arclength_continuation.txt", "a")

	f.write("%f" % float(lambda1)) #Write lambda1
	f.write("\t %f" % float(max_u1)) # Write max_u1
	f.write("\n") #New line
	f.close()

	#update u0 and lambda0
	u0 = u1;
	lambda0 = lambda1;



# Plot the solution:
#Import plotting modules:
from pylab import figure, grid, hold, legend, savefig
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
fig = plt.figure() #Create a figure instance
ax = fig.add_subplot(111)

with open("data/Bratu_max_u1_vs_lambda_arclength_continuation.txt") as f:
	rows = f.readlines()
	f.close()

parameter = []
l_inf_norm = []

for row in rows:
	p = row.split()
	parameter.append(float(p[0]));
	l_inf_norm.append(float(p[1]));

x1 = np.array(parameter);
y1 = np.array(l_inf_norm);

ax.plot(x1, y1, lw=3)
mpl.pyplot.grid(b=None, which='major', axis='both')
ax.set_title("Bratu arclength continuation")    
ax.set_xlabel('$\lambda$', fontsize=18)
ax.set_ylabel('$max. u$', fontsize=18)
savefig("plots/Bratu_arclength_continuation.png", bbox_inches='tight', dpi=100)
plt.tight_layout()
plt.show()

