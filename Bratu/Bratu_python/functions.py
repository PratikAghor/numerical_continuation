import numpy as np  # Import numpy
import math

def getF(h, u0_cropped, lambda0):
    """
    function specifying the nonlinear system of equations
    
    """
    n = np.size(u0_cropped);
    #Generate the array to fill in:
    F = np.zeros(n)

    oneByhSquare = 1.0/(h*h)

    #Specify F:
    for i in range(1, n - 1):
        F[i] = oneByhSquare * (u0_cropped[i+1] - 2.0*u0_cropped[i] + u0_cropped[i-1]) + lambda0*math.exp(u0_cropped[i]);

    F[0] = oneByhSquare*(u0_cropped[1] - 2.0*u0_cropped[0]) + lambda0*math.exp(u0_cropped[0]);
    F[n - 1] = oneByhSquare*(-2.0*u0_cropped[n-1] + u0_cropped[n-2]) + lambda0*math.exp(u0_cropped[n-1]);
    return F

def getF_u(h, u0_cropped, lambda0):
    """
    function to find the matrix F_u = del(F_{i})/del(u_{j})
    
    """
    n = np.size(u0_cropped);
    #Generate the array to fill in:
    F_u = np.zeros((n,n))

    oneByhSquare = 1.0/(h*h)
  
    for i in range(1, n - 1):
        #Specify F_u:
        F_u[i, i-1] = oneByhSquare;
        F_u[i, i]   = -2*oneByhSquare + lambda0*math.exp(u0_cropped[i]);
        F_u[i, i+1] = oneByhSquare;
    F_u[0,0] = -2*oneByhSquare + lambda0*math.exp(u0_cropped[0]);
    F_u[0,1] = oneByhSquare; 
            
    F_u[n-1, n-2] = oneByhSquare;    
    F_u[n-1, n-1] = -2*oneByhSquare + lambda0*math.exp(u0_cropped[n-1]);
    return F_u

def getF_lambda(u0_cropped, lambda0):
    """
    function to get the derivative of F with respect to the parameter lambda0 = del(F)/del(lambda0)
    
    """
    n = np.size(u0_cropped);
    #Generate the array to fill in:
    F_lambda = np.zeros(n)

    for i in range(0, n):
        #Specify F_lambda:
        F_lambda[i] = math.exp(u0_cropped[i]);
    return F_lambda
