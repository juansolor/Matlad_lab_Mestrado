
from sys import path
path.append(r"C:\Users\altop\Documents\BCS\casadi-windows-py39-v3.5.5-64bit")
from casadi import *


# Codigos do Professor Marcio Martins

# In[30]:


import os
import control as ct
import numpy as np
import control.optimal as opt
import matplotlib.pyplot as plt
import scipy as sp
from scipy import sparse

# Add do_mpc to path. This is not necessary if it was installed via pip
#sys.path.append('../../../')

# Import do_mpc package:
#import do_mpc


# Sub-Rotinas

# In[50]:


sysc = ct.TransferFunction([1], [1, 2, 1])
sysc


# In[48]:


# Create a MIMO transfer function object
# The transfer function from the 2nd input to the 1st output is
# (3s + 4) / (6s^2 + 5s + 4).
num = [[[2.3], [-0.7]], [[4.7], [1.4]]]
den = [[[1,0], [1,0]], [[ 9.3, 1], [ 6.8, 1]]]
sys1 = ct.tf2ss(num, den)
sys1


#  Nominal model (used by MPC)

# In[55]:


Ts = 1
s = ct.TransferFunction.s
# Create a variable 's' to allow algebra operations for systems
#s = ct.tf('s') # defining the Laplace's variables
num = [[[2.3], [-0.7]], [[4.7], [1.4]]]
den = [[[1,0], [1,0]], [[ 9.3, 1], [ 6.8, 1]]]
sys2 = ct.tf(num, den)
#sysd = sp.signal.cont2discrete(sys3, Ts, method='zoh', alpha=None)
[A,B,C,D] = ct.ssdata(sys2)

A_til = sparse.csc_matrix([[1,0,0],[0,0.8981,0],[0,0,0.8632]])
B_til = sparse.csc_matrix([[1.15,-0.35],[0.5,0],[0,0.5]])
C_til = sparse.csc_matrix([[2,0,0],[0,0.9583,0.3829]])
D_til = sparse.csc_matrix([[0,0],[0,0]])
AA= np.round(A,2)
C


# model dismatch (used by dismatch)

# In[ ]:


num = [[[2.3], [-0.7]], [[4.7], [1.4]]]
den = [[[1,0], [1,0]], [[ 9.3, 1], [ 6.8, 1]]]
sys2 = ct.tf(num, den)
sys2
[A_tilp,B_tilp,C_tilp,D_tilp] = ct.ssdata(sys2);
A_tilp = sparse.csc_matrix([[1,0,0],[0,0.8981,0],[0,0,0.8632]])
B_tilp = sparse.csc_matrix([[1.15,-0.35],[0.5,0],[0,0.5]])
C_tilp = sparse.csc_matrix([[2,0,0],[0,0.9583,0.3829]])
D_tilp = sparse.csc_matrix([[0,0],[0,0]])
B_tilp


# Updating the dimensions of variables

# In[ ]:


nu = 2 # Number of manipulated inputs
ny = 2 # Number of controlled outputs
nx = 2 # Number of system states


# Tuning parameters of MPC

# In[ ]:


p = 120                  # Output prediction horizon
m = 3                    # Input horizon
nsim = 250               # Simulation time in sampling periods
q = np.array([0.5,1])    # Output weights
r =  np.array([.01,.01]) # Input weights


# Constraints

# In[ ]:


def immpc (A,B,C,nu,nx,ny):
    Ap = sparse.csc_matrix ([[A],[b],[np.zeros([nu,nx])],np.ones(nu)])
    Bp = sparse.csc_matrix ([[B],[np.ones(nu)]])
    Cp = sparse.csc_matrix ([[C],[np.zeros([ny,nu])]])
    return Ap,Bp,Cp
# Revisar esta funcion
A = sparse.csc_matrix ([[1,0,0,1.5,-0.35],[0,0.8981,0,0.5,0],[0,0,0.8632,0,0.5],[0,0,0,1,0],[0,0,0,0,1]])
B = sparse.csc_matrix ([[1.15,-0.35],[0.5,0],[0,0.5],[1,0],[0,1]])
C = sparse.csc_matrix ([[2,0,0,0,0],[0,0.9583,0.3829,0,0]])
[nx, nu] = B.shape


# In[ ]:


u0    = 10.5916
umin  = np.array([[0],[0]]) 
umax  = np.array([[10],[10]]) 
xmin  = np.array([[0],[0]])
xmax  = np.array([[20],[20]])
dumax = np.array([[1],[1]])


# In[ ]:


xmin.shape


# In[ ]:





# sub program MPC ssmpc

# In[ ]:


def ssmpc(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ap,Bp,Cp,umax,umin,dumax,ys,uss,yss,xss,y0,u0,x0):
    


#  Objective function

# In[ ]:


Q  = 2*sparse.eye(nu)
QN = Q
R  = 0.1*sparse.eye(nu)


# Initial and reference states

# In[ ]:


x0 = np.array([4.37,2.67])
xr = np.array([47,52.5])


# In[ ]:


x0


# Cost MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))

# In[ ]:


N = 10 # Horizonte controle
# - quadratic objective
P = sparse.block_diag([sparse.kron(sparse.eye(N), Q), QN,
                       sparse.kron(sparse.eye(N), R)], format='csc')
# - linear objective
q = np.hstack([np.kron(np.ones(N), -Q.dot(xr)), -QN.dot(xr),
               np.zeros(N*nu)])
# - linear dynamics
Ax = sparse.kron(sparse.eye(N+1),-sparse.eye(nx)) + sparse.kron(sparse.eye(N+1, k=-1), A)
Bu = sparse.kron(sparse.vstack([sparse.csc_matrix((1, N)), sparse.eye(N)]), B)
Aeq = sparse.hstack([Ax, Bu])
leq = sparse.hstack([-x0, np.zeros(N*nx)])
ueq = leq
# - input and state constraints
Aineq = sparse.eye((N+1)*nx + N*nu)
lineq = np.hstack([np.kron(np.ones(N+1), xmin), np.kron(np.ones(N), umin)])
uineq = np.hstack([np.kron(np.ones(N+1), xmax), np.kron(np.ones(N), umax)])
# - OSQP constraints
A_A = sparse.vstack([Aeq, Aineq], format='csc')
l = np.hstack([leq, lineq])
u = np.hstack([ueq, uineq])

# Create an OSQP object
prob = osqp.OSQP()


# In[ ]:


leq


# In[ ]:


lineq.shape


# In[ ]:





# In[ ]:




