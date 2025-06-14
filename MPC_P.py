## -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 17:30:09 2023

@author: Solorboy
#o
https://osqp.org/docs/examples/mpc.html#cvxpy

"""
#
import numpy as np
import scipy.signal
import scipy.optimize
from scipy import linalg
import matplotlib.pyplot as plt
import control as ct

#--------------------------------------------------------------------------
# subfuntion
#--------------------------------------------------------------------------

def calc_ss(A,C,nx,ny,ys):
    [V,Dd] = scipy.linalg.eig(A) 
    aux = nx-ny+1
    V2 = Dd[0:aux].reshape(3, 2) #matriz
    C= C.reshape(2, 3) #matriz
    zi = np.linalg.lstsq(C @ V2 , yss)
    zii = np.asarray(zi)
    xss = V2 @ zii[0]
    return xss

def FKalman(ny,A,C,it):
    V= .5
    # W= .5
    sM = len(A)
    PP = np.eye(sM)
    VV = np.eye(ny) * V 
    #WW = np.eye(sM) * W
    #for j in range(it):
        #PP = A * PP * A.T - A @ PP @ C.T @ np.linalg.inv(VV+C @ PP @ C.T) @ C @ PP @ A.T + WW;

    #A*PP*C.T*np.linalg.inv(VV+C*PP*C.T);
    kf = A @ PP @ C.T @ np.linalg.inv(VV+C @ PP @ C.T)
    return  kf

def ssmpc(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ap,Bp,Cp,umax,umin,dumax,yspp,uss,yss,xss,y0,u0,x0):
    #Defining the initial conditions (deviation variables)
    xpk = x0-xss  #(dimension: nx x 1)
    xmk = xpk     # (dimension: nx x 1)
    ypk = y0-yss  # (dimension: ny x 1)
    uk_1 = u0-uss # (dimension: nu x 1)
    
    Psi=[]
    ThA=[]
    for x in range(p):
        Psi = np.matrix(Psi, C*A**x)
        ThA = np.matrix(ThA, C*A**(x-1)*B)
    
    # Creating the Dynamic Matrix
    a=ThA;
    Dm=[a];
    if m >= 2:
        for iu in range(m-2):
            a = [np.zeros(ny,nu) , a[1:(p-1)*ny,:]]
            Dm = [Dm, a]
    
    b = C*B
    Ai = np.eye(nx)
    
    for kk in range(p-m):
        Ai = Ai+A**kk
        b  = [b, C*Ai*B]
    
    Theta=[Dm,[np.zeros(ny*(m-1),nu), b]]

    if m==1:
        b = C*B
        Ai = np.eye(nx)
        for kk in range(p-m):
            Ai = Ai+A**kk
            b = [b, C*Ai*B]
    
    Theta = b
    
    # Matrices Qbar and Rbar
    
    
    
    # 
    
    return ur,yr,Jk
    
    
    
    
#--------------------------------------------------------------------------
# Program
#--------------------------------------------------------------------------
# Create a MIMO transfer function object
# The transfer function from the 2nd input to the 1st output is
# g = [2.3/s         -0.7/s 
#      4.7/(9.3*s +1) 1.4/(6.8*s+1)]; %  transfer function matrix
num = [[[2.3], [-0.7]], [[4.7], [1.4]]]
den = [[[1,0], [1,0]], [[ 9.3, 1], [ 6.8, 1]]]
sys1 = ct.tf2ss(num, den)
[A,B,C,D] = ct.ssdata(sys1)

#dismash 
num2 = [[[2.3], [-0.7]], [[4.7], [1.4]]]
den2 = [[[1,0], [1,0]], [[ 9.3, 1], [ 6.8, 1]]]
sys2 = ct.tf2ss(num2, den2)
#sysd = sp.signal.cont2discrete(sys3, Ts, method='zoh', alpha=None)
[Ap,Bp,Cp,Dp] = ct.ssdata(sys2)


#
nx = len(A) # Number of system states
p = 120        # Output prediction horizon
m = 3          # Input horizon
nsim = 250     # Simulation time in sampling periods
q = np.array([0.5, 1])    # Output weights
# q=[1,5e4];   # Output weights
r = 1*np.array([1, 1])   # Input moves weights
umax = np.triu([10, 10]) # maximum value for inputs
umin = np.triu([0, 0])  # minimum value for inputs
dumax = np.triu([1, 1] ) # maximum variation for input moves
ny = 2

# Initial condition

uss = np.array([4.7, 2.65]) # steady-state of the inputs
yss = np.array([47, 52.5])  # steady-state of the outputs
xss = calc_ss(Ap,Cp,nx,ny,yss) # steady-state of the states
ys  = [43, 54]    # Set-point of the outputs
# ys=yss;


#--------------------------------------------------------------------------
# Starting simulation
#--------------------------------------------------------------------------

# Matrix H
H = Theta.T*Qbar*Theta+IM.T*Rbar*IM;
H = (H+H.T)/2

#  State observer
Kf = FKalman(ny,A,C,100);
# Auxiliary constraint matrix
Dumax=dumax;
Umax=umax-uss;
Umin=umin-uss;
for i in range(m-1):
    Umax = [Umax,umax-uss]
    Umin = [Umin,umin-uss]
    Dumax = [Dumax,dumax]

# Initial condition
u0 = [4,2].T
y0 = [40,50].T
x0 = calc_ss(Ap,Cp,nx,ny,y0) 

for i in range(nsim):
    
    ur[:,i] = uk_1 + uss
    yr[:,i] = ypk + yss
    if i<= 100:
        ys=yss
    else:
        ys=yspp;
    
    ysp=[];
    #for i=1:p:
     #   ysp = [ysp;(ys-yss)] # set-point vector (p.ny x 1)
    
    el = Psi*xmk-ysp;
    ct = el.T*Qbar*Theta-uk_1.T*Ibar.T*Rbar*IM
    c = (Psi*xmk-ysp).T*Qbar*(Psi*xmk-ysp)+uk_1.T*Ibar.T*Rbar*Ibar*uk_1
    
    # Including constraints on the input changes
    Ain = [IM,-IM]
    Bin = [Dumax+Ibar*uk_1,Dumax-Ibar*uk_1];
    options = optimoptions('quadprog','display','off');
    ukk = quadprog(H,ct,Ain,Bin,[],[],Umin,Umax,[],options);
    uk = ukk[1:nu] # receding horizon
    Jk[i] = ukk.T*H*ukk+2*ct*ukk+c
    
    # Correction of the last control input
    xmk = A*xmk+B*uk;
    ymk = C*xmk;
    if i>=101:
        #xpk = Ap*xpk+Bp*(uk+0.1*[1 .2].T)
        xpk = Ap*xpk+Bp*(uk);
        ypk = Cp*xpk;
    else:
        xpk = Ap*xpk+Bp*(uk);
        ypk = Cp*xpk;
    
    # Correction of the last measurement
    de = ypk-ymk;
    xmk = xmk+Kf*de;
    uk_1 = uk;
end








#--------------------------------------------------------------------------
# Plotting
#--------------------------------------------------------------------------

# plot results
plt.figure(1)
plt.xlabel("N Simulaçoes")
plt.ylabel("Cost Funct")
plt.title("Cost Funct")
for i in range(len(J_k)):
    plt.plot(np.arange(nsim),[pt[i] for pt in J_k],label = 'id %s'%i)
    #plt.plot(J_k[i][0],J_k[i][1])
    
plt.legend()
plt.show()



plt.figure(2)
plt.xlabel("N Simulaçoes")
plt.ylabel("Uk")
plt.title("Entrada Uk")
for i in range(len(uk)):
    plt.plot(np.arange(nsim),[pt[i] for pt in uk],label = 'id %s'%i)
plt.legend()
plt.show()


plt.figure(3)
plt.xlabel("N Simulaçoes")
plt.ylabel("Tempo Calculo")
plt.title("Tempo Calculo")
for i in range(len(tcalc)):
    plt.plot(np.arange(nsim),[pt[i] for pt in tcalc],label = 'id %s'%i)
plt.legend()
plt.show()

plt.figure(4)
plt.xlabel("N Simulaçoes")
plt.ylabel("Tempo Calculo")
plt.title("Tempo Calculo")
for i in range(len(tcalc)):
    plt.plot(np.arange(nsim),[pt[i] for pt in tcalc],label = 'id %s'%i)
plt.legend()
plt.show()



plt.figure(5)
plt.xlabel("N Simulaçoes")
plt.ylabel("Error")
plt.title("Error")
for i in range(len(e)):
    plt.plot(np.arange(nsim),[pt[i] for pt in e],label = 'id %s'%i)
plt.legend()
plt.show()


plt.figure(6)
plt.xlabel("N Simulaçoes")
plt.ylabel("Saidas")
plt.title("Saidas")
for i in range(len(yp)):
    plt.plot(np.arange(nsim),[pt[i] for pt in yp],label = 'id %s'%i)
plt.legend()
plt.show()

