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
import cvxopt
from cvxopt import matrix
import cvxpy as cp
# Import math Library
import math
import time
#--------------------------------------------------------------------------
# subfuntion
#--------------------------------------------------------------------------

def calc_ss(A,C,nx,ny,ys):
    [V,Dd] = scipy.linalg.eig(A) 
    aux = nx-ny+1
    V2 = Dd[0:aux].reshape(3, 2) #matriz
    C = C.reshape(2, 3) #matriz
    zi = np.linalg.lstsq(C @ V2 , ys)
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


def get_mat_psi(A, N):
    psi = np.matrix(np.zeros((0, A.shape[1])))

    for i in range(1, N + 1):
        psi = np.vstack((psi, A ** i))

    return psi


def get_mat_gamma(A, B, N):
    (nx, nu) = B.shape
    gamma = np.zeros((nx, nu)) + B

    for i in range(1, N):
        tmat = (A ** i) @ B + gamma[-nx:, :]
        gamma = np.vstack((gamma, tmat))

    return gamma


def get_mat_theta(A, B, N):
    AiB = B
    (nx, nu) = B.shape
    theta = np.kron(np.eye(N), AiB)

    tmat = np.zeros((nx, 0))

    for i in range(1, N):
        t = np.zeros((nx, nu)) + B
        for ii in range(1, i + 1):
            t += (A ** ii) @ B
        tmat = np.hstack((t, tmat))

    for i in range(1, N):
        theta[i * nx:(i + 1) * nx, :i] += tmat[:, -i:]

    return theta


def generate_du_constraints_mat(G, h, N, nu, mindu, maxdu):

    if maxdu is not None:
        tG = np.matrix(np.eye(N * nu))
        th = np.kron(np.ones((N * nu, 1)), maxdu)
        G = np.vstack([G, tG])
        #h = np.vstack([h, th])
        h =  th

    if mindu is not None:
        tG = np.matrix(np.eye(N * nu)) * -1.0
        th = np.kron(np.ones((N * nu, 1)), mindu * -1.0)
        G = np.vstack([G, tG])
        h = np.vstack([h, th])

    return G, h


def generate_u_constraints_mat(G, h, N, nu, u0, minu, maxu):

    # calc F
    F = np.matrix(np.zeros((nu, N * nu + 1)))
    for i in range(N * nu):
        if maxu is not None:
            tF = np.zeros((nu, N * nu + 1))
            tF[0, i] = 1.0 / maxu[0, 0]
            tF[1, i] = 1.0 / maxu[0, 1]
            tF[0, -1] = -1.0
            tF[1, -1] = -1.0
            F = np.vstack([F, tF])

        if minu is not None:
            tF = np.zeros((nu, N * nu + 1))
            tF[0, i] = 1.0 / minu[0, 0]
            tF[1, i] = 1.0 / minu[0, 1]
            tF[0, -1] = -1.0
            tF[1, -1] = -1.0
            F = np.vstack([F, tF])

    f1 = F[:, -1]
    f2 = F[:, -2]
    f = np.concatenate((f1, f2), axis=1)
    #  print(f)

    # calc SF
    SF = np.matrix(np.zeros((F.shape[0], F.shape[1] - 1)))
    for i in range(0, N * nu):
        for ii in range(i, N * nu):
            SF[:, i] += F[:, ii]
    #  print(SF)

    SF1 = SF[:, 0]
    #  print(SF1)

    th = -SF1 * u0.T - f
    #  print(th)

    G = np.vstack([G, SF])
    #h = np.vstack([h, th])
    h = th

    return G, h


def generate_x_constraints_mat(G, h, N, nx, u0, theta, kappa, minx, maxx):

    # calc G
    cG = np.matrix(np.zeros((0, N * nx + 1)))
    for i in range(N):
        if maxx is not None:
            for ii in range(maxx.shape[1]):
                tG = np.zeros((1, N * nx + 1))
                tG[0, i * nx + ii] = 1.0 / maxx[0, ii]
                tG[0, -1] = -1.0
                cG = np.vstack([cG, tG])

        if minx is not None:
            for ii in range(minx.shape[1]):
                tG = np.zeros((1, N * nx + 1))
                tG[0, i * nx + ii] = 1.0 / minx[0, ii]
                tG[0, -1] = -1.0
                cG = np.vstack([cG, tG])

    #  print(cG)

    tau = cG[:, :-1]
    g = cG[:, -1]
    #  print(tau)
    #  print(g)

    tmpG = tau * theta
    tmph = -tau * kappa - g
    #  print(tmpG)
    #  print(tmph)

    #  print(cG)

    G = np.vstack([G, tmpG])
    #h = np.vstack([h, tmph])
    h = tmph

    return G, h


def model_predictive_control(A, B, N, Q, R, T, x0, ysp,u0,mindu=None, maxdu=None, minu=None, maxu=None, maxx=None, minx=None):

    #(nx, nu) = B.shape

    #du = np.matrix([0.0] * N).T

    #psi = get_mat_psi(A, N)
    #gamma = get_mat_gamma(A, B, N)
    #theta = get_mat_theta(A, B, N)

    #QQ = scipy.linalg.block_diag(np.kron(np.eye(N), Q))
    #RR = scipy.linalg.block_diag(np.kron(np.eye(N), R))

    #H = theta.T * QQ * theta + RR
    #H =  (QQ @ theta.T)@theta + RR
    #g = - theta.T @ QQ @ (T - (psi @ x0) - (gamma @ u0))
    #g = - QQ @ theta.T @ ((psi @ x0) - (gamma @ u0))
    # print(H)
    # print(g)
    # print(u0)

    #G = np.zeros((0, nu * N))
    #h = np.zeros((0, nu))

    #G0, h0 = generate_du_constraints_mat(G, h, N, nu, mindu, maxdu)
    #G1, h1 = generate_u_constraints_mat(G, h, N, nu, u0, minu, maxu)

    #kappa = psi * x0 + gamma * u0
    #G2, h2 = generate_x_constraints_mat(G, h, N, nx, u0, theta, kappa, minx, maxx)

    #  print(H)
    #  print(g)
    #  print(G)
    #  print(h)

    #sol = pyecosqp.ecosqp(H, g, A=G0, B=h0)
    #sol = cvxopt.solvers.qp(matrix(H), matrix(g), A=matrix(G0), B=matrix(H), Aeq=matrix(A), Beq=matrix(B))
    # du = np.matrix(sol["x"]).T
    
    
    
    x = cp.Variable((nx, N + 1))
    u = cp.Variable((nu, N))

    cost = 0
    constr = []
    
    x01 = x0 - ysp.reshape(3,1)
    
    for t in range(N):
        cost += cp.sum_squares(x[:, t]) + cp.sum_squares(u[:, t])
        constr += [x[:, t + 1] == A @ x[:, t] + B @ u[:, t], cp.norm(u[:, t], "inf") <= 1]
        # sums problem objectives and concatenates constraints.
    constr += [x[:, N] == 0, x[:, 0] == x0[:, 0]]
    problem = cp.Problem(cp.Minimize(cost), constr)
    start = time.time()
    problem.solve()
    elapsed_time = time.time() - start
    print("calc time:{0} [sec]".format(elapsed_time))
    print(cost.value)
    

    if problem.status == cp.OPTIMAL:
        du = [u[0, :].value,u[1, :].value]
        x = [x[0, :].value,x[1, :].value,x[2, :].value]
        u1 =  u0 + np.matrix((du[0][0],du[1][0])).T
    else: 
        x = [x[0, :].value,x[1, :].value,x[2, :].value]
        du = [u[0, :].value,u[1, :].value]
        u1 =  u0 + np.matrix((du[0][0],du[1][0])).T
    print(du)
    return x, u1, du,elapsed_time,cost



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
Q = np.diag([0.5, 1])    # Output weights
# q=[1,5e4];   # Output weights
nu = 2
R = 1*np.eye(nu)   # Input moves weights
umax = np.matrix([10, 10]).T # maximum value for inputs
umin = np.matrix([0, 0]).T  # minimum value for inputs
dumax = np.matrix([1, 1] ).T # maximum variation for input moves
ny = 2

# Initial condition

uss = np.array([4.7, 2.65]).T # steady-state of the inputs
yss = np.array([47, 52.5]).T  # steady-state of the outputs
xss = calc_ss(Ap,Cp,nx,ny,yss) # steady-state of the states
ys  = np.array([43, 54]).T    # Set-point of the outputs
# ys=yss;

#--------------------------------------------------------------------------
# Starting simulation
#--------------------------------------------------------------------------

#  State observer
Kf = FKalman(ny,A,C,100)
# Auxiliary constraint matrix
Dumax = dumax
Umax = umax-uss
Umin = umin-uss


# Initial condition
u0 = np.matrix([4,2]).T
u00 = [4,2]
y0 = np.matrix([40,50]).T
y00 = [40,50]
x0 = np.array(calc_ss(Ap,Cp,nx,ny,y0))
uk_1 = u00 - uss # (dimension: nu x 1)
# Defining the initial conditions (deviation variables)
xpk = x0 - xss #(dimension: nx x 1)
xmk = xpk    # (dimension: nx x 1)
ypk = y00 - yss # (dimension: ny x 1)
#Create list and memory of variavel
ur = []
yr = []
Hp = 10
du = np.matrix([1.0, 0.25] * Hp).T

#--------------------------------------------------------------------------
# Simulation
#--------------------------------------------------------------------------

ur = []
yr = []
xr = []
dr = []
tcalr = []
Jr = []
dur = []

for i in range(nsim):
    
   
    if i<= 100:
        ys = yss
    else:
        ys = yss*1.1
    
    xsp = calc_ss(A,Cp,nx,ny,ys)
    x, u, du,tcal,Jcust = model_predictive_control(A, B, Hp, Q, R, du, x0, xsp, u0, mindu=None, maxdu = Dumax, minu = Umin, maxu = Umax, maxx=None, minx=None)
    
    # Correction of the last control input
    x1 = np.delete(x, [-1,-2,-3,-4,-5,-6,-7,-8],-1)
    
    xmk = A @ x1 + B @ u
    ymk = C @ xmk
    
    if i>=101:
        #xpk = Ap*xpk+Bp*(uk+0.1*[1 .2].T)
        xpk = Ap @ xpk+Bp @ u
        ypk = Cp @ xpk
    else:
        xpk = Ap@xpk+Bp @ u
        ypk = Cp@xpk
    
    # Correction of the last measurement
    de = ypk - ymk
    xmk = xmk + Kf@de
    uk_1 = u
    x0 = np.array(np.delete(xmk, [-1,-2],-1))
    
    # memoria das variavel
    ur += [uk_1]
    yr += [ypk]
    xr += [xmk]
    dr += [de]
    tcalr += [tcal]
    Jr += [Jcust.value]
    dur += [du]






#--------------------------------------------------------------------------
# Plotting
#--------------------------------------------------------------------------

# plot results
plt.figure(1)
plt.xlabel("N Simulaçoes")
plt.ylabel("Cost Funct")
plt.title("Cost Funct")
#for i in range(len(J_k)):
#    plt.plot(np.arange(nsim),[pt[i] for pt in J_k],label = 'id %s'%i)
    #plt.plot(J_k[i][0],J_k[i][1])
    
plt.legend()
plt.show()



plt.figure(2)
plt.xlabel("N Simulaçoes")
plt.ylabel("Uk")
plt.title("Entrada Uk")
aux3 = []
aux4 = []
aux5 = []
for i in range(len(ur)):
    aux3.append(ur[0])
    aux4.append(ur[i][0])
    aux5.append(ur[i][1])
plt.plot(np.arange(nsim),aux4)
plt.plot(np.arange(nsim),aux5)
plt.legend()
plt.show()


plt.figure(3)
plt.xlabel("N Simulaçoes")
plt.ylabel("Tempo Calculo")
plt.title("Tempo Calculo")
#for i in range(len(tcalc)):
#    plt.plot(np.arange(nsim),[pt[i] for pt in tcalc],label = 'id %s'%i)
plt.legend()
plt.show()

plt.figure(4)
plt.xlabel("N Simulaçoes")
plt.ylabel("Tempo Calculo")
plt.title("Tempo Calculo")
#for i in range(len(tcalc)):
#    plt.plot(np.arange(nsim),[pt[i] for pt in tcalc],label = 'id %s'%i)
plt.legend()
plt.show()



plt.figure(5)
plt.xlabel("N Simulaçoes")
plt.ylabel("Error")
plt.title("Error")
#for i in range(len(e)):
 #   plt.plot(np.arange(nsim),[pt[i] for pt in e],label = 'id %s'%i)
plt.legend()
plt.show()


plt.figure(6)
plt.xlabel("N Simulaçoes")
plt.ylabel("Saidas")
plt.title("Saidas")
#for i in range(len(yp)):
#    plt.plot(np.arange(nsim),[pt[i] for pt in yp],label = 'id %s'%i)
plt.legend()
plt.show()

