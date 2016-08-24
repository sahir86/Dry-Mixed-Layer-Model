from sympy import *
x = Symbol('x',real=True)
l = Symbol('lambda',real=True)

# Import relevant modules
import numpy as np
import matplotlib.pyplot as pl
pl.ioff() # Stop outputs from coming up on the display.
from matplotlib import rc, rcParams
from numpy import size # Only for debuggin purposes.

# Close any existing figures
pl.close('all')

# Local coordinates
B0 = 1-x
B1 = x

# Lagrange polynomials
def L1_0(x):
    return 1-x
def L1_1(x):
    return x
def L2_0(x):
    return (2*x-1)*(x-1)
def L2_1(x):
    return 4*x*(1-x)
def L2_2(x):
    return (2*x-1)*x

# local basis functions
N = [L2_0(x),L2_1(x),L2_2(x)]

def cell_integrate(F):
    return integrate(F,(x,0,1))

nN = len(N)

# P1 local mass matrix

M = Matrix(nN,nN, lambda i,j: 0)
for a in range(nN):
    for b in range(nN):
        M[a,b] = cell_integrate(N[a]*N[b])

# P1 local derivative matrix
D = Matrix(nN,nN, lambda i,j: 0)
for a in range(nN):
    for b in range(nN):
        D[a,b] = cell_integrate(diff(N[a],x)*N[b])

#P1 local diffusion matrix
K = Matrix(nN,nN, lambda i,j: 0)
for a in range(nN):
    for b in range(nN):
        K[a,b] = cell_integrate(diff(N[a],x)*diff(N[b],x))

p = Symbol('phi',real=True)

assert(nN==3)
S = Matrix([[exp(-I*p/2),0],[0,Rational(1,1)],[exp(I*p/2),0]])

# Use np.linspace to create lvec in order to generate an array object correctly.
# lvec = [10,15,20,25,30,35,40,45]
# lvec = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
# lvec = np.linspace(0.8,0.8,num=1)

# lvec arrays for the full binary search - only use once bugs are ironed out.
# lvec1 = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
lvec1 = np.array(np.linspace(0.1,1.0,num=25))
lvec2 = np.array(np.linspace(1.0,4.0,num=10))
lvec3 = np.array([5.0,6.0,7.0,8.0,9.0])
lvec4 = np.array([10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0])
lvec5 = np.array([30.0,35.0,40.0,45.0,50.0,55.0,60.0])
# lvec3 = np.array([10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0])
lvec = np.concatenate([lvec1,lvec2,lvec3,lvec4,lvec5])
# lvec = lvec1
np.save("lvec_" + str(lvec[0]) + " to " + str(lvec[-1]) + ".npy",lvec)
lenlvec = len(lvec)

# Initialise parameters for the binary search
final_etavec = 0.*lvec # Initialise vector for holding final binary search values of eta.

# Create a text file to save the max |G| values
file = open("binary search |G| maxvals for sigma = " + str(lvec[0]) + " to " + str(lvec[-1]) + ".txt","w")
fig1 = pl.figure(1) # Generate first figure for the critical eta plots.
for lam in range(lenlvec):
    lval = lvec[lam]
    # print 'sigma','=',lval
    file.write("%s\n" % " ")

    # Reset and initialise binary search parameters before entering while loop.
    upper_eta = 0.55#0.47275 # The critical value of eta above which the scheme is stable for all lambda.
    lower_eta = 0.35
    delta_eta = (upper_eta - lower_eta)
    eta = delta_eta/2. # Initialise eta
    maxg = 1.2

    # Begin conditional loop.
    while delta_eta > 1e-6 or maxg > 1.0:
        # print 'eta','=',eta
        sn = 1
        c1 = 0.5*(1 + sn*(-1./3.+8*eta)**0.5)
        m10 = c1
        m20 = 0.5*(3-1./c1)
        m21 = 0.5*(1./c1-1)
        n10 = 0.5*c1**2-eta
        n20 = 0.25*(3*c1-1)-eta
        n21 = 0.25*(1-c1)

        # Z_i - eta*dt**2 Z_{i,tt} =
        # Z_0 + dt*sum_{j=0}^{i-1} m_{ij}Z_{j,t} +
        #    dt**2*sum_{j=0}^{i-1} n_{ij}Z_{j,tt}

        # Z_t = -Z_x so MZ_t = DZ
        # Z_{tt} = -Z_{xt} = Z_{xx} so MZ_{tt} = -KZ

        MLHS = S.H*(M + eta*l*l*K)*S
        M1RHS = S.H*(M + l*m10*D - l*l*n10*K)*S
        M21RHS = S.H*(M + l*m20*D - l*l*n20*K)*S
        M22RHS = S.H*(l*m21*D - l*l*n21*K)*S

        #Simplification
        for a in range(2):
            for b in range(2):
                MLHS[a,b] = expand(MLHS[a,b],complex=True)
                M1RHS[a,b] = expand(M1RHS[a,b],complex=True)
                M21RHS[a,b] = expand(M21RHS[a,b],complex=True)
                M22RHS[a,b] = expand(M22RHS[a,b],complex=True)


        #MLHS*Z^{n+1} = M21RHS*Z^n + M22RHS*Z_1^n
        #             = M21RHS*Z^n + M22RHS*MLHS^{-1}M1RHS

        G = Symbol('G')
        # print "Forming equation"
        EqMat = MLHS*G - M21RHS - M22RHS*MLHS.inv()*M1RHS

        # print "Computing roots"
        n_divs = 10
        # a = np.arange(-np.pi,np.pi,np.pi/n_divs)
        a = np.linspace(-np.pi,np.pi,num=n_divs*2)
        Gmag1 = 0.*a
        Gmag2 = 0.*a
        Garg1 = 0.*a
        Garg2 = 0.*a
        argExact = 0.*a

        # fig1 = pl.figure(1) # Generate first figure for the |G| plots.
        # fig2 = pl.figure(2) # Generate second figure for the arg(G) plots.
        # fig3 = pl.figure(3) # Generate third figure for dispersion error plots.
        for i in range(a.size):
            # print i
            TDet = EqMat.subs({p:a[i],l:lval}).det()
            roots = solve(TDet,G)
            if(abs(roots[1]).evalf()>abs(roots[0]).evalf()):
                tmproot = roots[0]
                roots[0] = roots[1]
                roots[1] = tmproot

            Gmag1[i], Gmag2[i] = abs(roots[0]),abs(roots[1])
            Garg1[i], Garg2[i] = arg(roots[0]),arg(roots[1])
            argExact[i] = -lval*a[i]

        argExact1 = argExact[n_divs:] + 2*lval*np.pi #argExact for -2pi < a < -pi
        argExact2 = argExact[0:n_divs] - 2*lval*np.pi #argExact for pi < a < pi

        lines = ["x","+","s","D","^","o","h","p"]
        colours = ["k","r","b","g","y","m","c","k"]

        rc('text',usetex=True)
        rc('font',**{'family':'serif','serif':['Computer Modern']})
        labelstr = "$\sigma =$ " + str(lval) # For |G| vs kh plot and arg(G) vs kh plot.
        maxg = max(max(Gmag1),max(Gmag2))
        maxstr = "Max |G| = " + str(maxg) + " for sigma = " + str(lval) + " and eta = " + str(eta) + " with delta_eta = " + str(delta_eta)
        if maxg > 1.0:
            lower_eta = eta
            delta_eta = (upper_eta - lower_eta)
            eta = lower_eta + delta_eta/2.
        elif maxg <= 1.0:
            upper_eta = eta
            delta_eta = (upper_eta - lower_eta)
            eta = upper_eta - delta_eta/2.

        # maxstr = "Max |G| = " + str(maxg) + " for lambda = " + str(lval) + " and eta = " + str(eta) + " with upper_eta = " + str(upper_eta) + " and lower_eta = " + str(lower_eta)
        file.write("%s\n" % maxstr) # Write the max |G| values to the text file.
        file.flush() # Command to update changes to the text file.
    else:
        final_etavec[lam] = eta
        np.save("final_etavec_" + str(lvec[0]) + " to " + str(lvec[-1]) + ".npy",final_etavec)
        finalstr = "    Critical eta = " + str(eta) + " for sigma = " + str(lval) + " with max |G| = " + str(maxg)

        file.write("%s\n" % finalstr)
        file.flush() # Command to update changes to the text file.

pl.figure(1)

pl.plot(lvec,final_etavec,'b')
pl.plot(lvec,final_etavec,'ro')
pl.title("Critical values of $\eta$ for " + str(lvec[0]) + " $\leq \sigma \leq$ " + str(lvec[-1]),fontsize=26)
pl.xlabel("$\sigma$", fontsize=24)
pl.ylabel("Critical $\eta$", fontsize=24)
pl.xlim(lvec[0]-1,lvec[-1]+1)

fig1.savefig("Critical eta vs sigma for sigma = " + str(lvec[0]) + " to " + str(lvec[-1]) + ".pdf")
pl.close(fig1)
# pl.show(block=False)

file.close() # Close maxvals text file.
