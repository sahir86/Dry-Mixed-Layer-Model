# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
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

# # Create a text file to save the max |G| values
# file = open('maxvals.txt','w')
#Local coordinates
B0 = 1-x
B1 = x

#Lagrange polynomials
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

#local basis functions
N = [L2_0(x),L2_1(x),L2_2(x)]

def cell_integrate(F):
    return integrate(F,(x,0,1))

nN = len(N)

#P1 local mass matrix
M = Matrix(nN,nN, lambda i,j: 0)
for a in range(nN):
    for b in range(nN):
        M[a,b] = cell_integrate(N[a]*N[b])

#P1 local derivative matrix
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

#etavec = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7]
# etavec = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
# etavec = np.linspace(0.43,0.5,num=21) # unstable to stable transition occurs between 0.472 and 0.4755
# etavec = np.linspace(0.470,0.4755,num=21) # Further zoomed in.
etavec = np.linspace(0.472,0.473,num=21)
letavec = len(etavec) #Length of etavec

lvec = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
# lvec = [1,2,3,4,5,6,7,8]
# lvec = [10,15,20,25,30,35,40,45]
lenlvec = len(lvec)

# Create a text file to save the max |G| values
file = open("abs(G) maxvals for eta = " + str(etavec[0]) + " to " + str(etavec[-1]) + " & sigma = " + str(lvec[0]) + " to " + str(lvec[-1]) + ".txt","w")

for e in range(letavec):
    eta = etavec[e]
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
    n_divs = 15
    # a = np.arange(-np.pi,np.pi,np.pi/n_divs)
    a = np.linspace(-np.pi,np.pi,num=n_divs*2)
    Gmag1 = 0.*a
    Gmag2 = 0.*a
    Garg1 = 0.*a
    Garg2 = 0.*a
    argExact = 0.*a
    # lvec = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    # lenlvec = len(lvec)
    fig1 = pl.figure(1) # Generate first figure for the |G| plots.
    fig2 = pl.figure(2) # Generate second figure for the arg(G) plots.
    fig3 = pl.figure(3) # Generate third figure for dispersion error plots.
    file.write("%s\n" % " ")
    for lam in range(lenlvec):
        lval = lvec[lam]
        # print 'sigma','=',lval
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

        # lines = ["x","+","s","D","^","o","h","p"]
        lines = ["x","+","^","d","o","s","h","p"]
        colours = ["k","r","b","g","y","m","c","k"]

        rc('text',usetex=True)
        rc('font',**{'family':'serif','serif':['Computer Modern']})
        labelstr = "$\sigma =$ " + str(lval) # For |G| vs kh plot and arg(G) vs kh plot.
        maxstr = "Max |G| = " + str(max(max(Gmag1),max(Gmag2))) + " for eta = " + str(eta) + " and sigma = " + str(lval)
        file.write("%s\n" % maxstr) # Write the max |G| values to the text file.
        file.flush() # Command to update changes to the text file.

    #     pl.figure(1)

    #     pl.plot(-2*np.pi+a[n_divs:],Gmag2[n_divs:],colours[lam]+lines[lam])
    #     pl.plot(a,Gmag1,colours[lam]+lines[lam],label=labelstr)
    #     pl.plot(2*np.pi+a[0:n_divs],Gmag2[0:n_divs],colours[lam]+lines[lam])

    #     pl.title("Magnitude of G vs $kh$ for $\eta$ = " + str(eta), fontsize=26)
    #     pl.xlabel("$kh$", fontsize=24)
    #     pl.ylabel("$|G|$", fontsize=24)
    #     pl.xlim(-7,7)
    #     pl.legend(ncol=lenlvec/2,loc = 'lower center')

    # 	pl.figure(2)

    #     pl.plot(-2*np.pi+a[n_divs:],Garg2[n_divs:],colours[lam]+lines[lam])
    #     pl.plot(a,Garg1,colours[lam]+lines[lam],label=labelstr)
    #     pl.plot(2*np.pi+a[0:n_divs],Garg2[0:n_divs],colours[lam]+lines[lam])

    #     pl.plot(-2*np.pi+a[n_divs:],argExact1,colours[lam])
    #     pl.plot(a,argExact,colours[lam])
    #     pl.plot(2*np.pi+a[0:n_divs],argExact2,colours[lam])

    #     pl.title("Argument of G vs $kh$ for $\eta$ = " + str(eta), fontsize=26)
    #     pl.xlabel("$kh$", fontsize=24)
    #     pl.ylabel("$arg(G)$", fontsize=24)
    #     pl.xlim(-7,7)
    #     pl.legend(ncol=lenlvec/2,loc='lower center')

        # pl.figure(3)

        # pl.plot(-2*np.pi+a[n_divs:],Garg2[n_divs:] - argExact1,colours[lam]+lines[lam])
        # pl.plot(a,Garg1-argExact,colours[lam]+lines[lam],label=labelstr)
        # pl.plot(2*np.pi+a[0:n_divs],Garg2[0:n_divs] - argExact2,colours[lam]+lines[lam])

        # pl.title("Dispersion Error vs $kh$ for $\eta$ = " + str(eta), fontsize=26)
        # pl.xlabel("$kh$",fontsize=24)
        # pl.ylabel(r"$\varepsilon_{\Phi}$",fontsize=24)
        # pl.xlim(-7,7)
        # pl.legend(ncol=lenlvec/2,loc='lower center')

    # fig1.savefig("G vs kh for eta = " + str(eta) + ".pdf")
    # pl.close(fig1)

    # fig2.savefig("arg(G) vs kh for eta = " + str(eta) + ".pdf")
    # pl.close(fig2)

    # fig3.savefig("Dispersion error vs kh for eta = " + str(eta) + ".pdf")
    # pl.close(fig3)

# pl.show(block=False)
file.close()
