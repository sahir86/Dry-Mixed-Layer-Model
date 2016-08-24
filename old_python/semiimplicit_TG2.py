from sympy import *
x = Symbol('x',real=True)
l = Symbol('lambda',real=True)

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

eta = 0.5
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
print "Forming equation"
EqMat = MLHS*G - M21RHS - M22RHS*MLHS.inv()*M1RHS

print "Computing roots"
import numpy as np
Np = 20
a = np.arange(-np.pi,np.pi+1.e-7,np.pi/Np)
Gmag1 = 0.*a
Gmag2 = 0.*a
Garg1 = 0.*a
Garg2 = 0.*a
argExact = 0.*a
lval = 0.1
for i in range(a.size):
    print i
    TDet = EqMat.subs({p:a[i],l:lval}).det()
    roots = solve(TDet,G)
    if(abs(roots[1]).evalf()>abs(roots[0]).evalf()):
        tmproot = roots[0]
        roots[0] = roots[1]
        roots[1] = tmproot

    Gmag1[i], Gmag2[i] = abs(roots[0]),abs(roots[1])
    Garg1[i], Garg2[i] = arg(roots[0]),arg(roots[1])
    argExact[i] = -lval*a[i]



import pylab as pl
pl.figure()
pl.plot(a,Gmag1,'k')
pl.plot(2*pl.pi+a[0:Np+1],Gmag2[0:Np+1],'r')
pl.plot(-2*pl.pi+a[Np+1:],Gmag2[Np+1:],'r')
pl.xlabel('kh')
pl.ylabel('|G|')

pl.figure()
pl.plot(a,Garg1,'.-k')
pl.plot(2*pl.pi+a[0:Np+1],Garg2[0:Np+1],'.-r')
pl.plot(-2*pl.pi+a[Np+1:],Garg2[Np+1:],'.-r')
pl.plot(a,argExact,'.y-')
pl.plot(2*pl.pi+a[0:Np+1],argExact[0:Np+1]-2*lval*pl.pi,'b.-')
pl.plot(-2*pl.pi+a[Np+1:],argExact[Np+1:]+2*lval*pl.pi,'b.-')
pl.xlabel('kh')
pl.ylabel('arg(G)')

pl.figure()
pl.plot(a,Garg1-argExact,'.-k')
pl.plot(2*pl.pi+a[0:Np+1],Garg2[0:Np+1]-(argExact[0:Np+1]-2*lval*pl.pi),'.-r')
pl.plot(-2*pl.pi+a[Np+1:],Garg2[Np+1:]-(argExact[Np+1:]+2*lval*pl.pi),'.-r')
pl.xlabel('kh')
pl.ylabel('Dispersion error')

pl.show()

