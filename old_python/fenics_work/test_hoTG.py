from dolfin import *
parameters["linear_algebra_backend"] = "PETSc"
mesh = UnitIntervalMesh(50)

nintervals = 50
mesh = UnitIntervalMesh(nintervals)

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 1.0

pbc = PeriodicBoundary()
degree = 2
V = FunctionSpace(mesh,"CG",degree, constrained_domain=pbc)
V_out = FunctionSpace(mesh,"DG",1)

u = TrialFunction(V)
v = TestFunction(V)


#u^{n+1} = u^n + dt L(u^n) + dt^2/2 L(L u^n)

Courant = 0.3
Dx = 1.0/nintervals
dt = Dx*Courant
Dt = Constant(dt)

eta = 0.6
sn = 1
c1 = 0.5*(1 + sn*(-1./3.+8*eta)**0.5)
m10 = c1
m20 = 0.5*(3-1./c1)
m21 = 0.5*(1./c1-1)
n10 = 0.5*c1**2-eta
n20 = 0.25*(3*c1-1)-eta
n21 = 0.25*(1-c1)

a = u*v*dx + eta*Dt*Dt*u.dx(0)*v.dx(0)*dx
M = assemble(a)
arhs1  = u*v*dx + m10*Dt*u*v.dx(0)*dx - n10*Dt*Dt*u.dx(0)*v.dx(0)*dx
arhs20 = u*v*dx + m20*Dt*u*v.dx(0)*dx - n20*Dt*Dt*u.dx(0)*v.dx(0)*dx
arhs21 = m21*Dt*u*v.dx(0)*dx - n21*Dt*Dt*u.dx(0)*v.dx(0)*dx

class MyExpression0(Expression):
    def eval(self, value, x):
        if(x[0]>0.25 and x[0]<0.75):
            value[0] = 1.0
        else:
            value[0] = 0.0

f0 = MyExpression0()

u = interpolate(MyExpression0(),V)
u1 = Function(V)

t = 0.0
T = 3.0
n = 0
file = File("test_TG_0.xml")
uplot = project(u,V_out)
file << uplot
xplot = interpolate(Expression("x[0]"),V_out)
file = File("test_TD_X.xml")
file << xplot
while(t<(T-dt/2)):
    t+=dt
    n+= 1
    L = assemble(action(arhs1,u))
    solve(M,u1.vector(),L)
    L = assemble(action(arhs20,u)+action(arhs21,u1))
    solve(M,u.vector(),L)
    file = File("test_TG_"+str(n)+".xml")
    uplot = project(u,V_out)
    file << uplot


