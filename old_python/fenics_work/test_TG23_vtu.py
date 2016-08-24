from dolfin import *
parameters["linear_algebra_backend"] = "PETSc"

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
V_out = FunctionSpace(mesh,"DG",degree)

u = TrialFunction(V)
v = TestFunction(V)


#u^{n+1} = u^n + dt L(u^n) + dt^2/2 L(L u^n)

Courant = 0.6 #1.0
Dx = 1.0/nintervals
dt = Dx*Courant
Dt = Constant(dt)

eta = 0.431180320501
root = 1
c1 = 0.5*(1+root*(-1.0/3.0+8*eta))**0.5

alpha = Constant(0.5)

a = u*v*dx + 0.5*alpha*Dt*Dt*u.dx(0)*v.dx(0)*dx
M = assemble(a)
arhs = u*v*dx + Dt*u*v.dx(0)*dx - (1-alpha)*Dt*Dt*0.5*u.dx(0)*v.dx(0)*dx

class MyExpression0(Expression):
    def eval(self, value, x):
        if(x[0]>0.25 and x[0]<0.75):
            value[0] = 1.0
        else:
            value[0] = 0.0

f0 = MyExpression0()

u = interpolate(MyExpression0(),V)

t = 0.0
T = 10.0
n = 0
file = File("test_TG_0.xml")
uplot = project(u,V_out)
file << uplot
xplot = interpolate(Expression("x[0]"),V_out)
file = File("test_TD_X.xml")
file << xplot
# while(t<(T-dt/2)):
#     t+=dt
#     n+= 1
#     L = assemble(action(arhs,u))
#     solve(M,u.vector(),L)
#     file = File("test_TG_"+str(n)+".xml")
#     uplot = project(u,V_out)
#     file << uplot


