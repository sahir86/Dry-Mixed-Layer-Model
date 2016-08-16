import numpy as np
import matplotlib.pyplot as pl
#from itertools import cycle
#lines = ["x","+","--",":","-."]
#colours = ["k","r","k","r","k"]

#from matplotlib import rc, rcParams
#rc('text',usetex=True)
#rc('font',**{'family':'serif','serif':['Computer Modern']})

def h_star_calc(delta_f, D, delta_s):
    h_star = delta_f / (D * delta_s)
    return h_star

def s_0_calc(g,Cp,T,z):
    s_0 = Cp*T + g*z
    return s_0

def sigma_calc(V,delta_f,delta_s):
    sigma = (V * delta_s )/ delta_f
    return sigma

def s_infty_calc(alpha, sigma, delta_s, s_0):
    s_infty = ((alpha - 1)/sigma) * delta_s + s_0
    return s_infty

def h_infty_calc(alpha, sigma, h_star):
    h_infty = h_star * ((alpha * sigma)/ (1 + sigma - alpha))
    return h_infty


# Define constants:

g = 9.80665 #m/s^2
Cp = 1005.0 #J/Kg/K

alpha = float(raw_input("Please enter your value for alpha:> "))
D = 4e-6 #1/s
delta_f = 40.0 #W/m^2
delta_s = 12.5e3 #J/kg
V = 0.008 #m/s
T = 273.0 + 25
z = 0.0

h_star = h_star_calc(delta_f, D, delta_s)
sigma = sigma_calc(V, delta_f, delta_s)
s_0 = s_0_calc(g,Cp,T,z)

s_infty = s_infty_calc(alpha, sigma, delta_s, s_0)

h_infty = h_infty_calc(alpha, sigma, h_star)

print "The value of alpha is: ", alpha
#print "The value of sigma is: ", sigma
print "The value of s_infty is: ", s_infty, "Joules per Kilogram."
print "The value of h_infty is: ", h_infty, "Metres."


