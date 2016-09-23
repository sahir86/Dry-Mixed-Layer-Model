import numpy as np
import matplotlib.pyplot as pl
from itertools import cycle

#lines = ["x","+","s","D","^","o","h","p"]
colours = ["k","r","b","g","y","m","c","k"]

from matplotlib import rc, rcParams
# Use LaTex to format figure axes labels. 
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})

#import params

# Import my 'common' module that includes relevant constants etc
from common import *

# Close any existing figures
pl.close('all')

# Define constants:
# This allows you to access stored dictionary constants using the dot method. 
#params = AttrDict(constants) # Try to find a more standard way of doing this.

# Define the main functions for computing the mathematical expressions.
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

# Ann-Kristin model functions.
ak_params = AttrDict(ak_consts)

def ak_delta_f_calc(Q_rad,ak_h_infty):
    ak_delta_f = - Q_rad * ak_h_infty
    return ak_delta_f

def ak_phi_infty_calc(phi_0,ak_delta_f,V,beta):
    ak_phi_infty = phi_0 - (ak_delta_f)/(V*(beta+1))
    return ak_phi_infty

def h_quad_solver(a,b,c):
    ak_plus_h_infty = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    ak_minus_h_infty = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    return ak_plus_h_infty, ak_minus_h_infty

# Note that the coefficients belwo correspond to delta_f = -Q_rad*h_infty. This is only true for potential temperature. 
def h_quad_coefficients_theta(D,gamma,Q_rad,V,beta,phi_0_ft,phi_0):

    a = D*(gamma - (Q_rad)/(V*(beta + 1)))
    b = D*(phi_0_ft - phi_0) + (beta*Q_rad)/(beta+1)
    c = 0.0
    return a,b,c

def ak_unit_converter(gamma,Q_rad):
    ak_gamma_per_metre = gamma / (1000.0)
    ak_qrad_per_second = Q_rad / (24.0 * 60.0 * 60.0)        
    return ak_gamma_per_metre, ak_qrad_per_second 

# Compute h_infty and s_infty using Ann-Kristin model functions.
ak_gamma_per_metre, ak_qrad_per_second = ak_unit_converter(ak_params.gamma,ak_params.Q_rad)

ak_plus_h_infty, ak_minus_h_infty = h_quad_solver(*h_quad_coefficients_theta(ak_params.D,ak_gamma_per_metre,ak_qrad_per_second,ak_params.V,ak_params.beta,ak_params.theta_0_ft,ak_params.theta_0))

ak_delta_f = ak_delta_f_calc(ak_qrad_per_second,ak_plus_h_infty)
ak_theta_infty = ak_phi_infty_calc(ak_params.theta_0,ak_delta_f,ak_params.V,ak_params.beta)

# Check to see what effective value of D is being used:

ak_D = - ak_qrad_per_second / (ak_gamma_per_metre * ak_plus_h_infty)
ak_wft = - ak_qrad_per_second / (ak_gamma_per_metre) 

print "Plus h_infty equals:> ", ak_plus_h_infty 
print "Minus h_infty equals:> ", ak_minus_h_infty

print "Delta_f equals:> ", ak_delta_f
print "Theta_infty equals:> ", ak_theta_infty
print "Effective D equals:> ", ak_D
print "Effective wft equals:> ", ak_wft

"""
# Lowest and highest values of alpha. Note that alpha should be > 0 and < 3.5. Never equal to 3.5. 
# Consider writing an if statement with an error message if alpha >= 3.5 is entered. 
start_alpha = float(raw_input("Please enter the lowest value for alpha:> ") or "0.025") # If the user enters nothing, use 0.025 as the default value.  
end_alpha = float(raw_input("Please enter the highest value for alpha:> ")or "2.5") # Use 2.5 as the default value if the user enters nothing. 


alpha = np.linspace(start_alpha,end_alpha,100) # Create an array of alpha values using the input values taken from the user. Note that 100 divisions are used. 

delta_f = np.linspace(30.0,50.0,5)

fig1 = pl.figure("h_vs_alpha")
#pl.figure(1)

fig2 = pl.figure("s_vs_alpha")
#pl.figure(2)

for i in range(delta_f.size):

    # Compute variables using functions defined above. 

    h_star = h_star_calc(delta_f[i], params.D, params.delta_s)
    
    sigma = sigma_calc(params.V, delta_f[i], params.delta_s)
    
    s_0 = s_0_calc(params.g,params.Cp,params.T,params.z)
    
    s_infty = s_infty_calc(alpha, sigma, params.delta_s, s_0)
    
    h_infty = h_infty_calc(alpha, sigma, h_star)
    

    labelstr_sigma =str(r"$\sigma =$ ") + "%.2f" % sigma
    labelstr_delta_f =str(r"$\quad \Delta F =$ ") + "%.2f" % delta_f[i] # This is incomplete, need to see how to use a second legend.
    combined_labelstr = labelstr_sigma + "  " + labelstr_delta_f 
    #pl.figure(1) # Switch to figure 1.
    pl.figure("h_vs_alpha")
        
    pl.plot(alpha,h_infty,colours[i],label=combined_labelstr)
    pl.title("Steady State Height vs alpha" ,fontsize=26)
    pl.xlabel(r"$\alpha$", fontsize=24)
    pl.ylabel(r"$\mathrm{h}_\infty \quad \mathrm{[m]}$", fontsize=24)
    pl.legend(ncol=1,loc = 'upper left')
    #pl.xlim(lvec[0]-1,lvec[-1]+1)
    
    #pl.figure(2) # Switch to figure 2. 
    pl.figure("s_vs_alpha") 
    
    s_infty_kJ = s_infty / 1000.0 # Convert values in Joules to kilo Joules. 
    pl.plot(alpha,s_infty_kJ,colours[i],label=combined_labelstr)
    pl.title("Steady State Dry Static Energy vs alpha", fontsize=26)
    pl.xlabel(r"$\alpha$", fontsize=24)
    pl.ylabel(r"$\mathrm{s}_\infty \quad \mathrm{[kJ kg^{-1}]}$", fontsize=24)
    pl.legend(ncol=1,loc = 'upper left')

#print "The value of alpha is: ", alpha
#print "The value of sigma is: ", sigma
#print "The value of s_infty is: ", s_infty, "Joules per Kilogram."
#print "The value of h_infty is: ", h_infty, "Metres."

#print "alpha_vec = ", alpha
#print " h_infty_vec = ", h_infty

 
 
#fig1.savefig("height_vs_alpha.pdf")

pl.show(fig1)
#pl.close(fig1)

#fig2.savefig("s_infinity_vs_alpha.pdf")
        
pl.show(fig2)
#pl.close(fig2)

"""
