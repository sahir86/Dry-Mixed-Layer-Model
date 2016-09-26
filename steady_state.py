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
params = AttrDict(constants) # Try to find a more standard way of doing this.

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

def ak_delta_f_calc_theta(Q_rad,ak_h_infty):
    ak_delta_f = - Q_rad * ak_h_infty
    return ak_delta_f

def ak_delta_f_calc_s(Cp,Q_rad,ak_h_infty):
    ak_delta_f = - Cp * Q_rad * ak_h_infty
    return ak_delta_f

def ak_phi_infty_calc(phi_0,ak_delta_f,V,beta):
    ak_phi_infty = phi_0 - (ak_delta_f)/(V*(beta+1))
    return ak_phi_infty

def h_quad_solver(a,b,c):
    ak_plus_h_infty = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    ak_minus_h_infty = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    return ak_plus_h_infty, ak_minus_h_infty

# Note that the coefficients below correspond to delta_f = -Q_rad*h_infty. This is only true for potential temperature. 
def h_quad_coefficients_theta(D,gamma,Q_rad,V,beta,phi_0_ft,phi_0):

    a = D*(gamma - (Q_rad)/(V*(beta + 1)))
    b = D*(phi_0_ft - phi_0) + (beta*Q_rad)/(beta+1)
    c = 0.0
    return a,b,c

def h_quad_coefficients_s(D,gamma,Q_rad,V,beta,phi_0_ft,phi_0,Cp):

    a = D*(Cp * gamma - (Cp * Q_rad)/(V*(beta + 1)))
    b = D*(phi_0_ft - phi_0) + (beta * Cp * Q_rad)/(beta+1)
    c = 0.0
    return a,b,c

def ak_unit_converter(gamma,Q_rad):
    ak_gamma_per_metre = gamma / (1000.0)
    ak_qrad_per_second = Q_rad / (24.0 * 60.0 * 60.0)        
    return ak_gamma_per_metre, ak_qrad_per_second 

def rho_calc(p_0, R_air, T):
    ak_rho = p_0 / (R_air * T)
    return ak_rho

#ak_rho = rho_calc(ak_params.p_0,ak_params.R_air,ak_params.theta_0)


# Compute h_infty and s_infty using Ann-Kristin model functions.
#ak_gamma_per_metre, ak_qrad_per_second = ak_unit_converter(ak_params.gamma,ak_params.Q_rad)
#
#ak_plus_h_infty, ak_minus_h_infty = h_quad_solver(*h_quad_coefficients_theta(ak_params.D,ak_gamma_per_metre,ak_qrad_per_second,ak_params.V,ak_params.beta,ak_params.theta_0_ft,ak_params.theta_0))
#
#ak_delta_f = ak_delta_f_calc_theta(ak_qrad_per_second,ak_plus_h_infty)
#ak_theta_infty = ak_phi_infty_calc(ak_params.theta_0,ak_delta_f,ak_params.V,ak_params.beta)
#
## Check to see what effective value of D is being used:
#
#ak_D = - ak_qrad_per_second / (ak_gamma_per_metre * ak_plus_h_infty)
#ak_wft = - ak_qrad_per_second / (ak_gamma_per_metre) 
#
#print "Plus h_infty equals:> ", ak_plus_h_infty 
#print "Minus h_infty equals:> ", ak_minus_h_infty
#
#print "Delta_f equals:> ", ak_delta_f
#print "Theta_infty equals:> ", ak_theta_infty
#print "Effective D equals:> ", ak_D
#print "Effective wft equals:> ", ak_wft
#

# Ask user to choose between theta or S.

choice = raw_input("Would you like to compute for potential temperature (t) or static energy (s)? (default is t):> ") or "t"

# Lowest and highest values of beta. Note that beta should never equal -1.0.
start_beta = float(raw_input("Please enter the lowest value for beta (default is 0.025):> ") or "0.025") # If the user enters nothing, use 0.025 as the default value.  
end_beta = float(raw_input("Please enter the highest value for beta (default is 2.5):> ")or "2.5") # Use 2.5 as the default value if the user enters nothing. 


beta = np.linspace(start_beta,end_beta,1000) # Create an array of beta values using the input values taken from the user. Note that 1000 divisions are used. 

Q_rad = np.linspace(-1.0,-6.0,6)

if choice == "t":

    fig1 = pl.figure("h_vs_beta")
    #pl.figure(1)

    fig2 = pl.figure("theta_vs_beta")
    #pl.figure(2)
    for i in range(Q_rad.size):

        # Compute variables using functions defined above. 

        ak_gamma_per_metre, ak_qrad_per_second = ak_unit_converter(ak_params.gamma,Q_rad[i])

        ak_plus_h_infty, ak_minus_h_infty = h_quad_solver(*h_quad_coefficients_theta(ak_params.D,ak_gamma_per_metre,ak_qrad_per_second,ak_params.V,beta,ak_params.theta_0_ft,ak_params.theta_0))

        ak_delta_f = ak_delta_f_calc_theta(ak_qrad_per_second,ak_plus_h_infty)
        ak_theta_infty = ak_phi_infty_calc(ak_params.theta_0,ak_delta_f,ak_params.V,beta)

        # Check to see what effective value of D is being used:

        #ak_D = - ak_qrad_per_second / (ak_gamma_per_metre * ak_plus_h_infty)
        #ak_wft = - ak_qrad_per_second / (ak_gamma_per_metre) 
        

    #    labelstr_sigma =str(r"$\sigma =$ ") + "%.2f" % sigma
        labelstr_Q_rad =str(r"$\quad Q_\mathrm{rad}  =$ ") + "%.2f" % Q_rad[i] + str(r" $\mathrm{K/Day}$") # This is incomplete, need to see how to use a second legend.
        combined_labelstr =  labelstr_Q_rad 
        #pl.figure(1) # Switch to figure 1.
        pl.figure("h_vs_beta")
            
        pl.plot(beta,ak_plus_h_infty,colours[i],label=combined_labelstr)
        pl.title("Steady State Height vs beta" ,fontsize=26)
        pl.xlabel(r"$\beta$", fontsize=24)
        pl.ylabel(r"$\mathrm{h}_\infty \quad \mathrm{[m]}$", fontsize=24)
        pl.legend(ncol=1,loc = 'upper left')
        #pl.xlim(lvec[0]-1,lvec[-1]+1)
        
        #pl.figure(2) # Switch to figure 2. 
        pl.figure("theta_vs_beta") 
        
    #    s_infty_kJ = s_infty / 1000.0 # Convert values in Joules to kilo Joules. 
        pl.plot(beta,ak_theta_infty,colours[i],label=combined_labelstr)
        pl.title("Steady State Potential Temperature vs beta", fontsize=26)
        pl.xlabel(r"$\beta$", fontsize=24)
        pl.ylabel(r"$\theta_\infty \quad \mathrm{[K]}$", fontsize=24)
        pl.legend(ncol=1,loc = 'upper left')

elif choice == "s":
    
    fig1 = pl.figure("h_vs_beta")
    #pl.figure(1)

    fig2 = pl.figure("s_vs_beta")
    #pl.figure(2)
    for i in range(Q_rad.size):

        # Compute variables using functions defined above. 

        ak_gamma_per_metre, ak_qrad_per_second = ak_unit_converter(ak_params.gamma,Q_rad[i])
        
        s_0 = s_0_calc(ak_params.g,ak_params.Cp,ak_params.T,ak_params.z)
        s_0_ft = s_0_calc(ak_params.g,ak_params.Cp,ak_params.theta_0_ft,ak_params.z) 

        ak_plus_h_infty, ak_minus_h_infty = h_quad_solver(*h_quad_coefficients_s(ak_params.D,ak_gamma_per_metre,ak_qrad_per_second,ak_params.V,beta,s_0_ft,s_0,ak_params.Cp))
        ak_delta_f = ak_delta_f_calc_s(ak_params.Cp,ak_qrad_per_second,ak_plus_h_infty)
        ak_s_infty = ak_phi_infty_calc(s_0,ak_delta_f,ak_params.V,beta)
        ak_min_delta_f = np.amin(ak_delta_f)
        ak_max_delta_f = np.amax(ak_delta_f)

        # Check to see what effective value of D is being used:

        #ak_D = - ak_qrad_per_second / (ak_gamma_per_metre * ak_plus_h_infty)
        #ak_wft = - ak_qrad_per_second / (ak_gamma_per_metre) 
        

    #    labelstr_sigma =str(r"$\sigma =$ ") + "%.2f" % sigma
        labelstr_Q_rad =str(r"$\quad Q_\mathrm{rad}  =$ ") + "%.2f" % Q_rad[i] # This is incomplete, need to see how to use a second legend.
        labelstr_range_delta_f = str(r"$\quad$") + "%.2f" % ak_min_delta_f + str(r"$ \leq \Delta \mathrm{F} \leq $ ") + "%.2f" % ak_max_delta_f 
        combined_labelstr =  labelstr_Q_rad + labelstr_range_delta_f 
        #pl.figure(1) # Switch to figure 1.
        pl.figure("h_vs_beta")
            
        pl.plot(beta,ak_plus_h_infty,colours[i],label=combined_labelstr)
        pl.title("Steady State Height vs beta" ,fontsize=26)
        pl.xlabel(r"$\beta$", fontsize=24)
        pl.ylabel(r"$\mathrm{h}_\infty \quad \mathrm{[m]}$", fontsize=24)
        pl.legend(ncol=1,loc = 'upper left')
        #pl.xlim(lvec[0]-1,lvec[-1]+1)
        
        #pl.figure(2) # Switch to figure 2. 
        pl.figure("s_vs_beta") 
        
        ak_s_infty_kJ = ak_s_infty / 1000.0 # Convert values in Joules to kilo Joules. 
        pl.plot(beta,ak_s_infty_kJ,colours[i],label=combined_labelstr)
        pl.title("Steady State Dry Static Energy vs beta", fontsize=26)
        pl.xlabel(r"$\beta$", fontsize=24)
        pl.ylabel(r"$\mathrm{s}_\infty \quad \mathrm{[kJ kg^{-1}]}$", fontsize=24)
        pl.legend(ncol=1,loc = 'upper left')
#fig1.savefig("height_vs_beta.pdf")

pl.show(fig1)
#pl.close(fig1)

#fig2.savefig("s_infinity_vs_beta.pdf")
        
pl.show(fig2)
#pl.close(fig2)


