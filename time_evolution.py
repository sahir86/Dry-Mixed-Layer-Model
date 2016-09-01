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

def s_0_calc(g,Cp,T,z):
    s_0 = Cp*T + g*z
    return s_0

# An initial guess for s_n when t = 0. Try to find a better way of doing this.
def s_init_calc(s_0,delta_s):
    s_init = s_0 + delta_s / 2.0
    return s_init

def h_n_plus_one_calc(alpha,delta_f,delta_s,s_n,s_0,D,h_n,delta_t):
    h_n_plus_one = ((alpha * delta_f) / (delta_s - s_n + s_0) - D*h_n) * delta_t + h_n 
    return h_n_plus_one

def s_n_plus_one_calc(delta_f,alpha,V,s_n,s_0,h_n,delta_t):
    s_n_plus_one = ((delta_f*(alpha - 1) - V*(s_n - s_0)) / h_n) * delta_t + s_n
    return s_n_plus_one

def time_steps_calc(t_max,delta_t):
    time_steps = int( t_max / delta_t)
    return time_steps
# Lowest and highest values of alpha. Note that alpha should be > 0 and < 3.5. Never equal to 3.5. 

# For time evolving problem. Ask users for alpha and initial h values.
alpha = float(raw_input("Please enter a value for alpha:> ") or "1.0") # If the user enters nothing, use 1.5 as the default value.  
h_init = float(raw_input("Please enter an initial value for h:> ") or "800.0") # If nothing is entered, use 10m for the initial height
t_max = int(raw_input("Please enter a maximum time value:> ") or "2000000")
delta_t = float(raw_input("Please enter a time step size:> ") or "50.0")
#alpha = np.linspace(start_alpha,end_alpha,100) # Create an array of alpha values using the input values taken from the user. Note that 100 divisions are used. 

#delta_f = np.linspace(30.0,50.0,5)

# Compute number of time_steps
time_steps = time_steps_calc(t_max,delta_t)
# Arrays for storing the computed solutions at each timestep
h_array = np.zeros(time_steps)
s_array = np.zeros(time_steps)
time_array = np.linspace(0,t_max,num=time_steps)
fig1 = pl.figure("h_vs_time")
#pl.figure(1)

fig2 = pl.figure("s_vs_time")
#pl.figure(2)

# Compute initial values
s_0 = s_0_calc(params.g,params.Cp,params.T,params.z)
s_init = s_init_calc(s_0, params.delta_s)

h_array[0] = h_init
s_array[0] = s_init

for i in range(1,time_steps):

    # Compute variables using functions defined above. 
    h_array[i] = h_n_plus_one_calc(alpha,params.delta_f,params.delta_s,s_array[i-1],s_0,params.D,h_array[i-1],delta_t) 
   
    s_array[i] = s_n_plus_one_calc(params.delta_f,alpha,params.V,s_array[i-1],s_0,h_array[i-1],delta_t)   

# For debugging purposes - print out the solution arrays for inspection.
np.set_printoptions(precision=3) # Controls the precision of values in the numpy arrays.
print "h_array equals: " , h_array
print "s_array equals: " , (s_array / 1000.0) # Presents values in kJ/kg
    
# Graphing commands
labelstr_alpha =str(r"$\alpha =$ ") + "%.2f" % alpha
labelstr_delta_f =str(r"$\quad \Delta F =$ ") + "%.2f" % params.delta_f 
combined_labelstr = labelstr_alpha + "  " + labelstr_delta_f 
#pl.figure(1) # Switch to figure 1.
pl.figure("h_vs_time")
    
pl.plot(time_array,h_array,label=combined_labelstr)
pl.title("Steady State Height vs Time" ,fontsize=26)
pl.xlabel(r"Time", fontsize=24)
pl.ylabel(r"$\mathrm{h} \quad \mathrm{[m]}$", fontsize=24)
pl.legend(ncol=1,loc = 'upper left')

#pl.figure(2) # Switch to figure 2. 
pl.figure("s_vs_time") 

s_array_kJ = s_array / 1000.0 # Convert values in Joules/kg to kilo Joules/kg. 
pl.plot(time_array,s_array_kJ,label=combined_labelstr)
pl.title("Steady State Dry Static Energy vs Time", fontsize=26)
pl.xlabel(r"Time", fontsize=24)
pl.ylabel(r"$\mathrm{s} \quad \mathrm{[kJ kg^{-1}]}$", fontsize=24)
pl.legend(ncol=1,loc = 'upper left')

 
#fig1.savefig("height_vs_time.pdf")

pl.show(fig1)
#pl.close(fig1)

#fig2.savefig("s_vs_time.pdf")
        
pl.show(fig2)
#pl.close(fig2)

