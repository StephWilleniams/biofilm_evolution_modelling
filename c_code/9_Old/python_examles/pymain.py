
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

## REQUIRED FUNCTIONS FOR ODE SOLUTION

# Calculate the RHS of the nutrient dynamics
def dCdt(t,Y,params):

    # Extract parameters
    r2, r3, H2, H3, chi, a, chimax, chimin, r4, r5, eff4, eff5, H4, H5, = params

    M12 = r2*Y[1]*Y[0]/(H2+Y[0]) # Nutrient loss to growth of free
    M13 = r3*Y[2]*Y[0]/(H3+Y[0]) # Nutrient loss to growth of biofilm

    return -(M12 + M13)
# Calculate the RHS of the free bacteria dynamics
def dFdt(t,Y,params):

    # Extract parameters
    r2, r3, H2, H3, chi, a, chimax, chimin, r4, r5, eff4, eff5, H4, H5, = params

    M12 = r2*Y[1]*Y[0]/(H2+Y[0]) # Nutrient loss to growth of free
    G2 = M12 # Growth of free by bacteria consumption
    M24 = r4*Y[3]*Y[1]/(H4+Y[1]) # Consumption of free bacteria by predators
    T23 = (a*chimax+chimin*Y[2])/(a+Y[2])*Y[1] # Attachment rate 
    T32 = chi*Y[2] # Dettachment rate 

    return (G2 + T32) - (M24 + T23)
# Calculate the RHS of the biofilm bacteria dynamics
def dBdt(t,Y,params):

    # Extract parameters
    r2, r3, H2, H3, chi, a, chimax, chimin, r4, r5, eff4, eff5, H4, H5, = params

    # Calculate terms
    M13 = r3*Y[2]*Y[0]/(H3+Y[0]) # Nutrient loss to growth of biofilm
    G3 = M13 # Growth of biofilm by bacteria consumption
    M35 = r5*Y[4]*Y[2]/(H5+Y[2]) # Biofilm loss to growth of predator
    T23 = (a*chimax+chimin*Y[2])/(a+Y[2])*Y[1] # Attachment rate 
    T32 = chi*Y[2] # Dettachment rate 

    # Output
    return (G3 + T23) - (M35 + T32)
# Calculate the RHS of the free predator dynamics
def dSdt(t,Y,params):

    # Extract parameters
    r2, r3, H2, H3, chi, a, chimax, chimin, r4, r5, eff4, eff5, H4, H5, = params

    # Calculate terms
    M24 = r4*Y[3]*Y[1]/(H4+Y[1]) # Consumption of free bacteria by predators
    G4 = eff4*M24 # Growth of free predator by free consumption

    # Output
    return G4
# Calculate the RHS of the biofilm predator dynamics
def dTdt(t,Y,params):

    # Extract parameters
    r2, r3, H2, H3, chi, a, chimax, chimin, r4, r5, eff4, eff5, H4, H5, = params

    M35 = r5*Y[4]*Y[2]/(H5+Y[2]) # Biofilm loss to growth of predator
    G5 = eff5*M35 # Growth of biofilm predator by biofilm consumption

    return G5
# Calculate the RHS of the full system
def dXdt(t,Y,params):
    return [dCdt(t,Y,params),dFdt(t,Y,params),dBdt(t,Y,params),dSdt(t,Y,params),dTdt(t,Y,params)]
# Solve the ODE model, single function
def solveModel(func,tspan,IC,params):
    sol = solve_ivp(lambda t,Y: func(t,Y,params), [tspan[0],tspan[-1]], IC, t_eval=tspan)
    return sol
# Calculate the log-likelihood of a given data set
def llhood(params,data,dmu,dsig,tspan,IC):
    ideal = solveModel(dXdt,tspan,IC,params) # Get model solution
    output = 0 # Store for the log-likelihood
    # Main loop
    for i in range(5): # For species
        for n in range(len(data[i,0,:])): # For synthetic data repeat
            for t in range(len(ideal.t)):
                if ideal.y[i,t] > 0:
                    resid = data[i,t,n]/ideal.y[i,t] # Calculate the residual of the data
                    # print(resid)
                    npdf = np.exp(-(np.log(resid)-dmu)**2/(2*dsig**2)) / (resid*dsig*np.sqrt(2*np.pi))
                    output += np.log(npdf) # Add this to the llh
    return -output

## PERFORM A BASIC RUN OF THE ODE SYSTEM

# Set the nominal parameters
r20 = 0.21; r30 = 0.007; r40 = 0.12; r50 = 0.09 # Growth rates
e1 = 0.5; e2 = 0.5 # Efficiency
H20 = 1; H30 = 1; H40 = 1; H50 = 0.1 # Half saturation of type-2 consumption
chi0 = 0.005; a0 = 0.01; chimax0 = 0.05; chimin0 = 0.005 # Attachment parameters
params = [r20,r30,H20,H30,chi0,a0,chimax0,chimin0,r40,r50,e1,e2,H40,H50]

# Set the initial conditions
C0=1; F0=1; B0=0.1; S0=0.01; T0=0.01 # Initial population sizes
IC = [C0,F0,B0,S0,T0] # Initial conditions for the ODE system

# Define the time parameter for the solver
tmin=0; tmax=24; Ntimes=24 # Time bounds
tspan = np.linspace(tmin,tmax,Ntimes) # Timepoints for solution

# Solve a basic ODE system as an example
sol = solveModel(dXdt,tspan,IC,params)

# Optional: plot the output
plot_basic_sol = False
if plot_basic_sol == True:
    for spec in range(5):
        plt.plot(sol.t,sol.y[spec,:])
    plt.legend(['0','1','2','3','4'])
    plt.show()

## GENERATE A SYNTHETIC DATASET

# Constants
mu = 0; sigma = 0.2 # LN parameters
N_samples = 5000 # Number of synthetic samples
modified_data = np.zeros([5,Ntimes,N_samples]) # Synth data store

# Loop through the various species to generate all the synthetic data-points
print("Now generating the synthetic data set.")
for spec in range(5):
    for n in range(N_samples):
        for s in range(Ntimes):
            modified_data[spec,s,n] = sol.y[spec, s] * np.exp(mu + sigma*np.random.normal(0,1))
print("Generated the synthetic data set.")

# Optional plot for the synthetic data
plot_synth = False
if plot_synth:
    species = 1
    index = 1
    plt.scatter(sol.t,modified_data[species,:,index])
    plt.show()

## CALCULATE THE LOGLIKELIHOOD OF THE SYNTHETIC DATASET

data = modified_data
dmu = 0; dsig = 1
# print(llhood(params,data,dmu,dsig,dXdt,tspan,IC))

params_guess = params
delta = 0.1
# Calculate bounds as ±10% of initial guesses
bounds = [(param * (1-delta), param * (1+delta)) for param in params_guess]

## OPTIMISE THIS DISTRIBUTION
print("Beginning optimisation.")
results = minimize(llhood,params_guess,args=(data,0,1,tspan,IC),bounds=bounds)
print("Optimisation complete.")

print( (results.x - params)/params )
