
# This script defines the BFP model and solves the ODEs for the model using the true parameters
# Author: Steve

##--------------------------------------------------##

# Load the required packages
using DifferentialEquations, Plots, LaTeXStrings, Distributions, Random

# Define the ODE system for the BFP model
function BFP_model(du, u, theta, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chioff = theta
    # Calculate the hill functions to output
    CHfunct = C/(C+H23)
    FHfunct = F/(F+H4)
    BHfunct = B/(B+H5)
    # Calculate the attatchment functions
    chion = (chi_on_max*interaction_strength + chi_on_min*B)/(interaction_strength + B)
    # Calculate the ODEs
    du[1] = - r2*CHfunct*F - r3*CHfunct*B
    du[2] = + e23*r2*CHfunct*F - r4*FHfunct*S - chion*F + chioff*B
    du[3] = + e23*r3*CHfunct*B - r5*BHfunct*A + chion*F - chioff*B
    du[4] = + e4*r4*FHfunct*S 
    du[5] = + e5*r5*BHfunct*A
end

### Constants

# Non-dimensionalised scaling factors
ND_T = 24.0 # Non dimension time - hours
ND_C = 1.0 # Non dimension concentration - ug/ml

# Efficiencies
e23 = 0.2 # Free bacteria growth efficiency
e4 = 0.5 # Bound bacteria growth efficiency
e5 = 0.5 # Predator growth efficiency

# Growth rates
r2 = (1/e23)*0.21*ND_T # Free bacteria growth rate - cells/hr * ND_hrs
r3 = (1/e23)*0.007*ND_T # Bound bacteria growth rate - cells/hr * ND_hrs
r4 = (1/e4)*0.12*ND_T # Predator growth rate - cells/hr * ND_hrs
r5 = (1/e5)*0.09*ND_T # Predator death rate - 1/hrs * ND_hrs

# Half saturation constants
H23 = 3.0/ND_C # Half saturations of bacteria
H4 = 3.0/ND_C # Half saturation of free predator
H5 = 0.8/ND_C # Half saturation of bound predator

# Attatchment parameters
chi_on_max = (8.0/3.0)*0.05*ND_T # Attatchment parameters
chi_on_min = (8.0/3.0)*0.0005*ND_T # Attatchment parameters
interaction_strength = (8.0/3.0)*0.01 # Interaction strength
chi_off = 0.005*ND_T

# Set the true values of the parameters
theta_true = (r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off) # Parameters

### Define and solve the ODE system

# Set the initial conditions
u0 = [ ND_C, 1.0*ND_C, 0.01, 3.0*0.01*ND_C, 8.0*0.001*ND_C ] # IVP - 2B

# Set the time span of the solution
tspan = (0.0, 1.0) # Time span of solution
sample_frequency = 0.01 # Sample frequency
output_times = tspan[1]:sample_frequency:tspan[2] # Times at which to output the solutions
N_times = length(output_times) # Number of output times
#ode_algo = AutoVern9(Rodas5()) # Optional stiff solver, required for optimisation
ode_algo = Tsit5() # Non-stiff solver

# Solve the ODE system
prob = ODEProblem(BFP_model, u0, tspan, theta_true) # Define the ODE problem
sol = solve(prob, ode_algo) # Solve the ODE problem
NofT = sol(output_times) # Extract the solution at the output times

# Define the noise level of the synthetic data
N_synthetics = 1000 # Number of synthetic data sets
synthetic_data = zeros(5, N_times, N_synthetics) # Array to store the synthetic data
μ = 1.0; σ = 0.3 # Mean and standard deviation of the noise
μ_for_mean(mu, sigma) = log( (mu^2)/sqrt(mu^2 + sigma^2) )
σ_for_sigma(mu, sigma) = log( 1 + (sigma/mu)^2 )
dist = LogNormal(μ_for_mean(μ, σ), σ_for_sigma(μ, σ)) # Log normal distribution for the noise

for i = 1:5
    for n = 1:N_synthetics
        synthetic_data[i,:,n] = NofT[i,:] .* rand(dist, N_times) # Add noise to the data
    end
end
