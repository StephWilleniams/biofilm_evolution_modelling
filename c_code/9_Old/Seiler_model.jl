
# This script defines the BFP model and solves the ODEs for the model using the true parameters
# Author: Steve

##--------------------------------------------------##

# Load the required packages
using DifferentialEquations, Plots

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
    # sigma
    sigma = 8.0/3.0
    # Calculate the ODEs
    du[1] = - (1/e23)*(r2*CHfunct*F - r3*CHfunct*B*sigma)
    du[2] = + r2*CHfunct*F - r4*FHfunct*S - chion*F*sigma + chioff*B*sigma
    du[3] = + r3*CHfunct*B - r5*BHfunct*A + chion*F - chioff*B
    du[4] = + e4*r4*FHfunct*S 
    du[5] = + e5*r5*BHfunct*A
end

### Constants

# Non-dimensionalised scaling factors
ND_T = 1.0 # Non dimension time - hours
ND_C = 1.0 # Non dimension concentration - ug/ml

# Efficiencies
e23 = 0.2 # Free bacteria growth efficiency
e4 = 0.5 # Bound bacteria growth efficiency
e5 = 0.5 # Predator growth efficiency

# Growth rates
r2 = 0.21*ND_T # Free bacteria growth rate - cells/hr * ND_hrs
r3 = 0.007*ND_T # Bound bacteria growth rate - cells/hr * ND_hrs
r4 = 0.12*ND_T # Predator growth rate - cells/hr * ND_hrs
r5 = 0.09*ND_T # Predator death rate - 1/hrs * ND_hrs

# Half saturation constants
H23 = 1.0 # Half saturations of bacteria
H4 = 1.0 # Half saturation of free predator
H5 = 0.1 # Half saturation of bound predator

# Attatchment parameters
chi_on_max = 0.05*ND_T # Attatchment parameters
chi_on_min = 0.0005*ND_T # Attatchment parameters
interaction_strength = 0.01 # Interaction strength
chi_off = 0.005*ND_T

# Set the true values of the parameters
theta_true = ( r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off ) # Parameters

### Define and solve the ODE system

# Set the initial conditions
u0 = [ 1.0, 1.0, 0.0, 0.01, 0.001 ] # IVP

# Set the time span of the solution
tspan = (0.0, 24.0) # Time span of solution
sample_frequency = 0.01 # Sample frequency
output_times = tspan[1]:sample_frequency:tspan[2] # Times at which to output the solutions
N_times = length(output_times) # Number of output times
#ode_algo = AutoVern9(Rodas5()) # Optional stiff solver, required for optimisation
ode_algo = Tsit5() # Non-stiff solver

# Solve the ODE system
prob = ODEProblem(BFP_model, u0, tspan, theta_true) # Define the ODE problem
sol = solve(prob, ode_algo) # Solve the ODE problem
NofT = sol(output_times) # Extract the solution at the output times

### Plot the solution

plot_solver = true

if plot_solver

    #Scale = [3.0,3.0,8.0,3.0,8.0];

    plt = plot()
    for i = 1:5
        plot!(plt,output_times, NofT[i,:], label="Species $i")
    end 
    display(plt)

end

##--------------------------------------------------##

function model(theta)

    # Fixed values needed for the solver to work
    u0 = [1.0,1.0,0.0,0.05,0.01] # IVP
    output_times = LinRange(0,24,25) # Solver times
    tspan = (output_times[1],output_times[end]) # End-points for the solver time

    # Run the solver for the system
    prob = ODEProblem(biofilm_predators, u0, tspan, theta) # Define the ODE problem
    sol = solve(prob, AutoVern9(Rodas5()), saveat=output_times) # Solve the ODE problem

    # Set the parameter of interest for our outputs to check sensitivity with respect to
    QOI = sol[2,end]

    # Return the quantity of interest
    return QOI

end

