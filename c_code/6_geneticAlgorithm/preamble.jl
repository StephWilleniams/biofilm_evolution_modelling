
# Import modules
include("GeneticAlgorithm.jl")
import .GeneticAlgorithm

using DifferentialEquations, Distributions, Plots

# Define constants
t_start = 0.0 # Sample schedule start time
t_end = 1.0 # Sample schedule maximum time
dt = 0.01 # Sample schedule minimum time step
n_times = 20 # Number of time points in the schedule
N_synthetics = 30 # Number of synthetic data sets (experimental repeat analog) to use for the profile likelihood optimisation
N_organisms = 100 # Number of organisms (schedules) in the population for selection
N_generations = 2 # Number of generations to run the genetic algorithm
cloning_noise = 0.1 # Cloning mutation noise level
index_of_interest = 11 # Index of the parameter of interest
PL_range = 0.3
n_PLs = 31 # Number of points in the profile likelihood optimisation
initial_varying_params = [0.05,0.05,0.05]
varying_indices = [12,13,14]

# Preallocate organisms
global organisms = Vector{GeneticAlgorithm.organism}(undef,N_organisms)

# Define the model
ode_algo = Tsit5()
ND_T = 24.0 # Non-dimension time - hours
ND_C = 1.0 # Non-dimension concentration - ug/ml
e23 = 0.2 # Free bacteria growth efficiency
e4 = 0.5 # Bound bacteria growth efficiency
e5 = 0.5 # Predator growth efficiency
r2 = (1/e23)*0.21*ND_T # Free bacteria growth rate - cells/hr * ND_hrs
r3 = (1/e23)*0.007*ND_T # Bound bacteria growth rate - cells/hr * ND_hrs
r4 = (1/e4)*0.12*ND_T # Predator growth rate - cells/hr * ND_hrs
r5 = (1/e5)*0.09*ND_T # Predator death rate - 1/hrs * ND_hrs
H23 = 3.0/ND_C # Half saturations of bacteria
H4  = 3.0/ND_C # Half saturation of free predator
H5  = 0.8/ND_C # Half saturation of bound predator
chi_on_max = (8.0/3.0)*0.05*ND_T # Attatchment parameters
chi_on_min = (8.0/3.0)*0.0005*ND_T # Attatchment parameters
interaction_strength = (8.0/3.0)*0.01 # Interaction strength
chi_off = 0.005*ND_T
theta_baseline = [r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off] # Parameters
u0 = [ ND_C, 1.0*ND_C, 0.01, 3.0*0.01*ND_C, 8.0*0.001*ND_C ] # IVP

# Define synthetic noise
μ = 1.0 # Mean of the noise in the measurements/synthetic data
σ = 0.1 # Standard deviation of the noise in the measurements/synthetic data
dist = LogNormal(GeneticAlgorithm.μ_for_mu(μ, σ), GeneticAlgorithm.σ_for_sigma(μ, σ))
