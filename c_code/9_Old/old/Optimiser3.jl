[println() for i in 1:10] # to make debugging easier

# SLIGHTLY REDUCED COMPLEXITY MODEL - STILL OVERPARAM

# Inclusions
using DifferentialEquations, Plots, Distributions, Random, NLopt#, ReverseDiff
using Base.Threads
println("Libraries Loaded.")

do_univariate = false

# Set the seed
Random.seed!(1234) # Use for repeatability

function biofilm_predators(du, u, theta, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    r2, r3, r4, r5, e4, e5, H23, H4, H5, chioff, chion = theta[1:end-1]
    # Calculate the ODEs to output
    du[1] = -r2*(C/(C+H23))*F/0.2 - r3*(C/(C+H23))*B/0.2
    du[2] = r2*(C/(C+H23))*F - r4*(F/(F+H4))*S - chion*F + chioff*B
    du[3] = r3*(C/(C+H23))*B - r5*(B/(B+H5))*A + chion*F - chioff*B
    du[4] = e4*r4*(F/(F+H4))*S 
    du[5] = e5*r5*(B/(B+H5))*A
end

# Set the values required for the ODE solver
theta_true = (0.21, 0.007, 0.12, 0.09, # growth rates
              0.5, 0.5 , # Efficiencies
              1.0, 1.0, 0.1, # Half saturations
              0.005,0.0005, # Attatchment parameters
              0.1) # Synthetic noise
u0 = [1.0,1.0,0.1,0,0] # IVP
tspan = (0.0, 24.0) # Time span of solution
sample_frequency = 0.1 # Time-step at which to sample
output_times = 0.0:sample_frequency:24.0 # Times at which to output the solutions
N_times = length(output_times)

# Solve the ODE system
prob = ODEProblem(biofilm_predators, u0, tspan, theta_true)
sol = solve(prob, Tsit5())
NofT = sol(output_times)
println("Basic solution calculated.")

# Constants
mu = 0 # LN "mean"
sigma = theta_true[end] # LN "variance"
N_samples = 10 # Number of synthetic samples
synthetic_data = zeros(5,N_times,N_samples)

# Loop through the various species to generate synthetic data
for type_species in 1:5
    for n in 1:N_samples
        for s in 1:N_times
            synthetic_data[type_species,s,n] = NofT[s][type_species] * exp.(mu.+sigma*randn())
        end
    end
end
println("Synthetic solution calculated.")

# Get the params on which to optimise
modification_scaling = 0.05 # Percent change range for uniform scaling of true parameters
scalings = (1 .+ modification_scaling*2*(rand(length(theta_true)).-0.5))
theta_guess = theta_true .* scalings # Initial guess for parameters
println("Parameter guess calculated")

# out = zeros(length(theta_guess))
# @sync @threads for ii in eachindex(theta_guess)
ii = 1
params_const = [theta_guess[ii]*1,ii]
params_vary = vcat(theta_guess[1:ii-1],theta_guess[ii+1:end])

# Define functions for the minimizer
function llhood(params_vary1,_)

    # Re-construct the parameter vector 
    ind = Int(params_const[2])
    params = vcat(params_vary1[1:ind-1], params_const[1])
    params = vcat(params, params_vary1[ind:end])

    # solve the ODE problem
    prob = ODEProblem(biofilm_predators, u0, tspan, params)
    predicted_sol = solve(prob, Tsit5())
    predicted_data = predicted_sol(output_times)

    # Calculate the loglikelihood
    output = BigFloat(0.0); e = BigFloat(0.0) # Set the llhood to output to 0)
    dist = LogNormal(0,params[end])

    for i in 2:5 # Species
        for j in 1:length(synthetic_data[1,1,:]) # Independent samples
            for k in 1:length(synthetic_data[1,:,1]) # Times
                if predicted_data[i,k] > 0
                    e = Distributions.loglikelihood(dist,synthetic_data[i,k,j]./predicted_data[i,k]) # Get pdf of given deviation given dist
                    #e = max(e,0)
                    output += e
                else
                    # Do nothing
                end
            end
        end
    end
    return -output
end

opts = NLopt.Opt(:LN_NELDERMEAD, length(params_vary))
NLopt.lower_bounds!(opts,1e-4)
NLopt.upper_bounds!(opts,5)
NLopt.xtol_rel!(opts,1e-8)
NLopt.min_objective!(opts, llhood)

min_f, min_x, ret = NLopt.optimize(opts, params_vary)

# STORE OUTPUTS
println(min_x./theta_true[2:end])
println(-llhood(min_x,[]))
