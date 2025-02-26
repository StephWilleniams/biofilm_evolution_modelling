
# Imports ############################
using DifferentialEquations, Plots, Distributions, Random, NLopt

# Functions ############################

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

function sample_model(model, u0, output_times, theta_true, ode_algo)

    tspan = (output_times[1], output_times[end]) # Time span of solution
    prob = ODEProblem(model, u0, tspan, theta_true) # Define the ODE problem
    sol = solve(prob, ode_algo) # Solve the ODE problem
    NofT = sol(output_times) # Extract the solution at the output times
    return NofT
end

function generate_synthetic_data(dist, theta, u0, N_samples, sample_times, ode_algo)
    NofT = sample_model(BFP_model, u0, sample_times, theta, ode_algo)
    N_times = length(sample_times)
    synthetic_data = zeros(5,N_times,N_samples)
    for i = 1:5
        for n = 1:N_samples
            synthetic_data[i,:,n] = NofT[i,:] .* rand(dist, N_times) # Add noise to the data
            synthetic_data[i,1,n] = NofT[i,1] # Set the initial condition exactly
        end
    end
    return synthetic_data
end

μ_for_mean(mu, sigma) = log( (mu^2)/sqrt(mu^2 + sigma^2) )
σ_for_sigma(mu, sigma) = log( 1 + (sigma/mu)^2 )

function estimate_parameters(model)

end

# Set the sample size and the time range ############################

N_samples = 3
N_times = 10
sample_times = LinRange(0,1,N_times)

# parameters ############################

ode_algo = Tsit5()
u0 = [1.0, 1.0, 0.1, 0.01, 0.01]

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
theta = (r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off) # Parameters

# parameter estimation ############################

Random.seed!(1)
vary_index = [false for i in eachindex(theta)]
params_to_vary = [1,3,5]
[vary_index[i] = true for i in params_to_vary]
modification_scaling = 0.1

theta_guess = zeros(length(theta))
for i in eachindex(vary_index)
    if vary_index[i]
        gamma = (1 + modification_scaling*2*(rand()-0.5))
        theta_guess[i] = gamma * theta[i]
    else
        theta_guess[i] = theta[i]
    end
end

# Generate Synthetic Data ############################

μ = 1.0 # Mean of noise distribution
σ = 0.1 # Standard deviation of noise distribution
dist = LogNormal(μ_for_mean(μ, σ), σ_for_sigma(μ, σ)) # Log-normal distribution for noise
synthetic_data = generate_synthetic_data(dist, theta, u0, N_samples, sample_times, ode_algo)

# plotcheck = false
# if plotcheck 
#     plot_sample_index = 1
#     plt = plot()
#     for i = 1:5
#         plot!(plt,sample_times, synthetic_data[i,:,plot_sample_index], label="Species $i")
#     end
#     xlabel!("Time")
#     ylabel!("Concentration")
#     title!("Synthetic Data")
#     display(Plots.plot!())
# end

# Estimate the parameters ############################

# Define functions for the minimizer
function llhood(params,_) 

    full_parameters = zeros(length(theta))  
    [full_parameters[i] = theta[i] for i in eachindex(theta)]
    [full_parameters[params_to_vary[i]] = params[i] for i in 1:3]

    NofT = sample_model(BFP_model, u0, sample_times, full_parameters, ode_algo)
    output = 0.0
    for i in 1:5
        for j = 1:N_samples
            rels = synthetic_data[i,:,j]./NofT[i,:]
            output += Distributions.loglikelihood(dist,rels)
        end
    end
    return output
end

opts = NLopt.Opt(:LN_NELDERMEAD, length(params_to_vary))
#NLopt.xtol_rel!(opts,1e-12)
NLopt.max_objective!(opts, llhood)

start_points = [theta_guess[i] for i in params_to_vary]
f_opt, x_opt, ret = NLopt.optimize(opts, start_points) 

values = 100 .* (1 .- (x_opt./theta[params_to_vary]))
println("Percent difference from true value: $values ")
