# Title: Simple parameter identifiablity
# Author: Stephen Williams

"""
Notes: 
The script calculates a simple logistic growth curve.
The values of this "true" dataset are then corrupted with some known noise.
We assume there is an unknown, but constant, noise level in the data sigma.
For a given parameter set (which contains the noise level) to log-likelihood.
By finding the set which maximise (minimising -llh) the llh we can find the optimal parameters.
Following this, we can look at the profile likelihood of a univariate model about the nominal.
This gives us an idea, for a given level of (synthetic) data, of the identifiability.
Outputs "figures/PI_*.png" showing the profile likelihoods.
"""

# Inclusions
using DifferentialEquations, Plots, Distributions, Random, Optimization, OptimizationOptimJL

# Set the seed
# Random.seed!(1234) # Use for repeatability

# Define required functions
function growth_equation(du,u,p,t) # Logistic growth
    du[1] = p[1] * u[1] * (1 - u[1]) - p[2] * u[1]
end

# Set parameters for the solver to create synthetic data
theta = [1.0, 0.5, 0.01] # True parameterss
u0 = [0.1] # Initial condition
tspan = (0.0, 20.0) # Time interval
N_times = 10 # Number of sample times
output_times = LinRange(tspan[1],tspan[2],N_times) # Sample times for the data 

# Get the baseline "true" solution
prob = ODEProblem(growth_equation, u0, tspan, theta)
sol = solve(prob, Tsit5()); data = sol(output_times)[1,:]

# Create the synthetic data
mu = 0
sigma_true = theta[end]
N_samples = 5
synthetic_data = zeros(N_samples,length(data)) # Create a store for the synthetic data
[synthetic_data[i,:] = [ data_point* exp.(mu.+sigma_true*randn()) for data_point in data] for i in 1:N_samples]

# Get some random guesses for the params
theta_guess = theta # Copy the baseline(/true) parameters 
modification_scaling = 0 # Percent change range for uniform scaling of true parameters
scalings = (1 .+ modification_scaling*2*(rand(length(theta)).-0.5))
theta_guess = theta_guess .* scalings # Initial guess for parameters

# Define functions for the minimizer
function llhood2!(params_vary,params_const)

    # Re-construct the parameter vector 
    ind = Int(params_const[2])
    params = vcat(params_vary[1:ind-1], params_const[1])
    params = vcat(params,params_vary[ind:end])

    # Simplifying hardcoding
    u0 = [0.1] 
    output_times = LinRange(0,20,N_times)
    tspan = (output_times[1], output_times[end])
    prob = ODEProblem(growth_equation, u0, tspan, params)
    predicted_ode = solve(prob, Tsit5()); 
    predicted_data = predicted_ode(output_times)[1,:]

    # Calculate the loglikelihood
    output = 0 # Set the llhood to output to 0

    # Make sure the noise in the distribtuion is physical
    if params[end] < 0
        return Inf
    else
        dist = LogNormal(0,params[end])
        for i in 1:length(synthetic_data[:,1]) # loop through the independent samples
            for n in 1:length(output_times) # loop through the times
                e = Distributions.loglikelihood(dist,synthetic_data[i,n]./predicted_data[n]) # Get pdf of given deviation given dist
                output = output + e 
            end
        end

        return -output
    end
    
end

# Loop through each of the parameters and check their practical identifiability
N_guesses = 100 # Number of samples to check
df=1; llstar=-quantile(Chisq(df),0.95)/2 # Get the 95% sig threshold
for ind in 1:3

    # For all the local values quantify the log-likelihoods
    local llh_univ = zeros(N_guesses) # Store for outputs for plotting
    local guess_values = LinRange(0.9*theta_guess[ind],1.1*theta_guess[ind],N_guesses) # Const parameters to loop through
    local theta_temp = vcat(theta_guess[1:ind-1],theta_guess[ind+1:end]) # Variable param values
    for i = 1:N_guesses
        local optf = OptimizationFunction(llhood2!)
        local prob2 = OptimizationProblem(optf,theta_temp,[guess_values[i],ind])
        local sol = solve(prob2, Optim.NelderMead())
        llh_univ[i] = -llhood2!(sol.u,[guess_values[i],ind])
    end

    # Plot the outputs
    plot(guess_values,llh_univ .- maximum(llh_univ),lw=4,color=:red)
    hline!([llstar],legend=false,lw=4,color=:gold)
    vline!([theta_guess[ind]],legend=false,lw=4,color=:blue)

    # Save the outputs
    savefig("c_code/figures/PI_$ind.png")

end
