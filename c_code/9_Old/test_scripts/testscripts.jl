# Title: Simple parameter identifiablity
# Author: Stephen Williams

# Inclusions
using DifferentialEquations, Plots, Distributions, Random, Optimization, OptimizationOptimJL

# Set the seed
Random.seed!(1234) # for repeatability

# Define required functions
function growth_equation(du,u,p,t) # Logistic growth
    du[1] = p[1] * u[1] * (1 - u[1]) - p[2] * u[1]
end

# Set parameters for the solver to create synthetic data
theta = [1.0, 0.5, 0.01]
u0 = [0.1]
tspan = (0.0, 20.0) # Time interval
N_times = 200
output_times = LinRange(tspan[1],tspan[2],N_times)
prob = ODEProblem(growth_equation, u0, tspan, theta)
sol = solve(prob, Tsit5()); data = sol(output_times)[1,:]
mu = 0
sigma_true = theta[end]
N_samples = 1
synthetic_data = zeros(N_samples,length(data)) # Create a store for the synthetic data
[synthetic_data[i,:] = [ data_point* exp.(mu.+sigma_true*randn()) for data_point in data] for i in 1:N_samples]

# Get some random guesses for the params
theta_guess = theta # Copy the baseline(/true) parameters 
modification_scaling = 0 # Percent change range for uniform scaling of true parameters
scalings = (1 .+ modification_scaling*2*(rand(length(theta)).-0.5))
theta_guess = theta_guess .* scalings # Initial guess for parameters

# Define functions for the minimizer
function llhood2!(params)

    # Simplifying hardcoding
    u0 = [0.1] 
    output_times = LinRange(0,20,N_times)
    tspan = (output_times[1], output_times[end])
    prob = ODEProblem(growth_equation, u0, tspan, params)
    predicted_ode = solve(prob, Tsit5()); 
    predicted_data = predicted_ode(output_times)[1,:]

    # Calculate the loglikelihood
    output = 0 # Set the llhood to output to 0

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

optf = OptimizationFunction(llhood2!)
prob2 = OptimizationProblem(optf,theta_guess,[])
sol = solve(prob2, Optim.NelderMead())

# Reconstruct the data
u0 = [0.1] 
output_times = LinRange(0,20,N_times)
tspan = (output_times[1], output_times[end])
prob = ODEProblem(growth_equation, u0, tspan, sol.u)
pred = solve(prob, Tsit5())
pred_vals = pred(output_times)[1,:]

plot(output_times,data)
plot!(output_times,synthetic_data[1,:])
plot!(output_times,pred_vals)
