[println() for i in 1:10] # to make debugging easier

# Inclusions
using DifferentialEquations, Plots, Distributions, Random, Optimization, OptimizationOptimJL, OptimizationNLopt, ReverseDiff
using Base.Threads
println("Libraries Loaded.")

do_univariate = false

# Set the seed
Random.seed!(1234) # Use for repeatability

function biofilm_predators(du, u, theta, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H2, H3, H4, H5, chioff, a, chimax, chimin = theta[1:end-1]
    # Other required values
    sig = 1 # Volume/surface area. 
    chion = (a*chimax + chimin*B)/(a+B) # Bind-on rate
    # Calculate the ODEs to output
    du[1] = -r2*(C/(C+H2))*F/e23 - sig*r3*(C/(C+H3))*B/e23
    du[2] = r2*(C/(C+H2))*F - r4*(F/(F+H4))*S - chion*F*sig + chioff*B*sig
    du[3] = r3*(C/(C+H3))*B - r5*(B/(B+H5))*A + chion*F - chioff*B
    du[4] = e4*r4*(F/(F+H4))*S 
    du[5] = e5*r5*(B/(B+H5))*A
end

# Set the values required for the ODE solver
theta_true = (0.21, 0.007, 0.12, 0.09, # growth rates
              0.2, 0.5, 0.5 , # Efficiencies
              1.0, 1.0, 1.0, 0.1, # Half saturations
              0.005,0.01 ,0.05,0.005, # Attatchment parameters
              0.05) # Synthetic noise
u0 = [1.0,1.0,0.0,0.05,0.05] # IVP
tspan = (0.0, 24.0) # Time span of solution
sample_frequency = 0.1 # Time-step at which to sample
output_times = 0.0:sample_frequency:24.0 # Times at which to output the solutions
N_times = length(output_times)

# Solve the ODE system
prob = ODEProblem(biofilm_predators, u0, tspan, theta_true)
NofT = solve(prob, Tsit5(),saveat=output_times).u
println("Basic solution calculated.")

# Constants
mu = 0 # LN "mean"
sigma = theta_true[end] # LN "variance"
N_samples = 1 # Number of synthetic samples
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

# Define functions for the minimizer
function llhood(params_vary,params_const)

    if length(params_const) == 0
        params = params_vary
    else
        # Re-construct the parameter vector 
        ind = Int(params_const[2])
        params = vcat(params_vary[1:ind-1], params_const[1])
    end

    # solve the ODE problem
    prob = ODEProblem(biofilm_predators, u0, tspan, params)
    predicted_data = solve(prob, Tsit5(),saveat=output_times).u

    # Calculate the loglikelihood
    output = 0.0 # Set the llhood to output to 0)

    if params[end] < 0
        return -Inf
    else
        dist = LogNormal(0,params[end])
        for i in 1:5 # Species
            for j in 1:length(synthetic_data[1,1,:]) # Independent samples
                for k in 1:length(synthetic_data[1,:,1]) # Times
                    if predicted_data[k][i] > 0
                        e = Distributions.loglikelihood(dist,synthetic_data[i,k,j]./predicted_data[k][i]) # Get pdf of given deviation given dist
                        output += e 
                    else
                        # Do nothing
                    end
                end
            end
        end
        return -output
    end
end

# Get the optimal parameters

# lbd = 0.5 .*theta_guess
# ubd = 1.5 .*theta_guess

println("Beginning optimisation")
optf = OptimizationFunction(llhood)
prob = OptimizationProblem(optf,theta_guess,[],lb=zeros(length(theta_guess)),ub=100*ones(length(theta_guess)))
sol = solve(prob, NLopt.LN_NELDERMEAD())
theta_predicted = sol.u

for i in eachindex(theta_predicted)
    println(theta_predicted[i])
    println(theta_predicted[i]./theta_true[i]) 
    println(theta_predicted[i]./theta_guess[i])
    println()
end
println("Optimal solution found.")

#######################################################
if do_univariate

    # Loop through each of the parameters and check their practical identifiability
    N_guesses = 50 # Number of samples to check
    df=1; llstar=-quantile(Chisq(df),0.95)/2 # Get the 95% sig threshold
    for ind in 1:length(theta_true)

        println("parameter $ind")
        println(ind)

        # For all the local values quantify the log-likelihoods
        local llh_univ = zeros(N_guesses) # Store for outputs for plotting
        local guess_values = LinRange(0.9*theta_predicted[ind],1.1*theta_predicted[ind],N_guesses) # Const parameters to loop through
        local theta_temp = vcat(theta_predicted[1:ind-1],theta_predicted[ind+1:end]) # Variable param values
        local lbd = 0.5 .* theta_temp
        local ubd = 1.5 .* theta_temp
        for i = 1:N_guesses
            local optf = OptimizationFunction(llhood,AutoFiniteDiff())
            local prob2 = OptimizationProblem(optf,theta_temp,(guess_values[i],ind),lb=lbd,ub=ubd)
            local sol2 = solve(prob2, NLopt.LN_NELDERMEAD())
            llh_univ[i] = -llhood(sol2.u,(guess_values[i],ind))
        end

        println(llh_univ)

        # Plot the outputs
        plot(guess_values,llh_univ.-maximum(llh_univ),lw=4,color=:red,ylims=(-3,0.1))
        hline!([llstar],legend=false,lw=4,color=:gold)
        vline!([theta_predicted[ind]],legend=false,lw=4,color=:blue)

        # Save the outputs
        savefig("c_code/figures/PI_$ind.png")

    end
end
