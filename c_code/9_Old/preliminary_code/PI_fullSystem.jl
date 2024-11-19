[println() for i in 1:10] # to make debugging easier
println("Now running full complexity model - Univariate parameter identifiability.")

# REDUCED COMPLEXITY BACTERIA ONLY MODEL

# Inclusions
using DifferentialEquations, Plots, Distributions, Random, NLopt#, ReverseDiff
using Base.Threads
println("Libraries Loaded.")

# PI factors
N_samples = 1 # Number of synthetic samples
sample_frequency = 0.005 # Time-step at which to sample
mu = 1; sigma = 0.1 # Noise parameters
modification_scaling = 0.5 # Percent change range for uniform scaling of true parameters
N_varies = 31 # Univariate PI range bins
uniRange = 0.5 # Univariate PI range for modifier
ode_algo = AutoVern9(Rodas5())
do_univariate = true

# Set the seed
Random.seed!(1234) # Use for repeatability

function biofilm_predators(du, u, theta, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H2, H3, H4, H5, chioff, a, chimax, chimin = theta
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
theta_true = (0.21,0.007,0.12,0.09, # growth rates
         0.2,0.5,0.5, # Efficiencies
         1.0,1.0,1.0,0.1, # Half saturations
         0.005,0.01,0.05,0.0005) # Attatchment parameters

u0 = [3.0,0.001,0.1,0.01,0.01] # IVP
tspan = (0.0, 24.0) # Time span of solution
output_times = tspan[1]:sample_frequency:tspan[2] # Times at which to output the solutions
N_times = length(output_times)

# Solve the ODE system
prob = ODEProblem(biofilm_predators, u0, tspan, theta_true)
sol = solve(prob, ode_algo)
NofT = sol(output_times)
println("Basic solution calculated.")

# Generate synthetic data
synthetic_data = zeros(5,N_times,N_samples)
μ_for_mean(mu, sigma) = log(mu) - sigma^2/2
dist = LogNormal(μ_for_mean(mu, sigma), sigma)
for type_species in 1:3
    for n in 1:N_samples
        for s in 1:N_times
            synthetic_data[type_species,s,n] = NofT[type_species,s] * rand(dist)
        end
    end
end
println("Synthetic solution calculated.")

# Get the params on which to optimise
scalings = (1 .+ modification_scaling*2*(rand(length(theta_true)).-0.5))
theta_guess = theta_true .* scalings # Initial guess for parameters
println("Parameter guess calculated")

# Define functions for the minimizer
function llhood1(params,_)
    # solve the ODE problem
    prob = ODEProblem(biofilm_predators, u0, tspan, params)
    predicted_sol = solve(prob, ode_algo)
    predicted_data = predicted_sol(output_times)

    output = 0.0
    for i in 1:3
        for j = 1:N_samples
            rels = synthetic_data[i,:,j]./predicted_data[i,:]
            output += Distributions.loglikelihood(dist,rels)
        end
    end
    return output
end

opts = NLopt.Opt(:LN_NELDERMEAD, length(theta_guess))
NLopt.lower_bounds!(opts,1e-6);NLopt.upper_bounds!(opts,3)
NLopt.xtol_rel!(opts,1e-14)
NLopt.max_objective!(opts, llhood1)

f_opt, x_opt, ret = NLopt.optimize(opts, theta_guess)

println(x_opt)
println(f_opt)

if do_univariate

    modifier = LinRange(1-uniRange,1+uniRange,N_varies)
    for ii in 1:15 # Loop on parameters

        out = zeros(N_varies)

        for jj in eachindex(out) # Loop on values of fixed param
            # Preconstruct problem parameter global values
            params_const = [x_opt[ii]*modifier[jj],ii]
            params_vary = vcat(x_opt[1:ii-1],x_opt[ii+1:end])

            # Define functions for the minimizer
            function llhood(params_vary1,_)
                #Re-construct the parameter vector 
                ind = Int(params_const[2])
                params = vcat(params_vary1[1:ind-1], params_const[1])
                params = vcat(params, params_vary1[ind:end])

                # solve the ODE problem
                prob = ODEProblem(biofilm_predators, u0, tspan, params)
                predicted_sol = solve(prob, ode_algo)
                predicted_data = predicted_sol(output_times)

                output = 0.0
                for i in 1:3
                    for j = 1:N_samples
                        rels = synthetic_data[i,:,j]./predicted_data[i,:]
                        output += Distributions.loglikelihood(dist,rels)
                    end
                end
                return output
            end

            local opts = NLopt.Opt(:LN_NELDERMEAD, length(params_vary))
            NLopt.lower_bounds!(opts,1e-6);NLopt.upper_bounds!(opts,3)
            NLopt.xtol_rel!(opts,1e-12)
            NLopt.max_objective!(opts, llhood)

            local max_f, max_x, ret = NLopt.optimize(opts, params_vary)
            out[jj] = max_f

        end
        println("Branch $ii calculated")

        plot(modifier,out .- maximum(out),ylims=(-3,1))
        df=1; llstar=-quantile(Chisq(df),0.95)/2 # Get the 95% sig threshold
        hline!([llstar],legend=false,lw=4,color=:gold)
        savefig("PI_basicTest_$ii.png")
        
    end

end
