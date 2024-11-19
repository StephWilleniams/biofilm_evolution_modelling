[println() for i in 1:10] # to make debugging easier
println("Now running reduced complexity model - Univariate parameter identifiability.")

# REDUCED COMPLEXITY BACTERIA ONLY MODEL

# Inclusions
using DifferentialEquations, Plots, Distributions, Random, NLopt#, ReverseDiff
using Base.Threads
println("Libraries Loaded.")

# Set the seed
Random.seed!(1234) # Use for repeatability

function biofilm_predators(du, u, theta, t)
    # Extract the different species
    C, F, B = u
    # Extract parameters
    r2, r3, H23, chioff, chion = theta
    # Calculate the ODEs to output
    du[1] = -r2*(C/(C+H23))*F/0.2 - r3*(C/(C+H23))*B/0.2
    du[2] = r2*(C/(C+H23))*F - chion*F + chioff*B
    du[3] = r3*(C/(C+H23))*B + chion*F - chioff*B
end

# Set the values required for the ODE solver
theta_true = (0.21, 0.007, # growth rates
              1.0, # Half saturations
              0.005,0.0005) # Attatchment parameters
u0 = [3.0,1.0,0.1] # IVP
tspan = (0.0, 24.0) # Time span of solution
sample_frequency = 0.05 # Time-step at which to sample
output_times = tspan[1]:sample_frequency:tspan[2] # Times at which to output the solutions
N_times = length(output_times)

# Solve the ODE system
prob = ODEProblem(biofilm_predators, u0, tspan, theta_true)
sol = solve(prob, Tsit5())
NofT = sol(output_times)
println("Basic solution calculated.")

# Generate synthetic data
N_samples = 5 # Number of synthetic samples
synthetic_data = zeros(5,N_times,N_samples)
μ_for_mean(mu, sigma) = log(mu) - sigma^2/2
epsilon = 1e-8
mu = 1; sigma = 0.05
dist = LogNormal(μ_for_mean(mu, sigma), sigma)
for type_species in 1:3
    for n in 1:N_samples
        for s in 1:N_times
            synthetic_data[type_species,s,n] = NofT[type_species,s] * rand(dist)
        end
    end
end
println("Synthetic solution calculated.")

plt = plot()
for i in 1:3
    for j = 1:N_samples
        rels = synthetic_data[i,:,j]./NofT[i,:]
        # scatter!(plt,1:N_times,rels)
    end
end
display(plt)

# # Get the params on which to optimise
# modification_scaling = 0.1 # Percent change range for uniform scaling of true parameters
# scalings = (1 .+ modification_scaling*2*(rand(length(theta_true)).-0.5))
# theta_guess = theta_true .* scalings # Initial guess for parameters
# println("Parameter guess calculated")

# # Define functions for the minimizer
# function llhood1(params,_)

#     # solve the ODE problem
#     prob = ODEProblem(biofilm_predators, u0, tspan, params)
#     predicted_sol = solve(prob, Tsit5())
#     predicted_data = predicted_sol(output_times)
#     output = 0.0
#     for i in 1:3 # Species
#         for k in 1:N_times # Independent samples
#             rels = [max(synthetic_data[i,k,j]./predicted_data[i,k],epsilon) for j in 1:N_samples]
            
#             matches = findall(x -> x <= 1e-5, rels)
#             if !isempty(matches)
#                 println(matches)
#             end

#             rels[matches].=epsilon

#             val = Distributions.logpdf(dist,rels)

#             matches = findall(x -> x == Inf, val)
#             if !isempty(matches)
#                 val[matches] .= 0
#             end

#             output += sum(val)
#         end
#     end
#     return output
# end

# opts = NLopt.Opt(:LN_NELDERMEAD, length(theta_guess))
# NLopt.lower_bounds!(opts,1e-6)
# NLopt.upper_bounds!(opts,3)
# NLopt.max_objective!(opts, llhood1)

# max_f_opt, min_x_opt, ret = NLopt.optimize(opts, theta_guess)
# max_f_opt, min_x_opt, ret = NLopt.optimize(opts, min_x_opt)
# max_f_opt, x_opt, ret = NLopt.optimize(opts, min_x_opt)

# println(x_opt)
# println(max_f_opt)

# do_univariate = false
# if do_univariate == true

#     N_varies = 101
#     modifier = LinRange(0.9,1.1,N_varies)
#     for ii in 1:5
#         out = zeros(N_varies)
#         for jj in eachindex(out)

#             params_const = [x_opt[ii]*modifier[jj],ii]
#             params_vary = vcat(x_opt[1:ii-1],x_opt[ii+1:end])

#             # Define functions for the minimizer
#             function llhood(params_vary1,_)

#                 # Re-construct the parameter vector 
#                 ind = Int(params_const[2])
#                 params = vcat(params_vary1[1:ind-1], params_const[1])
#                 params = vcat(params, params_vary1[ind:end])

#                 # solve the ODE problem
#                 prob = ODEProblem(biofilm_predators, u0, tspan, params)
#                 predicted_sol = solve(prob, Tsit5())
#                 predicted_data = predicted_sol(output_times)

#                 # Calculate the loglikelihood
#                 output = 0.0
#                 for i in 1:3 # Species
#                     for k in 1:N_times # Independent samples
#                         rels = [max(synthetic_data[i,k,j]./predicted_data[i,k],epsilon) for j in 1:N_samples]
#                         val = Distributions.logpdf(dist,rels)
#                         output += sum(val)
#                     end
#                 end
#                 return output
#                 # epsilon = 1e-12
#                 # output=0.0
#                 # for i in 1:3 # Species
#                 #     for k in 1:N_times # Independent samples
#                 #         #rels = [max(synthetic_data[i,k,j]./predicted_data[i,k],epsilon) for j in 1:N_samples]
#                 #         rels = synthetic_data[i,k,:]./predicted_data[i,k]
#                 #         val = Distributions.loglikelihood(dist,rels)
#                 #         #output += max(val,epsilon)
#                 #         output += val
#                 #     end
#                 # end
#                 # return output
#             end

#             local opts = NLopt.Opt(:LN_NELDERMEAD, length(params_vary))
#             NLopt.lower_bounds!(opts,1e-6)
#             NLopt.upper_bounds!(opts,3)
#             NLopt.max_objective!(opts, llhood)

#             local max_f, max_x, ret = NLopt.optimize(opts, params_vary)
#             out[jj] = max_f

#         end
#         println("Branch $ii calculated")
#         # println(out)
#         plot(modifier,out .- max_f_opt,ylims=(-10,1))
#         df=1; llstar=-quantile(Chisq(df),0.95)/2 # Get the 95% sig threshold
#         hline!([llstar],legend=false,lw=4,color=:gold)
#         savefig("PI_basicTest_$ii.png")
#     end

# end
