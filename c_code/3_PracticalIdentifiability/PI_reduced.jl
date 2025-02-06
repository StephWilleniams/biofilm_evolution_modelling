
# # Inclusions
# using Base.Threads
# println("Libraries Loaded.")
using DifferentialEquations, Plots, Distributions, Random, NLopt

# PI factors
N_samples = 50 # Number of synthetic samples
sample_frequency = 0.01 # Time-step at which to sample
μ = 1.0; σ = 0.1 # Mean and standard deviation of the noise
μ_for_mean(mu, sigma) = log((mu^2)/sqrt(mu^2 + sigma^2))
σ_for_sigma(mu, sigma) = log(1+(sigma/mu)^2)

modification_scaling = 0.1 # Percent change range for uniform scaling of true parameters
N_varies = 101 # Univariate PI range bins
uniRange = 0.2 # Univariate PI range for modifier
#ode_algo = AutoVern9(Rodas5())
ode_algo = Tsit5() # Non-stiff solver
do_univariate = true

# Set the seed
# Random.seed!(1234) # Use for repeatability

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
H4  = 3.0/ND_C # Half saturation of free predator
H5  = 0.8/ND_C # Half saturation of bound predator

# Attatchment parameters
chi_on_max = (8.0/3.0)*0.05*ND_T # Attatchment parameters
chi_on_min = (8.0/3.0)*0.0005*ND_T # Attatchment parameters
interaction_strength = (8.0/3.0)*0.01 # Interaction strength
chi_off = 0.005*ND_T

# Set the true values of the parameters
theta_true = (r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off) # Parameters

# Set the initial conditions
u0 = [ ND_C, 1.0*ND_C, 0.01, 3.0*0.01*ND_C, 8.0*0.001*ND_C ] # IVP - 2B

# Set the time span of the solution
tspan = (0.0, 1.0) # Time span of solution
sample_frequency = 0.01 # Sample frequency
output_times = tspan[1]:sample_frequency:tspan[2] # Times at which to output the solutions
N_times = length(output_times) # Number of output times

# Get the parameter guess to use as a startpoint to optimise
scalings = (1 .+ modification_scaling*2*(rand(length(theta_true)).-0.5))
theta_guess = theta_true .* scalings # Initial guess for parameters

dist = LogNormal(μ_for_mean(μ, σ), σ_for_sigma(μ, σ))

# Define functions for the minimizer
function llhood1(params,_) 

    params = [theta_true[1:10]..., params[1:4]...]

    prob = ODEProblem(BFP_model, u0, tspan, params)
    predicted_sol = solve(prob, ode_algo)
    NofT = predicted_sol(output_times)
    output = 0.0
    for i in 1:5
        for j = 1:N_samples
            rels = synthetic_data[i,:,j]./NofT[i,:]
            output += Distributions.loglikelihood(dist,rels)
        end
    end
    return output
end

opts = NLopt.Opt(:LN_NELDERMEAD, length(theta_guess[11:end]))
#NLopt.xtol_rel!(opts,1e-12)
NLopt.max_objective!(opts, llhood1)

f_opt, x_opt, ret = NLopt.optimize(opts, theta_guess[11:end])

[println(x_opt[i-10]/theta_true[i]) for i in 11:14]

do_univariate = true
if do_univariate
    modifier = LinRange(1-uniRange,1+uniRange,N_varies)
    for ii in 1:4 # Loop on parameters
        out = zeros(N_varies)
        for jj in eachindex(out) # Loop on values of fixed param
            # Preconstruct problem parameter global values
            params_const = [x_opt[ii]*modifier[jj],ii]
            params_vary = vcat(x_opt[1:ii-1],x_opt[ii+1:end])
            # Define functions for the minimizer
            function llhood(params_vary1,_)

                #Re-construct the parameter vector 
                ind = Int(params_const[2])
                params = [theta_true[1:10]...,params_vary1[1:ind-1]..., params_const[1]..., params_vary1[ind:end]...]

                # solve the ODE problem
                prob = ODEProblem(BFP_model, u0, tspan, params)
                predicted_sol = solve(prob, ode_algo)
                predicted_data = predicted_sol(output_times)
                output = 0.0
                for i in 1:5
                    for j = 1:N_samples
                        rels = synthetic_data[i,:,j]./predicted_data[i,:]
                        output += Distributions.loglikelihood(dist,rels)
                    end
                end

                return output
            end
            local opts = NLopt.Opt(:LN_NELDERMEAD, length(params_vary))
            NLopt.max_objective!(opts, llhood)
            local max_f, max_x, ret = NLopt.optimize(opts, params_vary)
            out[jj] = max_f
        end

        println("Branch $ii calculated")
        plot(modifier,out .- maximum(out),ylims=(-3,1)) #
        df=1; llstar=-quantile(Chisq(df),0.95)/2 # Get the 95% sig threshold
        hline!([llstar],legend=false,lw=4,color=:gold)
        savefig("c_code/figures/PI/PI_numsamples=50_$ii.png")
    end
end
