
#############################################################################################################

# # Inclusions
using DifferentialEquations, Plots, Distributions, Random, NLopt

#############################################################################################################

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

μ_for_mean(mu, sigma) = log( (mu^2)/sqrt(mu^2 + sigma^2) )
σ_for_sigma(mu, sigma) = log( 1 + (sigma/mu)^2 )

#############################################################################################################

# Set the seed
#Random.seed!(1234) # Use for repeatability

#############################################################################################################

# PI factors
N_samples = 1000 # Number of synthetic samples
sample_frequency = 0.01 # Time-step at which to sample
μ = 1.0; σ = 0.1 # Mean and standard deviation of the noise
modification_scaling = 0.1 # Percent change range for uniform scaling of true parameters
#ode_algo = AutoVern9(Rodas5())
ode_algo = Tsit5()
tspan = (0.0, 1.0) # Time span of solution
sample_frequency = 0.01 # Sample frequency
output_times = tspan[1]:sample_frequency:tspan[2] # Times at which to output the solutions
N_times = length(output_times) # Number of output times

# Non-dimensionalised scaling factors
ND_T = 24.0 # Non dimension time - hours
ND_C = 1.0 # Non dimension concentration - ug/ml
e23 = 0.2 # Free bacteria growth efficiency
e4 = 0.5 # Bound bacteria growth efficiency
e5 = 0.5 # Predator growth efficiency
r2 = (1/e23)*0.21*ND_T # Free bacteria growth rate - cells/hr * ND_hrs
r3 = (1/e23)*0.007*ND_T # Bound bacteria growth rate - cells/hr * ND_hrs
r4 = (1/e4)*0.12*ND_T # Predator growth rate - cells/hr * ND_hrs
r5 = (1/e5)*0.09*ND_T # Predator death rate - 1/hrs * ND_hrs
H23 = 3.0/ND_C # Half saturations of bacteria
H4 = 3.0/ND_C # Half saturation of free predator
H5 = 0.8/ND_C # Half saturation of bound predator
chi_on_max = (8.0/3.0)*0.05*ND_T # Attatchment parameters
chi_on_min = (8.0/3.0)*0.0005*ND_T # Attatchment parameters
interaction_strength = (8.0/3.0)*0.01 # Interaction strength
chi_off = 0.005*ND_T

#############################################################################################################

# Set the true values of the parameters
theta_true = (r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off) # Parameters
# scalings = (1 .+ modification_scaling*2*(rand(length(theta_true)).-0.5))
# theta_guess = theta_true .* scalings # Initial guess for parameters
u0 = [ ND_C, 1.0*ND_C, 0.01, 3.0*0.01*ND_C, 8.0*0.001*ND_C ] # IVP - 2B

#############################################################################################################

dist = LogNormal(μ_for_mean(μ, σ), σ_for_sigma(μ, σ))
opts = NLopt.Opt(:LN_NELDERMEAD, length(1))
NLopt.xtol_rel!(opts,1e-12)

redo_optimisation = true
if redo_optimisation
    N_repeats = 10
    outcomes = zeros(N_repeats,14,2)
    for repeats = 1:N_repeats
        for index_optimised = 1:14
            local scalings = (1 .+ modification_scaling*2*(rand(length(theta_true)).-0.5))
            local theta_guess = theta_true .* scalings # Initial guess for parameters
            # Define functions for the maximiser
            function llhood(params,_)   
                Param_to_optimise = params
                if index_optimised == 1
                    theta_start = [Param_to_optimise...,theta_true[2:end]...]
                elseif index_optimised == 14
                    theta_start = [theta_true[1:13]...,Param_to_optimise...]
                else
                    theta_start = [theta_true[1:index_optimised-1]...,Param_to_optimise...,theta_true[index_optimised+1:end]...]
                end
                prob = ODEProblem(BFP_model, u0, tspan, theta_start)
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
            NLopt.max_objective!(opts, llhood)
            local f_opt, x_opt, ret = NLopt.optimize(opts, [theta_guess[index_optimised]])
            x_opt = x_opt[1]
            outcomes[repeats,index_optimised,1] = x_opt/theta_true[index_optimised]
            outcomes[repeats,index_optimised,2] = x_opt/theta_guess[index_optimised]
        end
    end
end

plt = plot()
for index_optimised = 1:14
    local plt = plot()
    for ii = 1:N_repeats
        x = [ii,ii]
        y = [outcomes[ii,index_optimised,2],outcomes[ii,index_optimised,1]]
        GR.setarrowsize(1)
        plot!(plt,x,y, marker =:circle, arrow=(:closed, 0.3), label = false)
        savefig("c_code/figures/PI/OPAT_$index_optimised.png")
    end
end

