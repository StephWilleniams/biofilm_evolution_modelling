module GeneticAlgorithm  # Name of module

using DifferentialEquations, Plots, Distributions, 
Random, NLopt, CSV, DataFrames, Interpolations, Roots

# Random.seed!(1) 

# Log-normal distribution mean
function μ_for_mu(mu, sigma) 
    return log( (mu^2)/sqrt(mu^2 + sigma^2) )
end

# Log-normal distribution noise
function σ_for_sigma(mu, sigma)
    return log( 1 + (sigma/mu)^2 )
end

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

# Define organism structure
mutable struct organism
    schedule::Vector{Float64}
    ranking::Int
    fitness::Float64
end

# Initial schedule generation
function GenerateInitialSchedule(t_start, t_end, dt, n_times)
    schedule = zeros(n_times)
    schedule[1] = t_start # Set first time to t_start
    test_value = 0.0 # Initialize test value store
    for i in 2:n_times
        check = true
        reseed_attempts = 0 # Initialize reseed attempts
        while check
            test_value = t_start + (t_end - t_start)*rand() # Generate random time
            for j in 1:i-1
                if abs(test_value - schedule[j]) < dt
                    check = true
                    reseed_attempts += 1
                    if reseed_attempts > 100000000 # Attemps 100M reseeds before quitting.
                        error("Cannot create schedule, please lower number of sample points.")
                    end
                    break
                else
                    check = false
                end
            end
        end
        schedule[i] = test_value
    end
    schedule = sort(schedule)
    return schedule
end

# Manually set fitness
function SetFitness(organ::organism,fitness)
    organ.fitness = fitness
    return organ
end

# Solve the population ODE system for a given organism schedule and parameters
function SolvePopulation(model, u0, organ::organism, theta, solver_algorithm)
    tspan = [minimum(organ.schedule), maximum(organ.schedule)]
    prob = ODEProblem(model, u0, tspan, theta)
    solved_system = solve(prob, solver_algorithm)
    sol = solved_system(organ.schedule)
    return sol
end

# Generate synthetic data from the population solution
function SyntheticPopulation(sol, dist, N_synthetics)
    synthetic_data = zeros(size(sol,1), size(sol,2), N_synthetics)
    N_species = size(sol,1); N_times = size(sol,2)
    for i = 1:N_species
        for n = 1:N_synthetics
            synthetic_data[i,:,n] = sol[i,:] .* rand(dist, N_times) # Add noise to the data
        end
    end
    return synthetic_data
end

# Log-likelihood function for optimization
function LogLikelihood(params, ode_model, solver_algorithm, u0, dist, organ::organism, synthetic_data)
    llhood = 0.0 # Initialize the log-likelihood value
    localPoplulationSolution = SolvePopulation(ode_model, u0, organ, params, solver_algorithm)
    N_species, _ = size(localPoplulationSolution)
    N_repeats = size(synthetic_data)[3]
    for i = 1:N_species
        for j = 1:N_repeats
            rels = synthetic_data[i,:,j]./localPoplulationSolution[i,:]
            llhood += Distributions.loglikelihood(dist,rels)
        end
    end
    return llhood
end

# Overwrite fixed parameters with varying parameter guesses at specified indices
function ConstructParameters(varying_params, fixed_params, varied_param_indices)
    N_varied_params = length(varying_params)
    params = copy(fixed_params)
    for i in 1:N_varied_params
        params[varied_param_indices[i]] = varying_params[i]
    end
    return params
end

# Profile likelihood optimization over log-likelihood function for some varying parameters
function ProfileLikelihoodOptimisation(ode_model, solver_algorithm, u0, dist, organ::organism, fixed_paramaters, initial_varying_params, varying_indices, synthetic_data)
    function Objective_wrapper(varying_params::Vector{Float64}, _ )
        params = ConstructParameters(varying_params, fixed_paramaters, varying_indices)
        return LogLikelihood(params, ode_model, solver_algorithm, u0, dist, organ, synthetic_data)
    end
    optimisation_object = Opt(:LN_NELDERMEAD, length(initial_varying_params)) # Choose your algorithm and dimension
    NLopt.max_objective!(optimisation_object, Objective_wrapper)
    (best_llhood, best_params, ret) = optimize(optimisation_object, initial_varying_params)
    return best_llhood, best_params, ret
end

function UnivariateProfileLikelihood(index_of_interest, parameter_range_of_interest, n_PLs, ode_model, solver_algorithm, u0, dist, organ::organism, fixed_paramaters, initial_varying_params, varying_indices, synthetic_data)
    llhood_values = zeros(n_PLs)
    lower_bound = (1-parameter_range_of_interest) * fixed_paramaters[index_of_interest]
    upper_bound = (1+parameter_range_of_interest) * fixed_paramaters[index_of_interest]
    univariate_range = range(lower_bound, upper_bound, length = n_PLs)
    for i = 1:n_PLs
        local_fixed_paramaters = copy(fixed_paramaters)
        local_fixed_paramaters[index_of_interest] = univariate_range[i]
        (llhood_val, _ , _ ) = ProfileLikelihoodOptimisation(ode_model, solver_algorithm, u0, dist, organ::organism, local_fixed_paramaters, initial_varying_params, varying_indices, synthetic_data)
        llhood_values[i] = llhood_val
    end
    return llhood_values .- maximum(llhood_values) .+ quantile(Chisq(1),0.95)/2 # Shift by the 95% sig threshold (value = 0 now corresponds to 95% confidence interval)
end

function ZeroIntersectionsSpline(arr, parameter_range_of_interest, n_PLs)
    lower_bound = (1-parameter_range_of_interest)
    upper_bound = (1+parameter_range_of_interest)
    x = range(lower_bound, upper_bound, length = n_PLs)
    y = arr
    fit = CubicSplineInterpolation(x, y)
    function spline_zero(x)
        return fit(x)
    end

    if fit(lower_bound) * fit(1.0) > 0 || fit(1.0) * fit(upper_bound) > 0
        return NaN, NaN
    end

    root1 = find_zero(spline_zero,(lower_bound, 1.0))
    root2 = find_zero(spline_zero,(1.0, upper_bound))

    return root1, root2
end

function PLB(index_of_interest, parameter_range_of_interest, n_PLs, ode_model, solver_algorithm, u0, dist, organ::organism, fixed_paramaters, initial_varying_params, varying_indices, synthetic_data)
    output = UnivariateProfileLikelihood(index_of_interest, parameter_range_of_interest, n_PLs, ode_model, solver_algorithm, u0, dist, organ, fixed_paramaters, initial_varying_params, varying_indices, synthetic_data)
    root1,root2 = ZeroIntersectionsSpline(output, parameter_range_of_interest, n_PLs)
    return abs(root1 - root2)
end

# Order organisms by fitness metric
function OrderOrganisms(organisms::Vector{organism}, is_reversed::Bool = false)
    organisms = sort(organisms, by = x -> x.fitness, rev = is_reversed)
    for i in eachindex(organisms)
        organisms[i].ranking = i
    end
    return organisms
end

# Clone organisms according to fixed rule
function CloneOrganisms(organisms::Vector{organism}, sig, dt, t_start, t_end)
    highest_index_kept = Int(min(ceil(100/14)^2,length(organisms)))
    cloned_organisms = copy(organisms[1:highest_index_kept])
    for i in 1:length(organisms)
        new_index = Int(min(ceil(i/14)^2,highest_index_kept))
        if organisms[i].ranking == 1
            organisms[i].schedule = cloned_organisms[new_index].schedule
        else
            organisms[i].schedule[2:end] = cloned_organisms[new_index].schedule[2:end] + sig*randn(length(organisms[i].schedule)-1)
        end

        for j in 2:length(organisms[i].schedule)
            organisms[i].schedule[j] = min(max(organisms[i].schedule[j],t_start),t_end)
            if organisms[i].schedule[j] - organisms[i].schedule[j-1] < dt
                organisms[i].schedule[j] = organisms[i].schedule[j-1] + dt
            end
        end
        organisms[i].schedule[end] = min(organisms[i].schedule[end],t_end) 
    end
    return organisms
end

# Print organism information for debugging
function PrintOrganisms(organisms::Vector{organism})
    println("Organisms as follows:")
    for organ in organisms
        println("Schedule: ", organ.schedule)
        println("Ranking: ", organ.ranking)
        println("Fitness: ", organ.fitness)
        println(" ")
    end
end

# Save the worst organism
function SaveWorstOrganism(organisms::Vector{organism}) 
    filename = "worst_organism.csv"
    N_organisms = length(organisms)
    worst_index = findfirst(org -> org.ranking == N_organisms, organisms)
    df = DataFrame(organisms[worst_index], :auto)
    CSV.write(filename, df)
end

# Save the best organism
function SaveBestOrganism(organisms::Vector{organism})
    filename = "worst_organism.csv"
    worst_index = findfirst(org -> org.ranking == 1, organisms)
    df = DataFrame(organisms[worst_index], :auto)
    CSV.write(filename, df)
end

# Export functions
export organism, 
GenerateInitialSchedule, SetFitness, 
SolvePopulation, SyntheticPopulation, 
LogLikelihood, ConstructParameters, ProfileLikelihoodOptimisation, UnivariateProfileLikelihood, 
ZeroIntersectionsSpline, PLB, OrderOrganisms, CloneOrganisms, PrintOrganisms, 
SaveWorstOrganism, SaveBestOrganism

end # End of module GeneticAlgorithm