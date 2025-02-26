for i in 1:10; println(" "); end

# Define constants and modules
include("preamble.jl")

# Set initial organism states
for i in 1:N_organisms
    global organisms[i] = GeneticAlgorithm.organism(GeneticAlgorithm.GenerateInitialSchedule(t_start,t_end,dt,n_times),i,0.0)
end

# N_generations generation of model selection
for gen = 1:N_generations
    for i in 1:N_organisms
        synthetic_data = GeneticAlgorithm.SyntheticPopulation(GeneticAlgorithm.SolvePopulation(GeneticAlgorithm.BFP_model, u0, organisms[i], theta_baseline, ode_algo), dist, N_synthetics)
        organisms[i].fitness = GeneticAlgorithm.PLB(index_of_interest, PL_range, n_PLs, GeneticAlgorithm.BFP_model, ode_algo, u0, dist, organisms[i], theta_baseline, initial_varying_params, varying_indices, synthetic_data)
    end
    global organisms = GeneticAlgorithm.OrderOrganisms(organisms)
    global organisms = GeneticAlgorithm.CloneOrganisms(organisms,cloning_noise,dt,t_start,t_end)
end

GeneticAlgorithm.PrintOrganisms(organisms)
