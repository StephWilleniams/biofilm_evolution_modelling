# Title: Code for running the genetic algorithm for model selection
""" 
Description: This code uses a genetic algorithm to optimise the sampling time points 
for a model of a biological system. 
The model is a system of ordinary differential equations that describes the dynamics of free bacteria, 
bound bacteria, and predators. 
The genetic algorithm optimises the sampling schedule to maximise the information content of the data 
for a parameter of interest. 
The result minimises the width of the 95% confidence interval for the parameter estimates.
"""

for i in 1:3; println(" "); end # Blank lines for separation in terminal

# Define constants and modules
include("preamble.jl")

# Set initial organism states
for i in 1:N_organisms
    global organisms[i] = GeneticAlgorithm.organism(GeneticAlgorithm.GenerateInitialSchedule(t_start,t_end,dt,n_times),i,0.0)
end

# N_generations generation of model selection
for gen = 1:N_generations
    for i in 1:N_organisms
        println("Generation: ", gen, " Organism: ", i)
        synthetic_data = GeneticAlgorithm.SyntheticPopulation(GeneticAlgorithm.SolvePopulation(GeneticAlgorithm.BFP_model, u0, organisms[i], theta_baseline, ode_algo), dist, N_synthetics)
        global organisms[i].fitness = GeneticAlgorithm.PLB(index_of_interest, PL_range, n_PLs, GeneticAlgorithm.BFP_model, ode_algo, u0, dist, organisms[i], theta_baseline, initial_varying_params, varying_indices, synthetic_data)
    end
    global organisms = GeneticAlgorithm.OrderOrganisms(organisms)
    global organisms = GeneticAlgorithm.CloneOrganisms(organisms,cloning_noise,dt,t_start,t_end)
end

GeneticAlgorithm.PrintOrganisms(organisms)

plt = plot()
for i = eachindex(organisms)
    scatter!(plt,organisms[i].schedule,1:20,label=false,xlabel="Time",ylabel="Concentration")
end
