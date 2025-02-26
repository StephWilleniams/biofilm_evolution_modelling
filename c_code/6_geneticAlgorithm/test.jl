
# Define constants and modules
include("preamble.jl")

organo = GeneticAlgorithm.organism(GeneticAlgorithm.GenerateInitialSchedule(t_start,t_end,dt,n_times),1,0.0)
NofT = GeneticAlgorithm.SolvePopulation(GeneticAlgorithm.BFP_model, u0, organo, theta_baseline, ode_algo)
Synthetic_data = GeneticAlgorithm.SyntheticPopulation(NofT, dist, N_synthetics)

index_of_interest = 11
parameter_range_of_interest = 0.3
initial_params = [0.05,0.05,0.05]
varying_indices = [12,13,14]
n_PLs = 31

output = GeneticAlgorithm.UnivariateProfileLikelihood(index_of_interest, parameter_range_of_interest, n_PLs, GeneticAlgorithm.BFP_model, 
ode_algo, u0, dist, organo, theta_baseline, initial_params, varying_indices, 
Synthetic_data)

println(GeneticAlgorithm.ZeroIntersectionsSpline(output, parameter_range_of_interest, n_PLs))
