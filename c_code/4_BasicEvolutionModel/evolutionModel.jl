using DifferentialEquations, Plots

# Define the ODE of the system dynamics
function BFP_model(du, u, params, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chioff = params
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



# Define a function of one input set only for simplicity - used for the sensitivity analysis
function model(theta)

    # Extract info from the input theta
    C0, F0, B0, S0, A0 = theta[1:5]
    params = theta[6:end]

    # Fixed values needed for the solver to work
    u0 = [C0, F0, B0, S0, A0] # IVP

    output_times = LinRange(0,1,100) # Solver times
    tspan = (output_times[1],output_times[end]) # End-points for the solver time

    # Run the solver for the system
    algo = AutoVern9(Rodas5()) # Optional stiff solver, required for optimisation
    prob = ODEProblem(BFP_model, u0, tspan, params) # Define the ODE problem
    sol = solve(prob, algo, saveat=output_times) # Solve the ODE problem

    return sol

    # # Set the parameter of interest for our outputs to check sensitivity with respect to
    # QOI = sol[2,end]

    # # Return the quantity of interest
    # return QOI

end

### Constants
# Non-dimensionalised scaling factors
ND_T = 24.0 # Non dimension time - hours
ND_C = 1.0 # Non dimension concentration - ug/ml

# Initial values
C0 = 1.0/ND_C # Initial concentration of Free bacteria
F0 = 1.0/ND_C # Initial concentration of Free predator
B0 = 0.0/ND_C # Initial concentration of Biofilm bacteria
S0 = 0.05/ND_C # Initial concentration of Ciliate predator
A0 = 0.01/ND_C # Initial concentration of Biofilm predator

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
H23 = 3.0*1.0/ND_C # Half saturations of bacteria
H4 = 3.0*1.0/ND_C # Half saturation of free predator
H5 = 8.0*0.1/ND_C # Half saturation of bound predator

# Attatchment parameters
chi_on_max = (8.0/3.0)*0.05*ND_T # Attatchment parameters
chi_on_min = (8.0/3.0)*0.0005*ND_T # Attatchment parameters
interaction_strength = (8.0/3.0)*0.01 # Interaction strength
chi_off = 0.005*ND_T # Detachment rate

# Set the true values of the parameters
theta_true = ( C0, F0, B0, S0,A0, r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off ) # Parameters

plt = plot()
NconditionRuns = 30
colors = palette(cgrad([:green, :blue], NconditionRuns)); 
global colind=1  # 3 colors
N_bins = 1000 # Number of strucutral bins in parameter which is allowed to do evo.
sig = N_bins/20 # Width of the initial condition distribution

for condition in LinRange(0.0,1.0,NpredRuns)

    global F0 = ones(N_bins) # initialise the free bacteria distribution
    for i in 1:N_bins # Loop over the bins
        F0[i] = exp(-0.5*((i-N_bins/2)/sig)^2)/(sig*sqrt(2*pi)) # Center the IVP of bacteria in the center of a normal dist.
    end
    F0 = F0./sum(F0) # Normalise the distribution
    B0 = zeros(N_bins) # Initial condition of biofilm bacteria

    # Use the condition of the loop here.
    global A0 = condition # Initial condition of biofilm predator

    global u0 = vcat(C0,F0,B0,S0,A0) # IVP

    # Data for the solver to get appropaite outputs
    tspan = (0.0, 1.0) # Time span of solution
    sample_frequency = 0.01 # Time-step at which to sample
    output_times = 0.0:sample_frequency:1.0 # Times at which to output the solutions
    N_times = length(output_times)

    # Solve the ODE system
    bindRates = LinRange(0.000001,0.1,N_bins) # Set the range of the bind on rates

    RunNumber = 500
    EofRun = zeros(RunNumber)

    for runs in 1:RunNumber
        
        local prob = ODEProblem(BFP_model, u0, tspan, theta_true)
        local sol = solve(prob, Tsit5())
        global NofT = sol(output_times)
        global F0 = NofT[2+N_bins:2+2*N_bins-1,end]./sum(NofT[2+N_bins:2+2*N_bins-1,end])
        global u0 = vcat(C0,F0,B0,S0,A0) # IVP

        EofRun[runs] = sum(F0.*bindRates)

    end


    plot!(plt,EofRun,lw=3,legend=false, color=colors[colind])
    global colind += 1

end

display(plt)