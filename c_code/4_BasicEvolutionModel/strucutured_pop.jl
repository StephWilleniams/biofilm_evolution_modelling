[println() for i in 1:10] # to make debugging easier
"""
Title: Modeling the effects of selection under predation
Author: Stephen Williams

Description: Here, we will similate a stochiometric model of a community of bacteria. 
Within the population there will be an initially normally distributed value of their biofilm formation rate.
The bacteria themselves will not be able to modify their biofilm formation rate value.
Within the system at t=0 there will be predators of varying quantities, which predate specific phenotypes.
At t = N*24 (hours) the system will be periodically reset.
We will attempt to quantify the mean of this biofilm formation rate.
"""

using DifferentialEquations, Plots

# Define the ODE of the system dynamics
# function BFP_model(du, u, params, t)
#     # Extract the different species
#     C, F, B, S, A = u
#     # Extract parameters
#     r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chioff = params
#     # Calculate the hill functions to output
#     CHfunct = C/(C+H23)
#     FHfunct = F/(F+H4)
#     BHfunct = B/(B+H5)
#     # Calculate the attatchment functions
#     chion = (chi_on_max*interaction_strength + chi_on_min*B)/(interaction_strength + B)
#     # Calculate the ODEs
#     du[1] = - r2*CHfunct*F - r3*CHfunct*B
#     du[2] = + e23*r2*CHfunct*F - r4*FHfunct*S - chion*F + chioff*B
#     du[3] = + e23*r3*CHfunct*B - r5*BHfunct*A + chion*F - chioff*B
#     du[4] = + e4*r4*FHfunct*S 
#     du[5] = + e5*r5*BHfunct*A
# end

function BFP_model(du, u, theta, t)
    
    Nbins = Int((length(u)-3)/2)
    chionVals = LinRange(0.000001,0.1,Nbins)

    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H2, H3, H4, H5, chi_off = theta[1:end-1]
    
    # Calculate the ODEs to output
    du[1] = -r2*(u[1]/(u[1]+H2))*sum(u[2:1+N_bins]) - r3*(u[1]/(u[1]+H3))*sum(u[2+N_bins:2+2*N_bins-1]) # 
    for i in 1:Nbins
        du[1+i] = e23*r2*(u[1]/(u[1]+H2))*u[i+1] - r4*(u[i+1]/(u[i+1]+H4))*u[end-1] - chionVals[i]*u[i+1] + chi_off*u[1+Nbins+i]
        du[1+Nbins+i] = e23*r3*(u[1]/(u[1]+H3))*u[1+Nbins+i] - r5*(u[1+Nbins+i]/(u[1+Nbins+i]+H5))*u[end] + chionVals[i]*u[i+1] - chi_off*u[1+Nbins+i]
    end
    du[end-1] = e4*r4*(sum(u[2:1+N_bins])/(sum(u[2:1+N_bins])+H4))*u[end-1]
    du[end] = e5*r5*(sum(u[2+N_bins:2+2*N_bins-1])/(sum(u[2+N_bins:2+2*N_bins-1])+H5))*u[end]
    
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
chi_off = 0.005*ND_T # Detachment rate

# Set the true values of the parameters
theta_true = (r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_off ) # Parameters
NpredRuns = 30

plt = plot()
colors = palette(cgrad([:green, :blue], NpredRuns)); 
global colind=1 # 3 colors

for A0 in LinRange(0.0,1.0,NpredRuns)

    N_bins = 1000
    global F0 = ones(N_bins)./N_bins 
    sig = N_bins/20
    for i in 1:N_bins
        F0[i] = exp(-0.5*((i-N_bins/2)/sig)^2)/(sig*sqrt(2*pi))
    end
    F0 = F0./sum(F0)
    B0 = zeros(N_bins)

    global u0 = vcat(C0,F0,B0,S0,A0) # IVP
    tspan = (0.0, 1.0) # Time span of solution
    sample_frequency = 0.01 # Time-step at which to sample
    output_times = 0.0:sample_frequency:1.0 # Times at which to output the solutions
    N_times = length(output_times)

    # Solve the ODE system
    bindRates = LinRange(0.000001,0.1,N_bins)

    RunNumber = 500 # Number of consecutive experiments to perform
    EofRun = zeros(RunNumber) # Expected values of the bind on rate for this day

    for runs in 1:RunNumber # Number of consecutive experiments to perform
        local prob = ODEProblem(BFP_model, u0, tspan, theta_true)
        local sol = solve(prob, Tsit5())
        global NofT = sol(output_times)
        global F0 = NofT[2+N_bins:2+2*N_bins-1,end]./sum(NofT[2+N_bins:2+2*N_bins-1,end])
        global u0 = vcat(C0,F0,B0,S0,A0) # IVP
        EofRun[runs] = sum(F0.*bindRates) # Expected values of the bind on rate for this day
    end

    # Plot the results
    plot!(plt,EofRun,lw=3,legend=false, color=colors[colind])
    global colind += 1

end

display(plt)
