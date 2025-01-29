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
function biofilm_predators(du, u, theta, t)
    Nbins = Int((length(u)-3)/2)
    chionVals = LinRange(0.000001,0.1,Nbins)
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H2, H3, H4, H5, chioff, chion = theta[1:end-1]
    # Calculate the ODEs to output
    du[1] = -r2*(u[1]/(u[1]+H2))*sum(u[2:1+N_bins])/e23 - r3*(u[1]/(u[1]+H3))*sum(u[2+N_bins:2+2*N_bins-1])/e23
    for i in 1:Nbins
        du[1+i] = r2*(u[1]/(u[1]+H2))*u[i+1] - r4*(u[i+1]/(u[i+1]+H4))*u[end-1] - chionVals[i]*u[i+1] + chioff*u[1+Nbins+i]
        du[1+Nbins+i] = r3*(u[1]/(u[1]+H3))*u[1+Nbins+i] - r5*(u[1+Nbins+i]/(u[1+Nbins+i]+H5))*u[end] + chionVals[i]*u[i+1] - chioff*u[1+Nbins+i]
    end
    du[end-1] = e4*r4*(sum(u[2:1+N_bins])/(sum(u[2:1+N_bins])+H4))*u[end-1]
    du[end] = e5*r5*(sum(u[2+N_bins:2+2*N_bins-1])/(sum(u[2+N_bins:2+2*N_bins-1])+H5))*u[end]
end

# Set the values required for the ODE solver
theta_true = (0.21, 0.007, 0.12, 0.09, # growth rates
              0.2, 0.5, 0.5 , # Efficiencies
              1.0, 1.0, 1.0, 0.1, # Half saturations
              0.005, 0.0005, # Attatchment parameters
              0.001) # Synthetic noise
C0 = 1.0; 
S0 = 0.0; 
#A0 = 0.0;
NpredRuns = 30
plt = plot()
colors = palette(cgrad([:green, :blue], NpredRuns)); 
global colind=1  # 3 colors
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
    tspan = (0.0, 24.0) # Time span of solution
    sample_frequency = 0.1 # Time-step at which to sample
    output_times = 0.0:sample_frequency:24.0 # Times at which to output the solutions
    N_times = length(output_times)
    bindRates = LinRange(0.000001,0.1,N_bins)
    RunNumber = 500
    EofRun = zeros(RunNumber)
    for runs in 1:RunNumber
        local prob = ODEProblem(biofilm_predators, u0, tspan, theta_true)
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
savefig("f_figures/predA.png") 