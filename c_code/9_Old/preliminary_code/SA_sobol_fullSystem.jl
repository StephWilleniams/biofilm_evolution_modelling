
# Import the GS library
using GlobalSensitivity, QuasiMonteCarlo, DifferentialEquations

function biofilm_predators(du, u, theta, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H2, H3, H4, H5, chioff, a, chimax, chimin = theta
    # Other required values
    sig = 1 # Volume/surface area. 
    chion = (a*chimax + chimin*B)/(a+B) # Bind-on rate
    # Calculate the ODEs to output
    du[1] = -r2*(C/(C+H2))*F/e23 - sig*r3*(C/(C+H3))*B/e23
    du[2] = r2*(C/(C+H2))*F - r4*(F/(F+H4))*S - chion*F*sig + chioff*B*sig
    du[3] = r3*(C/(C+H3))*B - r5*(B/(B+H5))*A + chion*F - chioff*B
    du[4] = e4*r4*(F/(F+H4))*S 
    du[5] = e5*r5*(B/(B+H5))*A
end

# Define a function of one input set only for simplicity
function model(theta)

    # Fixed values needed for the solver to work
    u0 = [1.0,1.0,0.0,0.05,0.01] # IVP
    output_times = LinRange(0,24,25) # Solver times
    tspan = (output_times[1],output_times[end]) # End-points for the solver time

    # Run the solver for the system
    prob = ODEProblem(biofilm_predators, u0, tspan, theta) # Define the ODE problem
    sol = solve(prob, AutoVern9(Rodas5()), saveat=output_times) # Solve the ODE problem

    # Set the parameter of interest for our outputs to check sensitivity with respect to
    QOI = sol[2,end]

    # Return the quantity of interest
    return QOI

end

theta = [0.21,0.007,0.12,0.09, # growth rates
         0.2,0.5,0.5, # Efficiencies
         1.0,1.0,1.0,0.1, # Half saturations
         0.005,0.01,0.05,0.005] # Attatchment parameters

test = false
if test
    checker = model(theta)
    println(checker)
end

delta = 0.2
samples = 80000
lb = [(1-delta)*theta[i] for i in eachindex(theta)]
ub = [(1+delta)*theta[i] for i in eachindex(theta)]
sampler = GlobalSensitivity.SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)

res1 = gsa(model,Sobol(order=[0,1,2]),A,B)

plt = bar(1:15,res1.ST,legend=false, label="ST")
bar!(plt,1:15,res1.S1,legend=false, label="S1")
#bar!(plt,1:2:30,res1.S2,legend=false, label="S2")

display(plt)
