# Load the required packages
using DifferentialEquations, Plots, LatinHypercubeSampling, Statistics, LaTeXStrings

# Define the ODE system for the BFP model
function BFP_model(du, u, params, t)

    # Extract the different species
    C, F, B, S, A = u

    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chioff = params

    # Re-fix parameters not of interest to the sensitivity
    # r2 = 0.21*ND_T
    # r3 = 0.007*ND_T
    # r4 = 0.12*ND_T
    # r5 = 0.09*ND_T
    # e23 = 0.2

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
    u0 = theta[1:5]
    #u0 = [ ND_C, 1.0*ND_C, 0.0, 3.0*0.01*ND_C, 8.0*0.001*ND_C ] # Fixed values
    #u0 = [ ND_C, 0.3*ND_C, 0.0, 10.0*0.1*ND_C, 8.0*0.1*ND_C ] # IVP - 2C
    params = theta[6:end]

    # Fixed values needed for the solver to work
    output_times = LinRange(0,1,100) # Solver times
    tspan = (output_times[1],output_times[end]) # End-points for the solver time

    # Run the solver for the system
    algo = AutoVern9(Rodas5()) # Optional stiff solver, required for optimisation
    prob = ODEProblem(BFP_model, u0, tspan, params) # Define the ODE problem
    sol = solve(prob, algo, saveat=output_times) # Solve the ODE problem

    # Set the parameter of interest for our outputs to check sensitivity with respect to
    QOI = sol[2,end]

    # Return the quantity of interest
    return QOI

end

# Set parameters

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
H4  = 3.0*1.0/ND_C # Half saturation of free predator
H5  = 8.0*0.1/ND_C # Half saturation of bound predator

# Attatchment parameters
chi_on_max = (8.0/3.0)*0.05*ND_T # Attatchment parameters
chi_on_min = (8.0/3.0)*0.0005*ND_T # Attatchment parameters
interaction_strength = (8.0/3.0)*0.01 # Interaction strength
chi_off = 0.005*ND_T # Detachment rate

# Set the true values of the parameters
theta_true = ( C0, F0, B0, S0, A0, r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off ) # Parameters
theta_names = [L"C_0", L"F_0", L"B_0", L"S_0", L"A_0", L"r_F", L"r_B", L"r_S", L"r_A", L"e_{FB}", L"e_S", L"e_A", L"H_{FB}", L"H_S", L"H_A", L"\chi_{on}^{max}", L"\chi_{on}^{min}", L"a", L"\chi_{off}"]

# Get samples over the space of interest
# Define the number of samples and dimensions
n_samples = 100000  # Number of sample points
n_dims = 19 # Number of variables
delta = 0.5; # Proportional size of the space to sample around theta_zero
samples = 1.0 .+ ((randomLHC(n_samples, n_dims)./n_samples) .- 0.5) .* 2.0 .* delta

# Calculate the variance of the quantity of interest
output_values = zeros(n_samples)
for i in 1:n_samples
    output_values[i] = model(theta_true .* samples[i,:])
end

histogram_plot = true
if histogram_plot
    # Plot the histogram of the quantity of interest
    plt = histogram(output_values, bins=100, legend=false)
    xlabel!("Quantity of interest")
    ylabel!("Frequency")
    title!("Variance of the quantity of interest")
    display(plt)
    savefig("f_figures/4B_variance_histogram.png")
end

univariate_plot = false
if univariate_plot
    for quanitity_index in 1:19
        plt = scatter(samples[:,quanitity_index],output_values, alpha=0.1, legend=false)
        xlabel!(L"\Delta" * theta_names[quanitity_index])
        ylabel!("Quantity of interest")
        title!("Variance of the quantity of interest")
        display(plt)
    end
end

# Print the variance of the quantity of interest
println("Variance of the quantity of interest: ", var(output_values))
