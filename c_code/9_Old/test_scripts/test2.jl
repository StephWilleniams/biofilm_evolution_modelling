
using DifferentialEquations, Plots

function exponential_growth!(du, u, p, t)
    du[1] = p[1] * u[1]
end

# Parameters
k = 0.2  # Growth rate constant
y0 = [100]  # Initial population

# Define the ODE problem
tspan = (0.0, 10.0)  # Time interval
p = [k]  # Parameters for the ODE
prob = ODEProblem(exponential_growth!, y0, tspan, p)

# Solve the ODE
sol = solve(prob, RK4())

# Plot the results
plot(sol, label="Population", xlabel="Time", ylabel="Population")