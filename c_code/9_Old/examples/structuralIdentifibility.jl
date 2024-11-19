
# Import the library needed
using StructuralIdentifiability

# ODE model 1
# Measurements of the system are known for the planktonic,
# biofilm, ciliates and the ameoba.
ode1 = @ODEmodel(
    x1'(t) = -(1/eb)*( r2*x1(t)/(H2 + x1(t))*x2(t) + r3*x1(t)/(H3 + x1(t))*x3(t) ),
    x2'(t) = r2*x1(t)/(H2 + x1(t))*x2(t) + chi32*x3(t) - r4*x2(t)/(H4 + x2(t))*x4(t) - (a*chiMax + chiMin*x2(t))/(a+x2(t)),
    x3'(t) = r3*x1(t)/(H3 + x1(t))*x3(t) - chi32*x3(t) - r5*x3(t)/(H5 + x3(t))*x5(t) + (a*chiMax + chiMin*x2(t))/(a+x2(t)),
    x4'(t) = e4*r4*x2(t)/(H4 + x2(t))*x4(t),
    x5'(t) = e5*r5*x3(t)/(H5 + x3(t))*x5(t),

    y1(t) = x2(t),
    y2(t) = x3(t),
    y3(t) = x4(t),
    y4(t) = x5(t)
)

# ODE model 2
# Measurements of the system are known for the
# biofilm, ciliates and the ameoba.
ode2 = @ODEmodel(
    x1'(t) = -(1/eb)*( r2*x1(t)/(H2 + x1(t))*x2(t) + r3*x1(t)/(H3 + x1(t))*x3(t) ),
    x2'(t) = r2*x1(t)/(H2 + x1(t))*x2(t) + chi32*x3(t) - r4*x2(t)/(H4 + x2(t))*x4(t) - (a*chiMax + chiMin*x2(t))/(a+x2(t)),
    x3'(t) = r3*x1(t)/(H3 + x1(t))*x3(t) - chi32*x3(t) - r5*x3(t)/(H5 + x3(t))*x5(t) + (a*chiMax + chiMin*x2(t))/(a+x2(t)),
    x4'(t) = e4*r4*x2(t)/(H4 + x2(t))*x4(t),
    x5'(t) = e5*r5*x3(t)/(H5 + x3(t))*x5(t),

    y1(t) = x3(t),
    y2(t) = x4(t),
    y3(t) = x5(t)
)

# ODE model 2
# Measurements of the system are known for the
# biofilm only.
# ode3 = @ODEmodel(

#     x1'(t) = -(1/eb)*( r2*x1(t)/(H2 + x1(t))*x2(t) + r3*x1(t)/(H3 + x1(t))*x3(t) ),
#     x2'(t) = r2*x1(t)/(H2 + x1(t))*x2(t) + chi32*x3(t) - r4*x2(t)/(H4 + x2(t))*x4(t) - (a*chiMax + chiMin*x2(t))/(a+x2(t)),
#     x3'(t) = r3*x1(t)/(H3 + x1(t))*x3(t) - chi32*x3(t) - r5*x3(t)/(H5 + x3(t))*x5(t) + (a*chiMax + chiMin*x2(t))/(a+x2(t)),
#     x4'(t) = e4*r4*x2(t)/(H4 + x2(t))*x4(t),
#     x5'(t) = e5*r5*x3(t)/(H5 + x3(t))*x5(t),

#     y1(t) = x3(t)
# )

out1 = find_identifiable_functions(ode1)
println(out1)

out2 = find_identifiable_functions(ode2)
println(out2)

# out3 = find_identifiable_functions(ode3)
# println(out3)