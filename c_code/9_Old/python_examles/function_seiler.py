# Function to give the various dynamics in the system

# Calculate the RHS of the nutrient dynamics
def dCdt(t,Y,params,constants):

    # Extract parameters
    r2, r3, \
    eb, \
    H2, H3, \
    chi, a, chimax, chimin, \
    r4, r5, \
    eff4, eff5, \
    H4, H5, = params

    sigma = 8/3

    M12 = (1/eb)*r2*Y[1]*Y[0]/(H2+Y[0]) # Nutrient loss to growth of free
    M13 = (1/eb)*sigma*r3*Y[2]*Y[0]/(H3+Y[0]) # Nutrient loss to growth of biofilm

    return -(M12 + M13)

# Calculate the RHS of the free bacteria dynamics
def dFdt(t,Y,params,constants):

    # Extract parameters
    r2, r3, \
    eb, \
    H2, H3, \
    chi, a, chimax, chimin, \
    r4, r5, \
    eff4, eff5, \
    H4, H5, = params

    sigma = 8/3

    M12 = r2*Y[1]*Y[0]/(H2+Y[0]) # Nutrient loss to growth of free
    G2 = M12 # Growth of free by bacteria consumption
    M24 = r4*Y[3]*Y[1]/(H4+Y[1]) # Consumption of free bacteria by predators
    T23 = sigma*chi*Y[1] # Attachment rate 
    T32 = sigma*(a*chimax+chimin*Y[2])/(a+Y[2]) # Dettachment rate 

    return (G2 + T32) - (M24 + T23)

# Calculate the RHS of the biofilm bacteria dynamics
def dBdt(t,Y,params,constants):

    # Extract parameters
    r2, r3, \
    eb, \
    H2, H3, \
    chi, a, chimax, chimin, \
    r4, r5, \
    eff4, eff5, \
    H4, H5, = params

    M13 = r3*Y[2]*Y[0]/(H3+Y[0]) # Nutrient loss to growth of biofilm
    G3 = M13 # Growth of biofilm by bacteria consumption
    M35 = r5*Y[4]*Y[2]/(H5+Y[2]) # Biofilm loss to growth of predator
    T23 = chi*Y[1] # Attachment rate 
    T32 = (a*chimax+chimin*Y[2])/(a+Y[2]) # Dettachment rate 

    return (G3 + T23) - (M35 + T32)

# Calculate the RHS of the free predator dynamics
def dSdt(t,Y,params,constants):

    # Extract parameters
    r2, r3, \
    eb, \
    H2, H3, \
    chi, a, chimax, chimin, \
    r4, r5, \
    eff4, eff5, \
    H4, H5, = params

    M24 = r4*Y[3]*Y[1]/(H4+Y[1]) # Consumption of free bacteria by predators
    G4 = eff4*M24 # Growth of free predator by free consumption

    return G4

# Calculate the RHS of the biofilm predator dynamics
def dTdt(t,Y,params,constants):

    # Extract parameters
    r2, r3, \
    eb, \
    H2, H3, \
    chi, a, chimax, chimin, \
    r4, r5, \
    eff4, eff5, \
    H4, H5, = params

    M35 = r5*Y[4]*Y[2]/(H5+Y[2]) # Biofilm loss to growth of predator
    G5 = eff5*M35 # Growth of biofilm predator by biofilm consumption

    return G5

# Calculate the RHS of the full system
def dXdt(t,Y,params,constants):
    return [dCdt(t,Y,params,constants),dFdt(t,Y,params,constants),dBdt(t,Y,params,constants),dSdt(t,Y,params,constants),dTdt(t,Y,params,constants)]