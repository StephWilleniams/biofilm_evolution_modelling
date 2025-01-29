
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Super-Case 1 - Complex bind on dynamics

## State variables measured: C, F, B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y1=x1,
y2=x2,
y3=x3,
y4=x4,
y5=x5

RESULTS -- Fully Gloabally Identifiable.

## State variables measured: F, B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y2=x2,
y3=x3,
y4=x4,
y5=x5

RESULTS -- H23, e23, r2, r3, C(0) not identifiable.

## State variables measured: B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y3=x3,
y4=x4,
y5=x5

RESULTS -- H23, e23, r2, r3, C(0) not identifiable. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Super-Case 2 - Simple bind on dynamics

## State variables measured: C, F, B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - chi_on*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + chi_on*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y1=x1,
y2=x2,
y3=x3,
y4=x4,
y5=x5

RESULTS -- Fully Gloabally Identifiable.

## State variables measured: F, B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - chi_on*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + chi_on*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y2=x2,
y3=x3,
y4=x4,
y5=x5

RESULTS -- H23, e23, r2, r3, C(0) not identifiable. 

## State variables measured: B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - chi_on*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + chi_on*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y3=x3,
y4=x4,
y5=x5

RESULTS -- H23, e23, r2, r3, C(0) not identifiable.  
