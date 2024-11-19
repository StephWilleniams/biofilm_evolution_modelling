
using Plots

for i in 1:15
    scatter(A[i,:],Y)
    savefig("c_code/figures/fullSystem_sensitivity/reg_$i.png")
end
