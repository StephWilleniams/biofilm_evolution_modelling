
plt = plot()

for n = 1:10

    for ii = 1:5
        scatter!(plt,synthetic_data[ii,:,n],label=false)
    end

end

display(plt)
