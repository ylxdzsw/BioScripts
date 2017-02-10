using OhMyJulia
using HDF5

const X, y = h5open("data.hdf5", "r") do f
    read(f, "X"), read(f, "y")
end

const d = X[:, end-10] .+ X[:, end-11]
const f = X[:, end-9]

const maxd = maximum(d)

residual = 0

for (stop, start) in Δ(tuple, [-1, .1maxd, .2maxd, .3maxd, .4maxd, .5maxd, .6maxd, maxd])
    ind = start .< d .<= stop
    frequency = f[ind]
    issomatic = y[ind]
    ind = sortperm(frequency)
    frequency = frequency[ind]
    issomatic = issomatic[ind]
    cutoff = findmax(cumsum(issomatic .- .5)) |> cadr
    residual += sum(1 - issomatic[1:cutoff]) + sum(issomatic[cutoff+1:end])
    prt(round(Int, start), round(Int, stop), frequency[cutoff])
end

println("Overall correct rate: ", 1 - residual / length(y))
