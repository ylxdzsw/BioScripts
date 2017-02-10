#!/usr/bin/env julia

using OhMyJulia
using HDF5

const data = mapreduce(readdlm, vcat, ARGS)
const shuf = shuffle(1:nrow(data))

h5open("data.hdf5", "w") do f
    write(f, "X", data[shuf, 2:end])
    write(f, "y", data[shuf, 1])
end
