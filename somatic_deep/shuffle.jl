#!/usr/bin/env julia

using OhMyJulia
using HDF5

function get_data_shuffled(files)
    channel = Channel()

    data = map(files) do file
        x = h5open(file)
        x, length(x)÷2
    end

    N = sum(cadr, data)

    get_data(i) = for (file, len) in data
        if i > len
            i -= len
        else
            return read(file, "x$i"), read(file, "y$i")
        end
    end

    @async while true
        sam = get_data.(rand(1:N, 16))
        x = mapreduce(car, vcat, sam)
        y = mapreduce(cadr, vcat, sam)
        put!(channel, (x, y))
    end

    channel
end
