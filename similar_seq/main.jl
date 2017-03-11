using OhMyJulia
using BioDataStructures
using Fire
using HDF5

function bed_reader(bed, len=150)
    channel = Channel{Tuple{String, i32}}(256)

    @schedule begin
        for line in eachline(bed)
            chr, a, b = split(line, '\t')

            for i in parse(i32, a)-i32(74):parse(i32, b)+i32(75)
                put!(channel, (String(chr), i))
            end
        end
        close(channel)
    end

    channel
end

function find_N(seq::Bytes)
    s = IntRangeSet{Int}()

    for (i, b) in enumerate(seq)
        b == Byte('N') && push!(s, i)
    end

    s
end

a = rand(0x00:0x04, 65536*4 + 64)
b = rand(0x00:0x03, 128)
r = Array{Byte}(4)

ccall((:align, "ksw.so"), Void, (Ptr{Byte}, Ptr{Byte}, Ptr{Byte}, Cint), b, a, r, 4)
