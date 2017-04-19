include("/home/zhangsw/deepsomatic/src/OhMyJulia.jl")
include("/home/zhangsw/deepsomatic/src/Falcon.jl")
include("/home/zhangsw/deepsomatic/src/Fire.jl")

using OhMyJulia
using Falcon
using Fire

function get_clusters(bam)
    dict = Dict{Int64, Vector{String}}()
    for r in BamLoader(bam) @when r.flag & 0x900 == 0 && 0 <= r.refID <= 255 && r.pos >= 0
        code = UInt64(r.pos) << 8 | r.refID

        for base in SubString(r.qname, findlast(r.qname, ':')+1)
            code <<= 3
            code |= (UInt32(base) >> 1) & 0b0111
        end

        if code in keys(dict)
            push!(dict[code], r.seq)
        else
            dict[code] = [r.seq]
        end
    end
    mean(values(groupby(x->x>>>32, (x,y)->x+1, 0, keys(dict)))), collect(values(dict))
end
