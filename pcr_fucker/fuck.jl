@everywhere begin
    include("/home/zhangsw/deepsomatic/src/OhMyJulia.jl")
    include("/home/zhangsw/deepsomatic/src/Falcon.jl")
    include("/home/zhangsw/deepsomatic/src/Fire.jl")

    using OhMyJulia
    using Falcon
    using Fire
end

# |-----32-----|--8--|----24----|
#      pos       chr     code
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

@everywhere function get_distance(cluster)
    n, d, g = 1, 0, String[car(cluster)]

    for seq in cdr(cluster)
        n += 1

        if seq in g
            continue
        else
            d += minimum(ccall((:edit_distance, "edit_distance.so"), Cuint,
                               (Ptr{u8}, Cuint, Ptr{u8}, Cuint),
                               x, length(x), seq, length(seq)) for x in g)
            push!(g, seq)
        end
    end

    n, d, length(g)
end

@main function main(bam)
    s, clusters = get_clusters(bam)
    n, d, l = @parallel (x,y)->map(+, x, y) for cluster in clusters
        get_distance(cluster)
    end

    info("Average Error per Copy: $d ÷ $n = $(100d/n)%")
    info("Average Number of Unique Seq for Each DNA: $l ÷ $(length(clusters)) = $(l/length(clusters))")
    info("Average Number of Unique DNA of Each Position: $s")
end
