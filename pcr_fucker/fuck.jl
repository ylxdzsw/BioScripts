@everywhere begin
    include("/haplox/users/zhangsw/deepsomatic/src/OhMyJulia.jl")
    include("/haplox/users/zhangsw/deepsomatic/src/Falcon.jl")
    include("/haplox/users/zhangsw/deepsomatic/src/Fire.jl")

    using OhMyJulia
    using Falcon
    using Fire
end

function pair!(reads::Vector{Read})
    namedict = Dict{String, Read}()
    for r in reads
        if r.qname in keys(namedict)
            r.mate = namedict[r.qname]
            namedict[r.qname].mate = r
            delete!(namedict, r.qname)
        else
            namedict[r.qname] = r
        end
    end
    filter(x->isdefined(x, :mate), reads)
end

function overlap(r1, r2)
    if r1.refID != r2.refID || abs(r1.pos - r2.pos) > 145
        return 0, b""
    end

    start, stop = max(r1.pos, r2.pos), min(r1.pos + calc_distance(r1), r2.pos + calc_distance(r2)) - 1

    if start >= stop
        return 0, b""
    end

    seq = try
         r1.seq[car(calc_read_pos(r1, start)):car(calc_read_pos(r1, stop))]
    catch e
        show(STDERR, e)
        show(STDERR, r1)
        show(STDERR, r2)
        return 0, b""
    end

    if isempty(seq)
        println(STDERR, "fuck empty")
        show(STDERR, r1)
        show(STDERR, r2)
        return 0, b""
    end

    if seq != try
        r2.seq[car(calc_read_pos(r2, start)):car(calc_read_pos(r2, stop))]
    catch e
        show(STDERR, e)
        show(STDERR, r1)
        show(STDERR, r2)
        return 0, b""
    end
        println(STDERR, "overlap not the same")
        show(STDERR, r1)
        show(STDERR, r2)
        return 0, b""
    end

    start, seq
end

# |-----32-----|--8--|----24----|
#      pos       chr     code
function get_clusters(bam)
    dict = Dict{Int64, Vector{String}}()
    reads = collect(filter(r->r.flag & 0x90c == 0 && 0 <= r.refID <= 255 && r.pos >= 0, BamLoader(bam))) |> pair!

    for r in reads @when pointer_from_objref(r) < pointer_from_objref(r.mate)
        pos, seq = overlap(r, r.mate)

        pos == 0 && continue

        code = UInt64(pos) << 8 | r.refID

        for base in SubString(r.qname, findlast(r.qname, ':')+1)
            code <<= 3
            code |= (UInt32(base) >> 1) & 0b0111
        end

        if code in keys(dict)
            push!(dict[code], seq)
        else
            dict[code] = [seq]
        end
    end
    collect(dict)
end

@everywhere function get_distance(cluster)
    d, g = 0, String[car(cluster)]

    for seq in cdr(cluster)
        if seq in g
            continue
        else
            d += minimum(ccall((:edit_distance, "edit_distance.so"), Cuint,
                               (Ptr{u8}, Cuint, Ptr{u8}, Cuint),
                               x, length(x), seq, length(seq)) for x in g)
            push!(g, seq)
        end
    end

    length(cluster), d, sum(map(length, cluster)), length(g)
end

@main function main(bam)
    list = @parallel vcat for (s, cluster) in get_clusters(bam)
        (s, get_distance(cluster)...)
    end

    foreach(prt, list)
end
