#!/usr/bin/env julia

using OhMyJulia
using Insane
using IntRangeSets
using Base.Collections
import Base: start, next, done, iteratorsize, eltype,
             getindex, setindex!, show, ==, hash

include(rel"../falcon/mut.jl")
include(rel"../falcon/read.jl")
include(rel"../falcon/bam.jl")
include(rel"../falcon/sam.jl")
include(rel"../falcon/pair.jl")

const bam = Bam(STDIN)
const reads = collect(bam)
const thresholds = try parse(Int, ARGS[1]) catch 10 end

type State
    window::Vector{Int32} # end pos of current reads
    region::IntRangeSet{i32}
    dna::Set{Tuple{Int32, Int32}} # start, length; negative pos means single, negative length means reverse
    chr::Int32
    pos::Int32
    State() = new(Int32[], IntRangeSet{i32}(), Set{Tuple{Int32, Int32}}(), -1, -1)
end

function add_read!(p::State, r::Read)
    while length(p.window) >= thresholds && p.pos < r.pos
        produce_region!(p)
    end

    p.pos = r.pos

    while !isempty(p.window) && p.window[1] < p.pos
        heappop!(p.window)
    end

    push!(p.window, r.pos + calc_distance(r) - 1)
end

function produce_region!(p::State)
    push!(p.region, p.pos:p.window[1])
    p.pos = p.window[1] + 1
    while !isempty(p.window) && p.window[1] < p.pos
        heappop!(p.window)
    end
end

function flush!(p::State)
    while length(p.window) >= thresholds
        produce_region!(p)
    end

    empty!(p.window)
    empty!(p.dna)

    if p.chr >= 0
        refname = car(bam.refs[p.chr+1])
        foreach(p.region) do range
            prt(refname, range.start, range.stop)
        end
    end
end

function main(p::State)
    for r in reads
        contig = if isdefined(r, :mate)
            r.pos, r.tlen
        else
            factor = (r.flag & 0x0010 != 0) ? -1 : 1
            -r.pos, i32(factor * calc_distance(r))
        end

        contig in p.dna && continue

        if r.refID != p.chr
            flush!(p)
            p.chr = r.refID
            p.pos = r.pos
        end
        push!(p.dna, contig)
        add_read!(p, r)
    end
    flush!(p)
end

fast_pair!(reads)
main(State())
