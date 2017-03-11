#!/usr/bin/env julia

using OhMyJulia
using Insane
using StatsBase
using HDF5
using BinBlocks
import Base: start, next, done, iteratorsize, eltype,
             getindex, setindex!, show, ==, hash, write

include("/home/zhangsw/falcon/mut.jl")
include("/home/zhangsw/falcon/read.jl")
include("/home/zhangsw/falcon/bam.jl")
include("/home/zhangsw/falcon/sam.jl")
include("/home/zhangsw/falcon/pileup.jl")

const fbam, fpileup, fout = ARGS

function most_significant_mut_count(seq)
    i, counter = 1, fill(0, 6) # ATCGID

    while i < length(seq)
        c = uppercase(seq[i])

        m = findfirst("ATCG+-", c)

        if m > 4
            counter[m] += 1
            len, i = parse(seq, i+1, greedy=false)
            i += len
        elseif m > 0
            counter[m] += 1
            i += 1
        elseif c == '^'
            i += 2
        else
            i += 1
        end
    end

    maximum(counter)
end

const bam = Bam(open(fbam))
const ref = h5open("/haplox/users/zhangsw/hg19.h5")
const out = BinBlock(fout, "w")
STDERR << now() << " - gdna started" << '\n' << flush
const gdna = @with Dict{String, f32}() do x
    for line in eachline(fpileup)
        line = split(line)
        depth = parse(Int, line[4])

        if depth > 80
            freq = most_significant_mut_count(line[5]) / depth

            if freq < .001 && rand() < .01
                x[string(line[1], ':', line[2])] = 0.
            elseif freq > .24 || (depth > 200 && freq > .2)
                x[string(line[1], ':', line[2])] = 1.
            end
        end
    end
end
baseid(x) = x == Byte('A') ? 1 : x == Byte('T') ? 2 : x == Byte('C') ? 3 : x == Byte('G') ? 4 : 0

#= NOTE: both falcon and mpileup using 1-based indexing =#

current_chr, current_ref, current_pos = -2, Bytes(), -2

for (reads, chr, mut) in @task pileup(bam)
    if chr != current_chr
        STDERR << now() << " - $chr started" << '\n' << flush
        current_chr = chr
        current_ref = read(ref, car(bam.refs[chr+1]))
    end

    pos = isa(mut, Insertion) ? mut.pos-1 : mut.pos

    pos == current_pos && continue

    current_pos = pos

    label = let
        key = string(car(bam.refs[chr+1]), ':', mut.pos)
        if key in keys(gdna)
            gdna[key]
        else
            continue
        end
    end

    reads = collect(reads)
    if length(reads) < 80
        continue
    elseif length(reads) > 256
        reads = sample(reads, 256, replace=false, ordered=true)
    end

    image = fill(f32(0), 256, 64, 8)

    for (i, read) in enumerate(reads)
        seq, offset = let

            seq, readpos, refpos = [], 1, read.pos

            for cigar in read.cigar
                op, len = cigar & 0x0f, Int(cigar >> 4)
                if op == 0
                    while len > 0
                        push!(seq, (baseid(read.seq[readpos]),
                                    min(read.qual[readpos], 60) / 60,
                                    read.seq[readpos] != current_ref[refpos]))

                        if refpos == mut.pos
                            offset = length(seq)
                        end

                        refpos  += 1
                        readpos += 1
                        len     -= 1
                    end
                elseif op == 1
                    if !isempty(seq)
                        last = seq[end]
                        seq[end] = (last..., min(len, 20) / 20)
                    end
                    readpos += len
                elseif op == 2
                    while len > 0
                        push!(seq, (6, 1., false))

                        if refpos == mut.pos
                            offset = length(seq)
                        end

                        refpos  += 1
                        len     -= 1
                    end
                elseif op == 4
                    readpos += len
                elseif op == 5
                    # pass
                else
                    error("TODO")
                end
            end

            seq, offset
        end

        offset -= 31

        for j in max(1, 1-offset):min(64, length(seq) - offset)
            p = seq[j + offset]

            if car(p) != 0
                image[i, j, car(p)] = cadr(p)
                image[i, j, 7] = p[3]
            end

            if length(p) == 4
                image[i, j, 5] = p[4]
            end
        end

        image[i, :, 8] = read.flag & 0x0010 != 0
    end

    image[length(reads)+1:256, :, 8] = .5
    write(out, image, label)
end
