#!/usr/bin/env julia

using OhMyJulia
using Insane
using HDF5
using StatsBase
import Base: start, next, done, iteratorsize, eltype,
             getindex, setindex!, show, ==, hash, write

include("/home/zhangsw/falcon/mut.jl")
include("/home/zhangsw/falcon/read.jl")
include("/home/zhangsw/falcon/bam.jl")
include("/home/zhangsw/falcon/sam.jl")
include("/home/zhangsw/falcon/pileup.jl")

const fbam, fpileup, fout = ARGS

function filter_pileup_seq(x)
    i, buf = 1, IOBuffer()

    while i < length(x)
        c = uppercase(x[i])

        if c in "ATCG"
            buf << c
            i += 1
        elseif c == '^'
            i += 2
        elseif c == '+' || c == '-'
            buf << c
            len, i = parse(x, i+1, greedy=false)
            i += len
        else # '.' or ',' or '$' or 'N'
            i += 1
        end
    end

    takebuf_array(buf)
end

const bam = Bam(open(fbam))
const ref = h5open("/haplox/users/zhangsw/hg19.h5")
const out = h5open(fout, "w")
STDERR << now() << " - gdna started" << '\n' << flush
const gdna = @with Dict{String, Tuple{Int, Bytes}}() do x
    for line in eachline(fpileup)
        line = split(line)
        depth = parse(Int, line[4])
        if depth > 50
            x[line[1] * line[2]] = depth, filter_pileup_seq(line[5])
        end
    end
end
STDERR << now() << " - gdna finished" << '\n' << flush
baseid(x) = x == Byte('A') ? 1 : x == Byte('T') ? 2 : x == Byte('C') ? 3 : x == Byte('G') ? 4 : 0

#= NOTE: both falcon and mpileup using 1-based indexing =#

current_chr, current_ref, current_pos = -2, Bytes(), -2
bufx, bufy = Array{f32}(16, 64, 256, 8), Vector{f32}(16)
buf_idx, buf_batch = 1, 1

for (reads, chr, mut) in @task pileup(bam)
    if chr != current_chr
        STDERR << now() << " - $chr started" << '\n' << flush
        current_chr = chr
        current_ref = read(ref, car(bam.refs[chr+1]))
    elseif mut.pos == current_pos
        continue
    end

    current_pos = mut.pos

    pos, base = if isa(mut, SNP)
        mut.pos, mut.alt
    elseif isa(mut, Insertion)
        mut.pos-1, Byte('+')
    elseif isa(mut, Deletion)
        mut.pos, Byte('-')
    end

    label = let
        key = car(bam.refs[chr+1])*string(mut.pos)
        key in keys(gdna) || continue

        depth, seq = gdna[key]
        freq = count(x->x==base, seq) / depth
        if freq < .001
            0.
        elseif freq > .35
            1.
        else
            continue
        end
    end

    reads = collect(reads)
    if length(reads) < 50
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

    bufx[buf_idx, :] = image
    bufy[buf_idx]    = label

    if buf_idx == 16
        out[string('x', buf_batch), "compress", 8] = bufx
        out[string('y', buf_batch), "compress", 8] = bufy
        buf_idx = 1
        buf_batch += 1
    else
        buf_idx += 1
    end
end
