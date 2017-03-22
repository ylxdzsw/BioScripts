include("/home/zhangsw/falcon/Falcon.jl")
include("/haplox/users/zhangsw/ml/Keras.jl")

using OhMyJulia
using Fire
using HDF5
using BinBlocks
using StatsBase
using BioDataStructures
using Falcon

const model = let
    input = Keras.Input(shape=(256, 64, 8))

    output = input |>
        Keras.Convolution2D(64, 5, 5, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((2, 2), strides=(2, 2)) |>

        Keras.Convolution2D(128, 3, 3, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((2, 2), strides=(2, 2)) |>

        Keras.Convolution2D(256, 3, 3, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((4, 1), strides=(4, 1)) |>

        Keras.Convolution2D(256, 3, 3, activation="relu", border_mode="same") |>
        Keras.Convolution2D(256, 3, 3, activation="relu", border_mode="same") |>
        Keras.MaxPooling2D((2, 2), strides=(2, 2)) |>

        Keras.Flatten() |>
        Keras.Dense(1024, activation="sigmoid") |>
        Keras.Dense(1, activation="sigmoid")

    Keras.Model(input=input, output=output)
end

info(x) = STDERR << now() << ":\t" << x << '\n' << flush

baseid(x) = x == Byte('A') ? 1 : x == Byte('T') ? 2 : x == Byte('C') ? 3 : x == Byte('G') ? 4 : 0

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

function get_gdna(fpileup)
    x = Dict{String, Tuple{f64, Int}}()
    for line in eachline(fpileup)
        line = split(line)
        depth = parse(Int, line[4])

        if depth > 80
            freq = most_significant_mut_count(line[5]) / depth

            if (freq < .001 && rand() < .1) || freq > .24 || (depth > 200 && freq > .2)
                x[string(line[1], ':', line[2])] = freq, depth
            end
        end
    end
    x
end

function parse_vcf_line(line)
    chr, pos, ssc, info = split(line)[[1,2,end-3,end]]
    chr = String(chr)
    pos = parse(i32, pos)
    ssc = parse(i32, match(r"SSC=(.*?);", ssc).captures |> car)
    geno, freq, mrbam = split(info, ':')[[1, 6, end]]
    geno = String(geno)
    freq = parse(f64, freq[1:end-1]) / 100
    mrbam = map(x->parse(Int, x), split(mrbam, ','))
    uref, ualt = sum(mrbam[1:6]), sum(mrbam[7:12])
    udp = uref + ualt
    uaf = ualt / udp
    chr, pos, ssc, geno, freq, uaf, udp
end

function make_image(reads, pos, ref)
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
                                    read.seq[readpos] != ref[refpos]))

                        if refpos == pos
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

                        if refpos == pos
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
    image
end

@main function main(fbam, fpileup, fvcf, fref, fweight)
    info("getting dna")
    gdna = get_gdna(fpileup)

    info("parsing bam")
    bam = Bam(fbam)

    info("indexing bam")
    index = get_index(bam)

    info("loading model")
    model[:load_weights](fweight)

    info("loading ref")
    ref = h5open(fref) |> read

    info("start evaluation")

    i, bufx, bufy = 1, Array{f64}(64, 256, 64, 8), Array{Tuple}(64)

    for line in eachline(fvcf)
        chr, pos, ssc, geno, freq, uaf, udp = parse_vcf_line(line)

        udp < 80 && continue
        gaf, gdp = "$chr:$pos" in keys(gdna) ? gdna["$chr:$pos"] : continue

        reads = let idx = index[chr][pos]
            if length(idx) < 80
                continue
            elseif length(idx) > 256
                idx = sample(idx, 256, replace=false) |> sort
            end
            bam.reads[idx]
        end

        bufx[i, :] = make_image(reads, pos, ref[chr])
        bufy[i] = chr, pos, ssc, geno, freq, uaf, udp, gaf, gdp

        if i == 64
            pred = model[:predict](bufx)
            for i in 1:64
                prt(bufy[i]..., pred[i])
            end
            i = 1
        else
            i += 1
        end
    end

    if i != 1 # remaining
        i -= 1
        bufx = bufx[1:i, :, :, :]
        bufy = bufy[1:i]
        pred = model[:predict](bufx, batch_size=1)

        for i in 1:i
            prt(bufy[i]..., pred[i])
        end
    end
end
