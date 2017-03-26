using OhMyJulia
using BioDataStructures
using Fire
using HDF5

include("patch.jl")
include("align.jl")

function find_valid_interval(seq::Bytes)
    s = IntRangeSet{Int}()

    for (i, b) in enumerate(seq)
        b != Byte('N') && push!(s, i)
    end

    s
end

function origin_to_block(chr, pos, anchors)
    acc = 0
    for (c, r) in anchors
        if c == chr && pos in r
            return acc + pos - r.start + 1
        else
            acc += length(r) + 16
        end
    end
    error("fuck")
end

function block_to_origin(ind, anchors)
    for (c, r) in anchors
        if ind > length(r)
            ind -= length(r) + 16
        else
            return c, r.start + ind - 1
        end
    end
    error("fuck")
end

function encode_base(b)::Byte
    b == Byte('A') ? 0x00 :
    b == Byte('T') ? 0x01 :
    b == Byte('C') ? 0x02 :
    b == Byte('G') ? 0x03 : 0x04
end

function load_data()
    ref = h5read("cache.h5", "ref")
    anchors = h5read("cache.h5", "anchors")
    anchors = map(x->split(x, '\t'), split(anchors, '\n'))
    anchors = map(x->(String(x[1]), parse(Int, x[2]):parse(Int, x[3])), anchors)
    ref, anchors
end

function findsim(query, ref, batch_size)
    r = Array{Byte}(batch_size)
    @tcall((:align, "ksw.so"), Void,
           (Ptr{Byte}, Ptr{Byte}, Ptr{Byte}, Cint),
           query, ref, r, batch_size)
    find(r .>= 0x40)
end

function reverse_complement(seq)
    result, N = similar(seq), length(seq)

    for i in 1:N
        b = seq[N-i+1]
        result[i] = b == 0x00 ? 0x01 :
                    b == 0x01 ? 0x00 :
                    b == 0x02 ? 0x03 :
                    b == 0x03 ? 0x02 : 0x04
    end

    result
end

function parse_batch_line(line)
    line = split(line)
    bpos = parse(Int, car(line))
    list = cdr(line)
    star = findfirst(list, "*")
    forward = map(x->parse(Int, x), list[1:star-1])
    reverse = map(x->parse(Int, x), list[star+1:end])
    bpos, forward, reverse
end

@main function split_bed(bed)
    chr, anchors = load_data()
    f, n, batch = open("batch_00.task", "w"), 1, 0
    for line in eachline(bed)
        chr, a, b = split(line, '\t')

        for i in parse(Int, a)-63:64:parse(Int, b)+64
            write(f, origin_to_block(chr, i, anchors))

            if n == 2048
                close(f)
                batch += 1
                n = 1
                f = open(@sprintf("batch_%02d.task", batch), "w")
            else
                n += 1
            end
        end
    end
end

@main function prepare_data()
    hg19 = h5read("/haplox/users/zhangsw/hg19.h5", "/")
    ref, anchors = [], []
    for chr in keys(hg19)
        foreach(find_valid_interval(hg19[chr])) do x
            push!(anchors, (chr, x))
            push!(ref, map(encode_base, hg19[chr][x]), fill(0x04, 16))
        end
    end
    ref = Byte[ref...;]
    pad = 65536 - length(ref) % 65536 + 255
    info("pad = $pad")
    append!(ref, fill(0x04, pad))
    h5write("cache.h5", "ref", ref)
    h5write("cache.h5", "anchors", join(map(x->"$(x[1])\t$(x[2].start)\t$(x[2].stop)", anchors), '\n'))
end

@main function detect_similarity()
    ref, anchors = load_data()
    ref_rv       = reverse_complement(ref)
    batch_size   = (length(ref) - 255) ÷ 65536
    limit, cond  = Ref(128), Condition()
    for batch in readdir(".") @when startswith(batch, "batch") && endswith(batch, ".task")
        fin, fout = open(batch), open(car(splitext(batch)) * ".txt", "w")

        while !eof(fin)
            limit[] == 0 && wait(cond)
            limit[] -= 1
            let pos = fin >> Int
                @schedule begin
                    query = pointer(ref, pos)
                    prt(fout, pos, findsim(query, ref, batch_size)...,
                        '*', findsim(query, ref_rv, batch_size)...)
                    limit[] += 1
                    notify(cond)
                end
            end
        end

        while limit[] != 128 wait(cond) end

        close(fin)
        close(fout)
        STDERR << now() << '\t' << batch << " done" << '\n' << flush
        run(`rm $batch`)
    end
end

@main function report_alignments(batch)
    ref, anchors = load_data()
    ref_rv       = reverse_complement(ref)

    F, P = Matrix{Int}(256, 65536+255), Matrix{Byte}(256, 65536+255)

    for line in eachline(batch)
        bpos, forward, reverse = parse_batch_line(line)
        for (tp, tr, rev) in ((forward, ref, false), (reverse, ref_rv, true)), block in tp
            q = view(ref, bpos:bpos+255)
            t = view(ref, 65536(block-1)+1:65536block+255)

            for (pq, pt, var) in dp_report(q, t, F, P)
                gpq = qpos + pq - 1
                gpt = rev ? length(ref) - 65536(block-1) + pt - 1 : 65536(block-1) + pt

                !rev && gpq == gpt && continue # source

                println(join(block_to_origin(gpq, anchors), ':'), " ~ ",
                        join(block_to_origin(gpt, anchors), ':'),
                        rev ? " (reverse)\n" : "\n", var)
            end
        end
    end
end
