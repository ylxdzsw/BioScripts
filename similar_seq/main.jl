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

const base_code = b"ATCGN"

function encode_base(b)::Byte
    b == Byte('A') ? 0x00 :
    b == Byte('T') ? 0x01 :
    b == Byte('C') ? 0x02 :
    b == Byte('G') ? 0x03 : 0x04
end

function decode_base(b)::Byte
    base_code[b+1]
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

@main function split_bed(bed)
    chr, anchors = load_data()
    f, n, batch = open("batch_00.task", "w"), 1, 0
    for line in eachline(bed)
        chr, a, b = split(line, '\t')

        for i in parse(Int, a)-128:64:parse(Int, b)-63
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

@main function rebatch()
    i, c, f = 0, 0, open("rebatch_000.txt", "w")
    for line in eachline(STDIN)
        write(f, line)
        c += length(split(line, '\t')) - 2

        if c > 2048
            f = open(@sprintf("rebatch_%03d.txt", i += 1), "w")
            c = 0
        end
    end
end

@main function report_alignments(batch)
    ref, anchors = load_data()
    ref_rv       = reverse_complement(ref)

    F, P = Matrix{Int}(256, 65536+255), Matrix{Byte}(256, 65536+255)

    for line in eachline(batch)
        bpos, forward, reverse = parse_batch_line(line)
        q = view(ref, bpos:bpos+255)

        for (tp, tr, rev) in ((forward, ref, false), (reverse, ref_rv, true)), block in tp
            t = view(tr, 65536(block-1)+1:65536block+255)

            dp_fill!(q, t, F, P)
            for (pq, pt, var) in dp_report(q, t, F, P, 64, decode_base)
                gpq = bpos + pq - 1
                gpt = rev ? length(ref) - 65536(block-1) - pt + 1 : 65536(block-1) + pt

                !rev && gpq == gpt && continue # source

                try
                    println(join(block_to_origin(gpq, anchors), ':'), " ~ ",
                            join(block_to_origin(gpt, anchors), ':'),
                            rev ? " (reverse)\n" : "\n", var)
                catch
                    prt(STDERR, gpq, gpt, rev, var)
                end
            end
        end
    end

    flush(STDOUT)
    run(`rm $batch`)
end

@main function fucking_region(bed, alignments...)
    bed = let r = Dict{String, IntRangeSet{i64}}()
        for (chr, ps, pe) in eachline(split, bed)
            chr in keys(r) || (r[chr] = IntRangeSet{i64}())
            push!(r[chr], parse(Int, ps)+1:parse(Int, pe))
        end
        r
    end

    dict = Dict("chr$x"=>IntRangeSet{i64}() for x in (1:22...,:X,:Y))
    local cl, cr, pl, pr, rev, ll, lr
    for file in alignments, (i, line) in enumerate(eachline(chomp, file))
        if i % 5 == 1
            rev = contains(line, "reverse")
            left, _, right = split(line, ' ')
            cl, pl = split(left, ":")
            cr, pr = split(right, ":")
            pl, pr = parse(Int, pl), parse(Int, pr)
        elseif i % 5 == 2
            ll = count(x->x!='-', line)
        elseif i % 5 == 4
            lr = count(x->x!='-', line)
        elseif i % 5 == 0
            (ll < 50 || lr < 50) && continue
            push!(dict[cl], pl:pl+ll-1)
            push!(dict[cr], rev ? (pl-ll+1:pl) : (pl:pl+ll-1))
        end
    end

    for chr in intersect(keys(bed), keys(dict))
        foreach(intersect(bed[chr], dict[chr])) do x
            prt(chr, x.start-1, x.stop)
        end
    end

    # for (k, v) in dict
    #     foreach(v) do x
    #         prt(k, x.start-1, x.stop)
    #     end
    # end
end

"use `sort -u -k1,4 file` to get unique variants"
@main function find_variants(alignments...)
    local chr, pos, q, a, t
    for file in alignments, (i, line) in enumerate(eachline(chomp, file))
        if i % 5 == 1
            chr, pos = split(split(line, ' ')[1], ':')
            pos = parse(Int, pos)
        elseif i % 5 == 2
            q = line
        elseif i % 5 == 3
            a = line
        elseif i % 5 == 4
            t = line
        elseif i % 5 == 0
            p, l = 1, 0
            for s in 1:length(a)
                if a[s] == '|' # match
                    if l > 0
                        prt(chr, p, q[s-l:s-1], '-', "$(file[9:11]):$(i-4)")
                    elseif l < 0
                        prt(chr, p, '-', t[s+l:s-1], "$(file[9:11]):$(i-4)")
                    end
                    pos += 1
                    l = 0
                elseif q[s] == '-' # insertion
                    if l > 0
                        warn("deletion immediately followed by insertion")
                        prt(chr, p, q[s-l:s-1], '-', "$(file[9:11]):$(i-4)")
                        p = pos - 1
                        l = -1
                    elseif l == 0
                        p = pos - 1
                        l = -1
                    elseif l < 0
                        l -= 1
                    end
                elseif t[s] == '-' # deletion
                    if l > 0
                        l += 1
                    elseif l == 0
                        p = pos
                        l = 1
                    elseif l < 0
                        warn("insertion immediately followed by deletion")
                        prt(chr, p, '-', t[s+l:s-1], "$(file[9:11]):$(i-4)")
                        p = pos
                        l = 1
                    end
                    pos += 1
                else # snp
                    if l > 0
                        prt(chr, p, q[s-l:s-1], '-', "$(file[9:11]):$(i-4)")
                    elseif l < 0
                        prt(chr, p, '-', t[s+l:s-1], "$(file[9:11]):$(i-4)")
                    end
                    prt(chr, pos, q[s], t[s], "$(file[9:11]):$(i-4)")
                    pos += 1
                    l = 0
                end
            end
        end
    end
end

@main function find_fusion(alignments...)
    local chr, pos, l
    for file in alignments, (i, line) in enumerate(eachline(chomp, file))
        if i % 5 == 1
            chr, pos = split(split(line, ' ')[1], ':')
            pos = parse(Int, pos)
        elseif i % 5 == 2
            l = count(x->x!='-', line)
        elseif i % 5 == 3
            if count(x->x=='|', line) / length(line) > .98
                prt(chr, pos)
                prt(chr, pos + l - 1)
            end
        end
    end
end
