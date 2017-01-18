using OhMyJulia

if length(ARGS) < 2
    println(s"""
        $ setsid find_Kmer hg19.fa 4 3 3 > Kmers.bed
        表示输出重复至少4次的单碱基重复序列，重复至少3次的双碱基重复序列，重复至少为3次的三碱基序列，以此类推
        3G的hg19要跑大概900分钟
    """)
    exit(0)
end

const L = map(parse, cdr(ARGS))
const f = open(car(ARGS))

type Kmer{K}
    seq::NTuple{K, Byte}
    pos::Int
    len::Int
    next::Vector{Int}
    Kmer(seq) = new(NTuple{K, Byte}(seq), 1, 0, get_next(seq))
end

function get_next(seq)
    j = 1
    next = fill(1, length(seq))

    for i = 2:length(seq)
        while j > 1 && seq[i] != seq[j]
            j = next[j-1]
        end

        if seq[i] == seq[j]
            j += 1
        end

        next[i] = j
    end

    next
end

function feed!{K}(k::Kmer{K}, x::Byte, chr::String, p::Int)
    if k.seq[k.pos] == x
        if k.pos == K
            k.pos = 1
            k.len += 1
        else
            k.pos += 1
        end
    else
        if k.len >= L[K]
            pend = p-k.pos # 1-based, both side included
            println(chr, '\t', pend-k.len*K+1, '\t', pend, '\t', map(Char, k.seq)...)
        end

        if k.pos > 1
            k.pos, k.len = k.next[k.pos-1], 0
            feed!(k, x, chr, p)
        elseif k.len == 1
            k.pos, k.len = k.next[end], 0
            feed!(k, x, chr, p)
        else
            k.pos, k.len = 1, 0
        end
    end
end

cartesian_product(x) = ((x,) for x in x)
cartesian_product(x, y) = ((i...,j) for i in x for j in y)
cartesian_product(x...) = reduce(cartesian_product, x)

Kmers = collect(Kmer{i}(x) for i in 1:length(L) for x in cartesian_product(fill("ATCG", i)...) if i==1 || !all(Δ(==, collect(x))))

p = 0 # 1-based
chr = "bug"

while !eof(f)
    c = f >> Byte

    if c == Byte('>')
        chr, p = readline(f)[1:end-1], 0
        continue
    end

    if c == Byte('\n')
        continue
    end

    p += 1
    c > 0x60 && (c -= 0x20) # capital insencitive

    for i in Kmers
        feed!(i, c, chr, p)
    end
end
