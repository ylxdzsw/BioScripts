#!/usr/bin/env julia
using OhMyJulia

if isempty(ARGS)
    println("""\n
        example:
        \$ setsid julia ./intersection.jl batch?/00005.snp_MrBam.txt.newTxt > 00005.intersection.txt
    """)
    exit(0)
end

N = length(ARGS)
appear = Dict{String, Int}()

for i = 1:N
    STDOUT << ('@'+i) << ": " << ARGS[i] << '\n'
    for l in eachline(ARGS[i])
        d = findnext(l, '\t', 1)
        d = findnext(l, '\t', d+1)
        d = findnext(l, '\t', d+1)
        d = findnext(l, '\t', d+1)
        l = l[1:d-1]
        appear[l] = get(appear, l, 0) | (1 << i >> 1)
    end
end

counts = groupby(identity, (x,y)->x+1, 0, values(appear))

for combination = 1:1<<N-1
    for i = 1:N
        if combination & (1<<i>>1) != 0
            STDOUT << ('@'+i) << "  "
        else
            STDOUT << "   "
        end
    end
    STDOUT << sum(c for (a,c) in counts if a & combination == combination) << '\n'
end
