#=
    for x in $(du -sh /haplox/users/huanghuiqiang/project/tenThousantProject/breastProject_20160812/cf_result/*/*_sort.bam | grep -v G | cut -f 2); do
        gzip -cd $x | julia count.jl > $(basename "$x" .bam).mut &;
    done
=#

using OhMyJulia
using Insane
import Base: start, next, done, iteratorsize, eltype,
             getindex, setindex!, show, ==, hash

include(rel"../falcon/mut.jl")
include(rel"../falcon/read.jl")
include(rel"../falcon/bam.jl")
include(rel"../falcon/sam.jl")
include(rel"../falcon/pair.jl")

bam = Bam(STDIN)
reads = collect(bam)
fast_pair!(reads)

N = 0
M = Dict((x => y) => 0 for x in b"ATCG" for y in b"ATCG" if x != y)

read_base(r::Read, m::Mut) = r.seq[car(calc_read_pos(r, m.pos))]

for r in reads @when isdefined(r, :mate) && r.refID == r.mate.refID && (r.pos < r.mate.pos || (r.pos == r.mate.pos && pointer_from_objref(r) < pointer_from_objref(r.mate)))
    offset, status = calc_read_pos(r, r.mate.pos)
    status == 0x00 || continue
    len = min(length(r.seq) - offset + 1, length(r.mate.seq))
    any((isa(m, Insertion) || isa(m, Deletion)) && m.pos >= offset for m in r.muts) && continue
    any((isa(m, Insertion) || isa(m, Deletion)) && m.pos <= len for m in r.mate.muts) && continue
    N += 2len
    for m in r.muts @when offset <= m.pos <= offset + len - 1 && r.mate.seq[m.pos - offset + 1] == m.ref && m.alt != Byte('N')
        M[m.ref => m.alt] += 1
    end
    for m in r.mate.muts @when m.pos <= len && r.seq[m.pos + offset - 1] == m.ref && m.alt != Byte('N')
        M[m.ref => m.alt] += 1
    end
end

for (k, v) in M
    STDOUT << k.first << "=>" << k.second << ':' << round(100v / N, 4) << "%\n"
end