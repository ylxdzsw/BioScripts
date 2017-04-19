#=
/thinker/dstore/rawfq/160926_NS500713_0075_AHJTF2BGXY/S005_20160907208-2_cfdna_pan-cancer-v1_R1_001.fastq.gz
s005 - s009
=#

using OhMyJulia
using Libz

fq  = ZlibInflateInputStream(STDIN)
out = ZlibDeflateOutputStream(STDOUT)

for (i, line) in enumerate(eachline(fq))
    if i % 4 == 1
        sp   = findfirst(line, ' ')
        id   = SubString(line, 1, sp-1)
        info = SubString(line, sp)
        code = SubString(line, findlast(line, '+')+1, endof(line)-1)
        out << id << ':' << code << info
    else
        out << line
    end
end
