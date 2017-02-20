using OhMyJulia
using HDF5
const db = h5open("hg19.h5", "w")

chr, buffer = "", IOBuffer()

for line in eachline("hg19.fa")
    if car(line) == '>'
        if chr != ""
            write(db, chr, takebuf_array(buffer))
        end

        chr = line[2:end-1]
        println("importing $chr")
    else
        write(buffer, uppercase(chomp(line)))
    end
end

write(db, chr, takebuf_array(buffer))
