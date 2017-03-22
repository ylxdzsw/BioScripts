using OhMyJulia
using HDF5
using Fire

@main function main(fasta, hdf5)
    db = h5open(hdf5, "w")
    chr, buffer = "", IOBuffer()

    for line in eachline(fasta)
        if car(line) == '>'
            if chr != ""
                write(db, chr, takebuf_array(buffer))
            end

            chr = split(line[2:end-1], ' ') |> car |> String
            println("importing $chr")
        else
            write(buffer, uppercase(chomp(line)))
        end
    end

    write(db, chr, takebuf_array(buffer))
end
