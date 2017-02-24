using OhMyJulia
using DataFrames

files = ("stat_driver.txt", "stat_new.txt", "stat_sumed_driver.txt", "stat_sumed_new.txt")

data = DataFrame(id=readtable("fjj_stat_new.csv")[:, :id])

for i in ("stat_driver.txt", "stat_new.txt", "stat_sumed_driver.txt", "stat_sumed_new.txt")
    name = Symbol(i[6:end-4])
    data[name] = 0
    for line in eachline(i) @when startswith(line, "Cross-Validation Outlyzers ID: ")
        for sample in split(line[31:end-1], ", ")
            data[data[:id].==parse(Int, sample), name] += 1
        end
    end
end

writetable("table.csv", data)