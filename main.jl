var  = open(readlines, ARGS[1])
size = open(readlines, ARGS[2])

dict = Dict{ASCIIString, Dict{ASCIIString, ASCIIString}}()

for line in var
    name, _, t = split(line, '.')
    num = split(line, '\t')[2] |> Base.chomp!
    name in keys(dict) || (dict[name] = Dict{ASCIIString, ASCIIString}())
    dict[name][t] = num
end

for line in size
    s, t, name, num = match(r"^(.*)\t(gDNA|cfDNA)/(.*)_(1|2).*$", line).captures
    name in keys(dict) || (dict[name] = Dict{ASCIIString, ASCIIString}())
    dict[name][t*num] = s
end

@printf("%14s %6s %5s %5s %5s %5s %5s\n",
        "name","snp","indel","cf1","cf2","g1","g2")

for (n,v) in dict
    try
        @printf("%14s %6s %5s %5s %5s %5s %5s\n",
                n, 
                v["snp"], 
                v["indel"], 
                v["cfDNA1"], 
                v["cfDNA2"], 
                v["gDNA1"], 
                v["gDNA2"])
    catch
        continue
    end
end