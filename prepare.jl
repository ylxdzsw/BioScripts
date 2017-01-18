import JSON
using PlotlyJS

a = JSON.parsefile("cache.json")

cf = map(Bool, a["label"])

data = a["data"][cf]

g = map(1:4) do g
	d = [map(x->x[g+4i], data) for i in 0:9]
	d = map(d) do x x[0.05 .< x .< 0.95] end
	l = length(d[1])
	d = map(d) do x sort(x)[ceil(Int,.02l):ceil(Int,.98l)] end
end

open("processed.json", "w") do x
	write(x, JSON.json(g))
end