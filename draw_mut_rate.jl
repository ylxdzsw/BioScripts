using OhMyJulia
using JSON.json

function with_template(f::Function)
    STDOUT << open(read, rel"../falcon/template.html") << f << "</body></html>"
end

const M = Dict((x => y) => Tuple{String, Float64}[] for x in "ATCG" for y in "ATCG" if x != y)

const color = Dict('A' => "red", 'T' => "purple", 'C' => "blue", 'G' => "green")

for file in ARGS
    for line in eachline(file)
        push!(M[line[1]=>line[4]], (car(split(file, '_')), parse(line[6:end-2])))
    end
end

with_template() do f
    f << """<section class="scen">"""
    for ref = "ATCG"
        a,b,c = filter(x->x!=ref, "ATCG")
        f << """<div id="$ref" class="plotly-graph-div" style="float: left; width: 49%"></div>"""
        f << """<script>
            var a = $(json(map(cadr, M[ref=>a])))
            var b = $(json(map(cadr, M[ref=>b])))
            var c = $(json(map(cadr, M[ref=>c])))
            Plotly.newPlot("$ref", [
                { x: $(json(map(car, M[ref=>a]))), y: a, name: "$a", line: {color: "$(color[a])"} },
                { x: $(json(map(car, M[ref=>b]))), y: b, name: "$b", line: {color: "$(color[b])"} },
                { x: $(json(map(car, M[ref=>c]))), y: c, name: "$c", line: {color: "$(color[c])"} },
                { y: a, name: "$a", type: "box", marker: {color: "$(color[a])"}, boxmean: true },
                { y: b, name: "$b", type: "box", marker: {color: "$(color[b])"}, boxmean: true },
                { y: c, name: "$c", type: "box", marker: {color: "$(color[c])"}, boxmean: true }
            ], {
                margin: { r:80, l:80, b:80, t:100 },
                title: "$ref",
                xaxis: { title: "sample" },
                yaxis: { title: "frequency (%)" }
            })
        </script>"""
    end
    f << """</section>"""
end
