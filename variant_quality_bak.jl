using PyCall
using OhMyJulia
using JSON

@pyimport pysam

fvcf, fbam = ARGS

sam = pysam.AlignmentFile(fbam)

fout = open("fuck.html", "w")
fout << """<html><head><script>"""
fout << open(read, "/home/zhangsw/.julia/v0.5/PlotlyJS/deps/plotly-latest.min.js")
fout << """;window.PLOTLYENV={BASE_URL:"https://plot.ly"}"""
fout << """</script></head><body>"""

for line in eachline(fvcf)
    chr, pos, ref, alt, info = split(line, '\t')[[1, 2, 4, 5, end]]

    let
        rd, ad = split(info, ':')[[5, 6]] ~ map(x->parse(Int, x))

        if rd < 1000 || ad < 10
            continue
        end
    end

    pos = parse(Int, pos) - 1

    rl, al, ol = Int[], Int[], Int[]

    for r in sam[:fetch](chr, pos, pos+1)
        aligned_pairs = r[:get_aligned_pairs](matches_only=true)

        query_pos = let
            i = findfirst(x->cadr(x)==pos, aligned_pairs)
            i == 0 && continue
            aligned_pairs[i] |> car
        end

        base, qual = r[:query_sequence][query_pos+1], r[:query_qualities][query_pos+1] # julia use 1-base index

        if base in ref
            push!(rl, qual)
        elseif base in alt
            push!(al, qual)
        else
            push!(ol, qual)
        end
    end

    fout << """<div id="$chr-$pos" class="plotly-graph-div"></div>"""
    fout << """<script>
        Plotly.newPlot("$chr-$pos", [
            { y: $(json(rl)), type: "box", name: "ref" },
            { y: $(json(al)), type: "box", name: "alt" },
            { y: $(json(ol)), type: "box", name: "other" }
        ], {
            title: "$chr:$(pos+1)",
            margin: {r:80, l:80, b:80, t:100}
        },{
            showLink: false
        })
    </script>"""
end

fout << """</body></html>""" << close
