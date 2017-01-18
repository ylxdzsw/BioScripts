using JSON: json

samples, data = let
    lines = open(readlines, "depth.txt")

    samples = lines ~ car ~ split ~ cdr ~ map(String)

    data = map(cdr(lines)) do i
        i = i ~ chomp ~ :split('\t')
        car(i), map(parse, cdr(i))
    end

    samples, data
end

open("fuck.html", "w") do fout
    fout << """<html><head><script>"""
    fout << open(read, s"C:\Users\ylxdzsw\.julia\v0.5\PlotlyJS\deps\plotly-latest.min.js")
    fout << """;window.PLOTLYENV={BASE_URL:"https://plot.ly"}"""
    fout << """;window.data={samples:$(json(samples))}"""
    fout << """</script></head><body>"""

    fout << """<div id="heatmap" class="plotly-graph-div"></div>"""
    fout << """<script>
        Plotly.newPlot("heatmap", [{
            x: window.data.samples,
            y: $(json(map(car, data))),
            z: $(json(map(x->map(x->x>8000?8000:x,cadr(x)), data))),
            type: "heatmap"
        }], {
            title: "heatmap",
            height: 1080,
            margin: {r:80, l:120, b:80, t:80}
        },{
            showLink: false
        })
    </script>"""

    fout << """<script>
        window.data.scatter = {
            x: window.data.samples,
            y: $(json(map(mean, zip(map(cadr, data)...)))),
            type: "scatter",
            mode: "markers"
        }
    </script>"""

    for i in data
        chr, pos = i ~ car ~ :split(':')

        fout << """<div id="$chr-$pos" class="plotly-graph-div"></div>"""
        fout << """<script>
            Plotly.newPlot("$chr-$pos", [
                window.data.scatter,
                {
                    x: window.data.samples,
                    y: $(json(cadr(i))),
                    type: "bar"
                }
            ], {
                title: "$chr:$pos",
                margin: {r:80, l:80, b:80, t:100}
            },{
                showLink: false
            })
        </script>"""
    end

    fout << """</body></html>"""
end
