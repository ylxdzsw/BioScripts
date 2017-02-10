using OhMyJulia
using JSON: json
using HDF5

const X, y = h5open("data.hdf5", "r") do f
    read(f, "X"), read(f, "y")
end

const sample = rand(1:length(y), 5000)
const d = X[sample, end-10] .+ X[sample, end-11]
const f = X[sample, end-9]
const s = map(Bool, y[sample])

open("fuck.html", "w") do fout
    fout << """<html><head><script>"""
    fout << open(read, s"C:\Users\ylxdzsw\.julia\v0.5\PlotlyJS\deps\plotly-latest.min.js")
    fout << """</script></head><body>"""

    fout << """<div id="cutoff" class="plotly-graph-div"></div>"""
    fout << """<script>
        Plotly.newPlot("cutoff", [{
            x: $(json(min(d[s], 5000))),
            y: $(json(f[s])),
            name: "somatic",
            type: "scatter",
            mode: 'markers',
            marker: { size: 4 }
        }, {
            x: $(json(min(d[!s], 5000))),
            y: $(json(f[!s])),
            name: "germline",
            type: "scatter",
            mode: 'markers',
            marker: { size: 4 }
        }], {
            title: "dist. of depth & freq",
            height: 600,
            margin: {r:80, l:120, b:80, t:80},
            xaxis: {
                tickvals: [0,716,1433,2149,2866,3582,4298]
            },
            shapes: [
                { type: 'line', x0: 0, x1: 716, y0: 0.21223021582733814, y1: 0.21223021582733814 },
                { type: 'line', x0: 716, x1: 1433, y0: 0.08299319727891157, y1: 0.08299319727891157 },
                { type: 'line', x0: 1433, x1: 2149, y0: 0.051713395638629284, y1: 0.051713395638629284 },
                { type: 'line', x0: 2149, x1: 2866, y0: 0.08563194730341704, y1: 0.08563194730341704 },
                { type: 'line', x0: 2866, x1: 3582, y0: 0.0978225544361391, y1: 0.0978225544361391 },
                { type: 'line', x0: 3582, x1: 4298, y0: 0.08480086697371986, y1: 0.08480086697371986 },
                { type: 'line', x0: 4298, x1: 5000, y0: 0.0951362907536077, y1: 0.0951362907536077 }
            ]
        })
    </script>"""

    fout << """</body></html>"""
end







