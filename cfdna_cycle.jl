import JSON
using JsonBuilder
using OhMyJulia

data = open(readstring, "cfdna_data.json") |> JSON.parse

A = map(x->x[1:4:40], data)
T = map(x->x[2:4:40], data)
C = map(x->x[3:4:40], data)
G = map(x->x[4:4:40], data)

gA = map(1:10) do i @json """
    {
        x: $i,
        y: $(map(x->x[i], A)),
        type: "box",
        xaxis: "x2",
        yaxis: "y2",
        line: { color: "red" },
        boxpoints: false,
        showlegend: false,
        name: $("cycle$i"),
        x0: $i,
        legendgroup: "A"
    }
""" end

gT = map(1:10) do i @json """
    {
        x: $i,
        y: $(map(x->x[i], T)),
        type: "box",
        xaxis: "x3",
        yaxis: "y3",
        line: { color: "purple" },
        boxpoints: false,
        showlegend: false,
        name: $("cycle$i"),
        x0: $i,
        legendgroup: "T"
    }
""" end

gC = map(1:10) do i @json """
    {
        x: $i,
        y: $(map(x->x[i], C)),
        type: "box",
        xaxis: "x4",
        yaxis: "y4",
        line: { color: "blue" },
        boxpoints: false,
        showlegend: false,
        name: $("cycle$i"),
        x0: $i,
        legendgroup: "C"
    }
""" end

gG = map(1:10) do i @json """
    {
        x: $i,
        y: $(map(x->x[i], G)),
        type: "box",
        xaxis: "x5",
        yaxis: "y5",
        line: { color: "green" },
        boxpoints: false,
        showlegend: false,
        name: $("cycle$i"),
        x0: $i,
        legendgroup: "G"
    }
""" end

gL = [
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], A)), 1:10)), name: "A", line: { color: "gray" }, xaxis: "x2", yaxis: "y2", showlegend: false, legendgroup: "A"}"""),
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], T)), 1:10)), name: "T", line: { color: "gray" }, xaxis: "x3", yaxis: "y3", showlegend: false, legendgroup: "T"}"""),
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], C)), 1:10)), name: "C", line: { color: "gray" }, xaxis: "x4", yaxis: "y4", showlegend: false, legendgroup: "C"}"""),
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], G)), 1:10)), name: "G", line: { color: "gray" }, xaxis: "x5", yaxis: "y5", showlegend: false, legendgroup: "G"}""")
]

gM = [
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], A)), 1:10)), name: "A", line: { color: "red" },    showlegend: true, legendgroup: "A"}"""),
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], T)), 1:10)), name: "T", line: { color: "purple" }, showlegend: true, legendgroup: "T"}"""),
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], C)), 1:10)), name: "C", line: { color: "blue" },   showlegend: true, legendgroup: "C"}"""),
    @json("""{x: $(1:10), y: $(map(x->mean(map(i->i[x], G)), 1:10)), name: "G", line: { color: "green" },  showlegend: true, legendgroup: "G"}""")
]

layout = @json """{
    xaxis:  {anchor: "y1"},
    xaxis2: {anchor: "y2"},
    xaxis3: {anchor: "y3"},
    xaxis4: {anchor: "y4"},
    xaxis5: {anchor: "y5", title: "Cycles", dtick: 1},
    yaxis:  {domain: [0.8, 1.0], title: "Average"},
    yaxis2: {domain: [0.6, 0.8], title: "A"},
    yaxis3: {domain: [0.4, 0.6], title: "T"},
    yaxis4: {domain: [0.2, 0.4], title: "C"},
    yaxis5: {domain: [0.0, 0.2], title: "G"},
    height: 800,
    width: 600
}"""

open("fuck.html", "w") do fout
    fout << """<html><head><script>"""
    fout << open(read, s"C:\Users\ylxdzsw\.julia\v0.5\PlotlyJS\deps\plotly-latest.min.js")
    fout << """</script></head><body>"""
    fout << """<div id="plot" class="plotly-graph-div"></div>"""
    fout << """<script>
        Plotly.newPlot("plot", [
            $(join(gA, ',')),
            $(join(gT, ',')),
            $(join(gC, ',')),
            $(join(gG, ',')),
            $(join(gL, ',')),
            $(join(gM, ','))
        ], $layout)
    </script>"""
    fout << """</body></html>"""
end