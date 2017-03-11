using OhMyJulia

open("fuck.html", "w") do fout
    fout << """<html><head><script>"""
    fout << open(read, s"C:\Users\ylxdzsw\.julia\v0.5\PlotlyJS\deps\plotly-latest.min.js")
    fout << """;window.data=""" << open(read, "weight15.json")
    fout << """</script></head><body>"""

    # for i in 1:8
    #     fout << """<div id="filter-$i" class="plotly-graph-div"></div>"""
    #     fout << """<script>
    #         Plotly.newPlot("filter-$i", [
    #             $(join(map(1:8) do c
    #                 """{
    #                     x: data[$i-1][$c-1].x,
    #                     y: data[$i-1][$c-1].y,
    #                     z: Array(25).fill($c),
    #                     mode: "markers",
    #                     type: "scatter3d",
    #                     marker: {
    #                         size: data[$i-1][$c-1].z.map(x=>100*(x+0.1)),
    #                         opacity: 0.4
    #                     },
    #                     name: "$(("A", "T", "C", "G", "I", "D", "SNP", "Reverse")[c])"
    #                 }"""
    #             end, ','))
    #         ], {
    #             title: "Filter-$i"
    #         })
    #     </script>"""
    # end

    for i in 1:64
        fout << """<div id="filter-$i" class="plotly-graph-div"></div>"""
        fout << """<script>!function(){
            var mdzz = [[], [], [], [], [], [], [], []]
            for(var x = 0; x<25; x++){
                var best = 0
                for(var c = 0; c<8; c++) {
                    if (Math.abs(data[$i-1][c][x]) > Math.abs(data[$i-1][best][x])) {
                        best = c
                    }
                }
                mdzz[best].push(x)
            }
            Plotly.newPlot("filter-$i", [
                $(join(map(1:8) do c
                    """{
                        x: mdzz[$c-1].map(x=>Math.floor(x/5)),
                        y: mdzz[$c-1].map(x=>x%5),
                        mode: "markers",
                        marker: {
                            symbol: mdzz[$c-1].map(x=>data[$i-1][$c-1][x]>0?'circle':'diamond'),
                            size: mdzz[$c-1].map(x=>280*Math.abs(data[$i-1][$c-1][x]))
                        },
                        name: "$(("A", "T", "C", "G", "I", "D", "SNP", "Reverse")[c])"
                    }"""
                end, ','))
            ], {
                title: "Filter-$i",
                width: 600,
                height: 550
            })
        }()</script>"""
    end

    fout << """</body></html>"""
end
