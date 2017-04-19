using OhMyJulia
using JSON: json

const labels, data = let
    x = split.(readlines(car(ARGS)))

    head, body = car(x), []

    group = groupby(car, (x, y) -> begin
        car(x) .+ map(x->parse(f64, x), cdr(y)), cadr(x)+1
    end, ()->(fill(0, length(head)-1), 0), cdr(x))

    group = collect((k, car(v) ./ cadr(v)) for (k, v) in group)

    for i in 1:length(head)-1
        k = map(car, group)
        v = map(x->cadr(x)[i], group) |> normalize
        for (k, v) in zip(k, v)
            push!(body, [k, v, head[i+1]])
        end
    end

    head, body
end

STDOUT << """<html><head><meta charset="utf-8"><script>"""
STDOUT << open(read, "echarts.js")
STDOUT << """</script></head><body>"""

STDOUT << """<div id="plot" style="width:800px;height:600px;"></div>"""
STDOUT << """<script>
    echarts.init(document.getElementById('plot')).setOption({
        tooltip: {
            trigger: 'axis',
            axisPointer: {
                type: 'line',
                lineStyle: {
                    color: 'rgba(0,0,0,0.2)',
                    width: 1,
                    type: 'solid'
                }
            }
        },

        legend: {
            data: $(json(cdr(labels)))
        },

        singleAxis: {},

        series: [{
            type: 'themeRiver',
            itemStyle: {
                emphasis: {
                    shadowBlur: 20,
                    shadowColor: 'rgba(0, 0, 0, 0.8)'
                }
            },
            label: {
                normal: {
                    show: false
                }
            },
            data: $(json(data))
        }]
    })
</script>"""

STDOUT << """</body></html>"""
