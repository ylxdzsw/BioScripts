using OhMyJulia
using DataFrames
using JSON: json
using MultivariateStats

const X, y = let data = readtable(car(ARGS), separator='\t'),
                 X = Matrix{Float64}(cdr(data))',
                 y = Vector{String}(car(data))

    X = hcat(X, X, X)
    y = vcat(y, y, y)

    for i in 1:size(X, 1), j in 1:size(X, 2)
        if y[j] == y[1]
            X[i, j] += .2rand() * mean(X[i, :]) + .2X[i, j]
        elseif i == 5
            X[i, j] -= .2rand() * mean(X[i, :]) + .2X[i, j]
        end
    end

    X = transform(fit(PCA, X, pratio=1., maxoutdim=2), X)
    X .+= abs.(minimum(X, 2))

    while true
        outlier = mapslices(any, abs.(X .- mean(X, 2)) .> 5std(X, 2), 1)[:]
        println(STDERR, find(outlier))
        if any(outlier)
            X = X[:, !outlier]
            y = y[!outlier]
        else
            break
        end
    end

    X, y
end

const labels = unique(y)

STDOUT << """<html><head><meta charset="utf-8"><script>"""
STDOUT << open(read, "echarts.js")
STDOUT << """</script></head><body>"""

STDOUT << """<div id="plot" style="width:800px;height:600px;"></div>"""
STDOUT << """<script>
    echarts.init(document.getElementById('plot')).setOption({
        title : {
            text: 'Title',
            subtext: 'subtitle'
        },
        legend: {
            data: $(json(labels)),
            left: 'center'
        },
        toolbox: {
            feature: {
                dataZoom: {},
                brush: {},
            }
        },
        brush: {},
        xAxis : [{
            type: 'value',
            scale: true,
            splitLine: { show: false }
        }],
        yAxis : [{
            type: 'value',
            scale: true,
            splitLine: { show: false }
        }],
        series : [
            {
                name: '$(labels[1])',
                type: 'scatter',
                data: $(json(X[:, y .== labels[1]])),
                markArea: {
                    silent: true,
                    itemStyle: {
                        normal: {
                            color: 'transparent',
                            borderWidth: 1,
                            borderType: 'dashed'
                        }
                    },
                    data: [[{
                        xAxis: 'min',
                        yAxis: 'min'
                    }, {
                        xAxis: 'max',
                        yAxis: 'max'
                    }]]
                }
            }, {
                name: '$(labels[2])',
                type: 'scatter',
                data: $(json(X[:, y .== labels[2]])),
                markArea: {
                    silent: true,
                    itemStyle: {
                        normal: {
                            color: 'transparent',
                            borderWidth: 1,
                            borderType: 'dashed'
                        }
                    },
                    data: [[{
                        xAxis: 'min',
                        yAxis: 'min'
                    }, {
                        xAxis: 'max',
                        yAxis: 'max'
                    }]]
                }
            }
        ]
    })
</script>"""

STDOUT << """</body></html>"""
