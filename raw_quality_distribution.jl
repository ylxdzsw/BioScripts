collect_Q = s"""
#include <stdio.h>

int main(int argc, char** argv){
    int c, l, t = argv[1][0];
    while(getchar() != EOF){ // read a byte, but it will be ignored anyway
        while(getchar() != '\n');
        while(getchar() != '\n');
        while(getchar() != '\n');

        for(l = c = -1; c != '\n'; c = getchar()){
            l += c <= t;
        }

        printf("%d\n", l);
    }
    return 0;
}
"""

open("/tmp/collect_Q.c", "w") do f
    write(f, collect_Q)
end

run(`gcc /tmp/collect_Q.c -o /tmp/collect_Q`)

qual = parse(Int, ARGS[1]) + 33

@everywhere len = 152

@everywhere samples = let
    R1 = Set{String}()
    R2 = Set{String}()

    for i in readdir(".")
        endswith(i, "_R1.fq.gz") ?  push!(R1, i[1:end-9]) :
        endswith(i, "_R2.fq.gz") && push!(R2, i[1:end-9])
    end

    collect(R1 ∩ R2)
end

@sync for i in samples, r in (1, 2)
    @async run(`bash -c "gzip -cd $(i)_R$(r).fq.gz | /tmp/collect_Q $qual > /tmp/$(i)_R$(r).collect"`)
end

@everywhere function reverse_acc_percentage(x)
    result, acc = similar(x), 0

    for i in len:-1:1
        acc += x[i]
        result[i] = acc
    end

    result ./ sum(x)
end

plots = @parallel (*) for i in samples
    R1, R2, either = zeros(Int64, len), zeros(Int64, len), zeros(Int64, len)

    for (l1, l2) in zip(eachline("/tmp/$(i)_R1.collect"), eachline("/tmp/$(i)_R2.collect"))
        l1, l2 = parse(Int, l1), parse(Int, l2)
        R1[l1+1] += 1
        R2[l2+1] += 1
        either[max(l1, l2)+1] += 1
    end

    """
        <div id="$i" class="plotly-graph-div"></div>
        <script>
            Plotly.newPlot("$i", [{
                x: window.data,
                y: [$(join(R1, ','))],
                name: "r1",
                type: "bar"
            }, {
                x: window.data,
                y: [$(join(R2, ','))],
                name: "r2",
                type: "bar"
            }, {
                x: window.data,
                y: [$(join(R1.+R2, ','))],
                name: "both",
                type: "bar"
            }, {
                x: window.data,
                y: [$(join(either, ','))],
                name: "max",
                type: "bar"
            }, {
                x: window.data,
                y: [$(join(reverse_acc_percentage(R1), ','))],
                name: "r1_acc",
                yaxis: "y2",
                type: "scatter"
            }, {
                x: window.data,
                y: [$(join(reverse_acc_percentage(R2), ','))],
                name: "r2_acc",
                yaxis: "y2",
                type: "scatter"
            }, {
                x: window.data,
                y: [$(join(reverse_acc_percentage(R1.+R2), ','))],
                name: "both_acc",
                yaxis: "y2",
                type: "scatter"
            }, {
                x: window.data,
                y: [$(join(reverse_acc_percentage(either), ','))],
                name: "max_acc",
                yaxis: "y2",
                type: "scatter"
            }], {
                title: "$i - Q$qual",
                height: 800,
                margin: {r:180, l:120, b:80, t:80},
                xaxis: {title: 'number of bases whose quality <= $qual'},
                yaxis: {title: 'number of reads'},
                yaxis2: {
                    title: 'percentage_acc',
                    overlaying: 'y',
                    side: 'right'
                }
            })
        </script>
    """
end

open("/tmp/Q$qual.html", "w") do f
    print(f, """
        <html>
            <head>
                <script src="https://cdn.plot.ly/plotly-1.2.0.min.js" />
                <script> window.data = [$(join(0:len-1, ','))] </script>
            </head>
            <body>
                $plots
            </body>
        </html>
    """)
end