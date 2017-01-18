using OhMyJulia

@everywhere begin

#= user configuration start =#

const pA2T, pA2C, pA2G = .000_028, .000_038, .000_014
const pT2A, pT2C, pT2G = .000_027, .000_015, .000_038
const pC2A, pC2T, pC2G = .000_051, .000_011, .000_030
const pG2A, pG2T, pG2C = .000_010, .000_052, .000_030

const N = 200_000
const FREQSTEP = .0:.001:.2
const MAXDEPTH = 1000
const MAX_RAW_FREQ = .24 # must larger than max(FREQSTEP) but less than .5
const CONFIDENCE_LEVELS = (.90, .98)
const COLLISION_RATE = 4.17e-5 # caculated from an inset length distribution that centered at 151. the smaller the center, the higher the possibility of collision

#= user configuration end =#

const THRESHOLDS = map(x->floor(Int, N*x), CONFIDENCE_LEVELS)

function trans(origin::Char)
    c = rand()
    origin == 'A' ? c > pA2T+pA2C+pA2G ? 'A' : c > pA2T+pA2C ? 'G' : c > pA2T ? 'C' : 'T' :
    origin == 'T' ? c > pT2A+pT2C+pT2G ? 'T' : c > pT2A+pT2C ? 'G' : c > pT2A ? 'C' : 'A' :
    origin == 'C' ? c > pC2A+pC2T+pC2G ? 'C' : c > pC2A+pC2T ? 'G' : c > pC2A ? 'T' : 'A' :
    origin == 'G' ? c > pG2A+pG2T+pG2C ? 'G' : c > pG2A+pG2T ? 'C' : c > pG2A ? 'T' : 'A' :
    error("unknown base")
end

max_alt(depth) = ceil(Int, MAX_RAW_FREQ * depth * (1 + (depth < 30)))

collision(n) = rand() < COLLISION_RATE * n

end # everywhere

nprocs() == 1 && println(STDERR, "single process would be slow; use -p to start more processes")

println(join(map(x->100x, CONFIDENCE_LEVELS), ' '))

P = SharedArray(Float64, (length(THRESHOLDS), max_alt(MAXDEPTH), MAXDEPTH))

for ref = "ATCG", alt = "ATCG" @when ref != alt
    @sync @parallel for d in 1:MAXDEPTH
        P[:, 1:max_alt(d), d] = .0
        for r in FREQSTEP
            n = fill(0, d+1) # n[i+1] is the number of having i reads supporting alt
            for i in 1:N
                s, j = 0, 0
                while j < d
                    x = rand() < r ? alt : ref
                    if trans(x) == alt
                        collision(s) || (j += 1; s += 1)
                    else
                        collision(j-s) || (j += 1)
                    end
                end
                n[s+1] += 1
            end
            s, level = n[1], 0
            for i in 1:max_alt(d)
                s += n[i+1]
                while level != length(THRESHOLDS) && s > THRESHOLDS[level+1]
                    level += 1
                end
                level == 0 && continue
                P[1:level, i, d] = r
                if level == length(THRESHOLDS)
                    P[1:level, i+1:max_alt(d), d] = r
                    break
                end
            end
        end
    end

    println(ref, "->", alt, ':')

    for d in 1:MAXDEPTH
        print(max_alt(d), ':')
        for i in 1:max_alt(d)
            for x in 1:length(THRESHOLDS)
                @printf(" %.1f%%", 100P[x, i, d])
            end
        end
        println()
    end
end
