using OhMyJulia
using Fire

const CLIP, DIAG, UP, LEFT = 0x01, 0x02, 0x04, 0x08
const MAT,  SNP, INS,  DEL = 0x10, 0x20, 0x40, 0x80

function dp_fill!(q::AbstractBytes, t::AbstractBytes, F::Matrix{Int}, P::Matrix{Byte}, mm=12, go=18, ge=3)::Void
    F[1, :] = t .== q[1]
    P[1, :] = CLIP + DIAG * (t .== q[1])

    F[:, 1] = q .== t[1]
    P[:, 1] = CLIP + DIAG * (q .== t[1])

    # core dp loop: fill each cell with best score and path
    for i in 2:length(q), j in 2:length(t)
        score, case = findmax([
            0,
            F[i-1, j-1] - (q[i]      == t[j] ? -1 : mm),
            F[i-1, j]   - (P[i-1, j] == UP   ? ge : go),
            F[i, j-1]   - (P[i, j-1] == LEFT ? ge : go),
        ])

        F[i, j] = score
        P[i, j] = 0x01 << (case - 1)
    end

    # find optimum path for each start position
    for i in length(q):-1:1, j in length(t):-1:1
        if P[i, j] & CLIP != 0
            continue
        elseif P[i, j] & DIAG != 0
            if F[i, j] > F[i-1, j-1]
                P[i-1, j-1] &= 0x0f
                P[i-1, j-1] |= q[i] == t[j] ? MAT : SNP
                F[i-1, j-1] = F[i, j]
            end
        elseif P[i, j] & UP != 0
            if F[i, j] > F[i-1, j]
                P[i-1, j] &= 0x0f
                P[i-1, j] |= DEL
                F[i-1, j] = F[i, j]
            end
        elseif P[i, j] & LEFT != 0
            if F[i, j] > F[i, j-1]
                P[i, j-1] &= 0x0f
                P[i, j-1] |= INS
                F[i, j-1] = F[i, j]
            end
        else
            error("BUG")
        end
    end
end

function get_variants!(i::Int, j::Int, P::Matrix{Byte}, buf)
    l = 0

    if P[i, j] & DIAG != 0
        l += 1
        buf[1, l] = true
        buf[2, l] = false
        buf[3, l] = true
    end

    while P[i, j] >> 4 != 0
        println(bits(P[i,j]))
        l += 1
        if P[i, j] & (MAT | SNP) != 0
            i += 1
            j += 1

            if P[i-1, j-1] & SNP != 0
                buf[:, l] = true
            else
                buf[1, l] = true
                buf[2, l] = false
                buf[3, l] = true
            end
        elseif P[i, j] & INS != 0
            buf[1, l] = false
            buf[2, l] = true
            buf[3, l] = true
            j += 1
        elseif P[i, j] & DEL != 0
            buf[1, l] = true
            buf[2, l] = true
            buf[3, l] = false
            i += 1
        end
    end
    l
end

function write_alignment(io::IO, start_pos, q, t, l, buf, code_transform=identity)
    i, j = start_pos
    for k in 1:l
        if buf[1, k]
            io << code_transform(q[i])
            i += 1
        else
            io << '-'
        end
    end
    println(io)
    for k in 1:l
        io << (buf[2, k] ? ' ' : '|')
    end
    println(io)
    for k in 1:l
        if buf[3, k]
            io << code_transform(t[j])
            j += 1
        else
            io << '-'
        end
    end
    println(io)
end

"""
the return value is [start-position-in-q, start-position-in-t, str]
"""
function dp_report(q::AbstractBytes, t::AbstractBytes, F::Matrix{Int}, P::Matrix{Byte},
                   th::Int=64, code_transform::Function=identity)
    result, buf = [], Matrix{Bool}(3, length(q) + length(t))
    for i in 1:length(q)-th, j in 1:length(t)-th # skip those shorter than th
        if F[i, j] >= th && P[i, j] & CLIP != 0
            start_pos = P[i, j] & DIAG != 0 ? (i, j) : (i+1, j+1)
            l = get_variants!(i, j, P, buf)
            str = sprint(write_alignment, start_pos, q, t, l, buf, code_transform)
            push!(result, (start_pos..., str))
        end
    end
    result
end

"Align two sequences"
@main function main(a::String, b::String; mismatch::Int=12, gap_open::Int=18, gap_extend::Int=3, threshold::Int=64)
    a, b = a.data, b.data
    l1, l2 = length(a), length(b)
    F = Matrix{Int}(l1, l2)
    P = Matrix{Byte}(l1, l2)
    dp_fill!(a, b, F, P, mismatch, gap_open, gap_extend)
    v = dp_report(a, b, F, P, threshold)
    foreach(v->println(v[3]), v)
end
