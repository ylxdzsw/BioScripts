const DIAG, UP, LEFT, STOP = 0x01, 0x02, 0x03, 0x00

function dp_fill!(q::Bytes, t::Bytes, F::Matrix{Byte}, P::Matrix{Byte},
                  score_metric=(0x30, 0x50, 0x10), limit=0xa0)::Void
    SNP, GO, GE = score_metric

    F[1, :] = t .!= q[1]
    P[1, :] = t .== q[1]

    F[1:32, 1] = (1:32) .- (q[1:32] .== t[1])
    P[1:32, 1] = q[1:32] .== t[1]

    F[33:end, 1] = limit
    P[33:end, 1] = CLIP

    for i in 2:16, j in 2:65536

    end

    for i in 17:112, j in 2:65536
        score, case = findmin([
            0,
            F[i-1, j-1] + (q[i] == t[j] ? 0x00 : SNP),
            F[i-1, j] + (P[i-1, j] == UP   ? GE : GO),
            F[i, j-1] + (P[i, j-1] == LEFT ? GE : GO),
        ])

        if i <= 16 && score > i
            score, case = Byte(i), Int(CLIP)
        end

        F[i,j]   = min(score, limit)
        Ptr[i,j] = case - 1
    end

    for i in 113:128, j in 2:65536

    end

    nothing
end

"""
return: pos in q, variant(0: snp, positive: insertion, negative: deletion)
"""
function traceback(q::Bytes, t::Bytes, F::Matrix{Byte}, P::Matrix{Byte}, limit::Int=0xa0)::Vector{Int, Int}
    opt = findmax(F)

    a1, a2 = let
        a1, a2 = IOBuffer(), IOBuffer()

        function traceback(i, j)
            ptr = Ptr[i, j]

            if ptr == DIAG
                traceback(i-1, j-1)
                a1 << seq1[i-1]
                a2 << seq2[j-1]
            elseif ptr == UP
                traceback(i-1, j)
                a1 << seq1[i-1]
                a2 << '_'
            elseif ptr == LEFT
                traceback(i, j-1)
                a1 << '_'
                a2 << seq2[j-1]
            end

            nothing
        end

        traceback(ind2sub(size(F), cadr(opt))...)

        takebuf_array(a1), takebuf_array(a2)
    end

    a1, a2, car(opt)
end

"""
only output variants within q
for indel, the other should be "-"
"""
function extract_alignments()

end
