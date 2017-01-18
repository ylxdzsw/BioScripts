const N  = big(250_000)
const pm = big(0.1)
const r  = big(0.2)
const k  = big(50)

binomial_cdf(n, p, m) = begin
    mapreduce(+, 0:m) do i
        (i==0?1:(*(n-i+1:n...)/*(1:i...)))*p^i*(1-p)^(n-i)
    end
end

for a in big(5:10), b in big(5:10)
    p = (1 / 2^b) * (1 - binomial_cdf(a, pm, ceil(BigInt, a*r)))
    p = binomial_cdf(N, p, k)
    @printf("%d\t%d\t%.12f\n", a, b, p)
end
