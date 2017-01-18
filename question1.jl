x = try uppercase(ARGS[1]) catch "CCCATCTAGCCCCGCTAGATCCC" end
complement(a,b) = (a, b) in (('A','T'), ('T','A'), ('C','G'), ('G','C'))
foo(s) = s == length(x) ? (0, s, s) : max(bar(s, length(x)), foo(s+1))
bar(s, e) = e <= s ? (0, s, e) : max((baz(s,e), s, e), bar(s, e-1))
baz(s, e) = e <= s || !complement(x[s], x[e]) ? 0 : baz(s+1, e-1) + 1
l, s, e = foo(1)
println(x, '\n', " "^(s-1), "*"^l, " "^(e-s-2l+1), "*"^l)