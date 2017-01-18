#!/usr/bin/env julia
if length(ARGS) <= 1
    print("""

        Usage: run_with_limit limit [-"prefix"] task1 taks2 task3 ... > output

        Example1: setsid run_with_limit 64 *.sh > ~/output 2> ~/log

        Example2: setsid ding run_with_limit 64 -"julia extract_vcf" **/*_MrBam.txt > ~/output 2> ~/log

        use `grep -v "code 0" log` to check if all task succeeded

        """)
    exit(0)
end

cond, nlimit, nruning, tasks = Condition(), parse(Int, ARGS[1]), Ref(0), ARGS[2:end]

0 < nlimit <= length(tasks) || (nlimit = length(tasks))

while !isempty(tasks) || nruning[] != 0
    if isempty(tasks) || nruning[] == nlimit
        wait(cond)
        continue
    end

    let
        task = shift!(tasks)
        nruning[] += 1
        starttime = time()
        @schedule try
            p = spawn(`bash $task`, (DevNull, STDOUT, STDOUT))
            close(p.in)
            wait(p.exitnotify)
            duration = round(Int, time() - starttime)
            println(STDERR, "$task exit with code $(p.exitcode) in $duration seconds")
        catch e
            println(STDERR, "$task failed with Exception: $e")
        finally
            nruning[] -= 1
            notify(cond)
        end
    end
end

println(STDERR, "Finished~")