using Distributed
@everywhere using Distributions
using Test
using BenchmarkTools
using Dates

@everywhere struct Hack end
function fixRC()
    for p in workers()
        @fetchfrom p Hack()
    end
end
fixRC()
println("We have $(nprocs()-1) workers")
exit(0)
