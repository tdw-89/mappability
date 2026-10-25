#!/usr/bin/env -S julia --startup-file=no -t 1,0
#
# check_score.jl <simulate.toml> <score.toml>
#
# Every simulated pair must have been seen by `unimap score` exactly once, with
# both mates. score itself already refuses duplicated or separated pairs, but it
# counts a pair that is absent from the BAM as a fragment that failed to map --
# so a lost shard would come out as low mappability rather than as an error.
# This makes it one.
#
# Standard library only, so it runs with or without the UniMap project.

using TOML

length(ARGS) == 2 || (println(stderr, "usage: check_score.jl <simulate.toml> <score.toml>"); exit(2))
sim   = TOML.parsefile(ARGS[1])
score = TOML.parsefile(ARGS[2])

expected = sim["output"]["n_pairs_per_L"]
stats    = score["bam_stats"]

ok = true
for L in sort(collect(keys(expected)); by=x -> parse(Int, x))
    want = expected[L]
    if !haskey(stats, L)
        println(stderr, "L=$L: no BAM statistics in score.toml")
        global ok = false
        continue
    end
    got, incomplete = stats[L]["pairs"], stats[L]["incomplete_pairs"]
    if got != want || incomplete != 0
        println(stderr, "L=$L: simulated $want pairs, scored $got ($incomplete with a mate missing)")
        global ok = false
    else
        println(stderr, "L=$L: all $want simulated pairs scored")
    end
end

ok || (println(stderr, "error: the BAMs do not hold every simulated pair exactly once"); exit(1))
