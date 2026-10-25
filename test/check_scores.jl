#!/usr/bin/env julia
#
# julia --project=<UniMap.jl> test/check_scores.jl equal  <scoresA> <scoresB>
# julia --project=<UniMap.jl> test/check_scores.jl expect <scores>
#
# <scores> is a `unimap score` output directory with mappability.zarr unpacked.
#
# equal   every Zarr array and every Parquet table is identical between the two.
#         Tables are compared as row sets: bleed/site rows come out of Dicts,
#         whose order is not part of the result.
# expect  the fixture's planted features come out as they must (see
#         test/make_fixture.py).
#
# Prints one line per check; exits 1 if any failed.

using DataFrames, Parquet, Zarr

const FAILED = Ref(0)
check(ok::Bool, what) = (println(ok ? "  ok   " : "  FAIL ", what); ok || (FAILED[] += 1); ok)

table(dir, name) = DataFrame(read_parquet(joinpath(dir, "$name.parquet")); copycols=true)

function equal(a, b)
    za = zopen(joinpath(a, "mappability.zarr"); consolidated=true)
    zb = zopen(joinpath(b, "mappability.zarr"); consolidated=true)
    check(sort(collect(keys(za.arrays))) == sort(collect(keys(zb.arrays))), "same Zarr arrays")
    for name in sort(collect(keys(za.arrays)))
        haskey(zb.arrays, name) || continue
        x, y = za.arrays[name][:], zb.arrays[name][:]
        check(size(x) == size(y) && x == y, "Zarr $name identical ($(length(x)) values, $(sum(Float64, x)) total)")
    end
    check(za.attrs["loci"] == zb.attrs["loci"], "Zarr locus table identical")
    for name in ("misplacement_summary", "bleed_matrix", "misplaced_sites")
        ta, tb = table(a, name), table(b, name)
        check(names(ta) == names(tb) && nrow(ta) == nrow(tb) &&
              Set(map(Tuple, eachrow(ta))) == Set(map(Tuple, eachrow(tb))),
              "$name identical ($(nrow(ta)) rows)")
    end
end

function expect(dir)
    s = table(dir, "misplacement_summary")
    at(locus, cutoff) = s[(s.locus_id .== locus) .& (s.mapq_cutoff .== cutoff), :]
    frac(locus, cutoff) = (d = at(locus, cutoff); sum(d.n_correct_kept) / sum(d.n_fragments))

    loci = Set(s.locus_id)
    check(issubset(Set(["dupA", "uniqB", "nearC", "edgeD", "revE", "nblkF", "sameG"]), loci),
          "all 7 protein-coding loci scored")
    check(isempty(intersect(loci, Set(["lnc1", "pseudo1", "rna-lnc1", "gene-pseudo1"]))),
          "lncRNA and pseudogene not selected")

    for locus in ("uniqB", "nearC", "edgeD", "revE", "sameG")
        f = frac(locus, 2)
        check(f > 0.99, "$locus (unique sequence): Pr(map back, MAPQ>=2) = $(round(f; digits=4)) > 0.99")
    end
    f0, f2 = frac("dupA", 0), frac("dupA", 2)
    # Not ~0: fragments overhanging the ends of the copied segment carry unique
    # sequence and place correctly.
    check(f2 < 0.2, "dupA (paralog on chr2): Pr(map back, MAPQ>=2) = $(round(f2; digits=4)) < 0.2")
    check(f0 > 0.3 && f0 < 0.7, "dupA at MAPQ>=0 is a coin flip between the copies: $(round(f0; digits=4))")

    b = table(dir, "bleed_matrix")
    check(any((b.origin_locus .== "dupA") .& (b.mapq_cutoff .== 0) .& (b.n_fragments .> 0)),
          "dupA fragments land elsewhere (bleed_matrix)")

    nskip = sum(at("nblkF", 0).n_skipped_N)
    check(nskip > 0, "nblkF: $nskip fragments skipped for the N block")
    check(sum(at("uniqB", 0).n_skipped_N) == 0, "uniqB: no N-skipped fragments")
    check(sum(at("uniqB", 0).n_missing) == 0 && sum(at("dupA", 0).n_missing) == 0,
          "no simulated fragment missing from the BAMs")

    z = zopen(joinpath(dir, "mappability.zarr"); consolidated=true)
    check(sort(Int.(z.arrays["L"][:])) == [30, 50], "Zarr holds L = 30, 50")
end

mode = isempty(ARGS) ? "" : ARGS[1]
if mode == "equal" && length(ARGS) == 3
    equal(ARGS[2], ARGS[3])
elseif mode == "expect" && length(ARGS) == 2
    expect(ARGS[2])
else
    println(stderr, "usage: check_scores.jl equal <scoresA> <scoresB> | expect <scores>")
    exit(2)
end
FAILED[] == 0 || exit(1)
