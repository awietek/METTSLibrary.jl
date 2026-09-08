using METTSLibrary
using Test
using HDF5
using TOML
using Random

const ML = METTSLibrary

const TJ_BONDS = [("t", "HOP", :nn), ("J", "HB", :nn), ("t_prime", "HOP", :nnn)]

function synthetic_tj(; Lx=4, Ly=2, nsamples=10, nup=3, ndn=2, beta=2.0, seed=1)
    rng = MersenneTwister(seed)
    N = Lx * Ly
    states = zeros(UInt8, N, nsamples)
    for j in 1:nsamples
        perm = randperm(rng, N)
        states[perm[1:nup], j] .= 1          # Up
        states[perm[nup+1:nup+ndn], j] .= 2  # Dn
    end
    Ensemble(
        model="tJ", site_type="tJ",
        lattice=square_lattice_toml(Lx, Ly; yperiodic=true, bonds=TJ_BONDS),
        lattice_name="square.L$(Lx).W$(Ly).cyl",
        beta=beta,
        parameters=Dict("t" => 3.0, "J" => 0.4, "t_prime" => -0.3),
        sector=Dict("nup" => nup, "ndn" => ndn),
        algorithm=Dict{String,Any}("code" => "test", "tau" => 0.1, "maxdim" => 100),
        provenance=Dict{String,Any}("creator" => "runtests"),
        collapse_bases=["Z"],
        states=states,
        observables=Dict{String,Array{Float64}}("energy" => randn(rng, nsamples),
                                                "entropy" => rand(rng, N - 1, nsamples)),
    )
end

# copy of an ensemble with some fields replaced
with(e::Ensemble; kw...) = Ensemble(; (k => getfield(e, k) for k in fieldnames(Ensemble))..., kw...)

const MINI_TOML = """
# a comment
Coordinates = [
  [0.0, 0.0],
  [1.0, 0.0],
  [0.0, 1.0],
  [1.0, 1.0]
]

Interactions = [
  ['J', 'HB', 0, 1],
  ['J', 'HB', 2, 3],
  ['Jd', 'HB', 0, 3]
]

Symmetries = [
  [0, 1, 2, 3]
]
"""

const MINI_LAT_CPP = """
[Coordinates]
0 0
1 0
0 1
1 1
[Interactions]
HB J 0 1
HB J 2 3
HB Jd 0 3
"""

const MINI_LAT_XDIAG = """
# header comment
[Dimension]=2
[Sites]=4
0.0000000000000000 0.0000000000000000
1.0000000000000000 0.0000000000000000
0.0000000000000000 1.0000000000000000
1.0000000000000000 1.0000000000000000
[Interactions]=3
HB J 0 1
HB J 2 3
HB Jd 0 3
[SymmetryOps]=1
[S0] 0 1 2 3
"""

@testset "METTSLibrary" begin

@testset "lattice files" begin
    l = parse_lattice(MINI_TOML)
    @test nsites(l) == 4
    @test l.coordinates[:, 4] == [1.0, 1.0]
    @test length(l.interactions) == 3
    @test l.interactions[3] == ("Jd", "HB", [0, 3])
    @test lattice_couplings(l) == ["J", "Jd"]

    for txt in (MINI_LAT_CPP, MINI_LAT_XDIAG)
        l2 = parse_lattice(ML.lat_to_toml(txt))
        @test l2.coordinates == l.coordinates
        @test l2.interactions == l.interactions
    end
    @test startswith(ML.lat_to_toml(MINI_LAT_XDIAG), "# header comment")

    mktempdir() do dir
        ft = joinpath(dir, "a.toml"); write(ft, MINI_TOML)
        fl = joinpath(dir, "a.lat");  write(fl, MINI_LAT_CPP)
        @test read_lattice(ft) == MINI_TOML
        @test parse_lattice(read_lattice(fl)).interactions == l.interactions
        @test_throws ErrorException read_lattice(joinpath(dir, "missing.toml"))
    end

    @test_throws ErrorException parse_lattice("Interactions = []")
    @test_throws ErrorException parse_lattice("Coordinates = [[0,0],[1]]")
    @test_throws ErrorException parse_lattice("Coordinates = [[0,0]]\nInteractions = [['J']]")

    # generated square lattice, ITensorMPS ordering: n = x*Ly + y
    l3 = parse_lattice(square_lattice_toml(4, 2; yperiodic=true, bonds=[("t", "HOP", :nn)]))
    @test nsites(l3) == 8
    @test l3.coordinates[:, 1] == [0.0, 0.0]
    @test l3.coordinates[:, 2] == [0.0, 1.0]
    @test l3.coordinates[:, 3] == [1.0, 0.0]
    @test length(l3.interactions) == 3 * 2 + 4          # 6 x-bonds, 4 y-bonds (deduplicated wrap)
    @test all(all(0 .<= s .< 8) for (_, _, s) in l3.interactions)
    l4 = parse_lattice(square_lattice_toml(3, 3; yperiodic=false, bonds=[("t", "HOP", :nn), ("tp", "HOP", :nnn)]))
    @test count(i -> i[1] == "t", l4.interactions) == 12
    @test count(i -> i[1] == "tp", l4.interactions) == 8
end

@testset "schema and validation" begin
    e = synthetic_tj()
    @test nsites(e) == 8
    @test nsamples(e) == 10
    @test validate(e) === e
    @test lattice(e) isa Lattice
    @test nsites(lattice(e)) == 8
    @test occursin("square.L4.W2.cyl (8 sites)", sprint(show, e))

    @test_throws ArgumentError validate(with(e; sector=Dict("nup" => 1)))
    st = copy(e.states); st[1, 1] = 7
    @test_throws ArgumentError validate(with(e; states=st))
    @test_throws ArgumentError validate(with(e; observables=Dict{String,Array{Float64}}("energy" => zeros(3))))
    # X-basis samples are not sector checked
    @test validate(with(e; collapse_bases=["X"], sector=Dict("nup" => 1))) isa Ensemble
    # lattice consistency
    @test_throws ArgumentError validate(with(e; parameters=Dict("t" => 3.0, "J" => 0.4)))   # t_prime missing
    @test_throws ArgumentError validate(with(e; lattice=MINI_TOML))                           # 4 sites vs 8
    @test_throws ArgumentError validate(with(e; lattice="not toml ["))
    @test_throws ArgumentError validate(with(e; lattice_name=""))
    @test_throws ArgumentError validate(with(e; lattice_name="a/b"))
    # extra parameters not in the lattice are fine
    @test validate(with(e; parameters=merge(e.parameters, Dict("U" => 1.0)))) isa Ensemble
end

@testset "paths" begin
    e = synthetic_tj()
    rel = ML.relpath_for(e; tag="chain01")
    @test rel == joinpath("tJ", "square.L4.W2.cyl", "J=0.4_t=3.0_t_prime=-0.3", "ndn=2_nup=3",
                          "beta=2.0", "chain01.h5")
    @test_throws ArgumentError ML.relpath_for(e; tag="")
    @test_throws ArgumentError ML.relpath_for(e; tag="a/b")
    @test ML.unflatten_path(ML.flatten_path(rel)) == rel
    @test ML.lattice_path("/lib/tJ/sq/p/b.h5", "sq") == "/lib/tJ/sq/sq.toml"
end

@testset "write, read, lattice sharing" begin
    e = synthetic_tj()
    mktempdir() do root
        rel = write_ensemble(root, e; tag="a")
        @test rel == ML.relpath_for(e; tag="a")
        path = joinpath(root, rel)
        @test isfile(path)
        @test_throws ErrorException write_ensemble(root, e; tag="a")     # append-only

        # the lattice lives one directory above the HDF5 file, named after the lattice
        latf = joinpath(root, "tJ", "square.L4.W2.cyl", "square.L4.W2.cyl.toml")
        @test isfile(latf)
        @test read(latf, String) == e.lattice
        @test isempty(filter(endswith(".toml"), readdir(dirname(path))))

        r = read_ensemble(path)
        for k in fieldnames(Ensemble)
            k in (:algorithm, :provenance, :observables) && continue
            @test getfield(r, k) == getfield(e, k)
        end
        @test r.observables["energy"] ≈ e.observables["energy"]
        @test r.observables["entropy"] ≈ e.observables["entropy"]
        @test r.algorithm["code"] == "test"
        @test r.algorithm["maxdim"] == 100
        @test r.provenance["creator"] == "runtests"
        h5open(path, "r") do f
            @test read(f["coordinates"])[:, 2] == [0.0, 1.0]
            # schema 2 dropped /chain: one file is one run, seed lives in /algorithm
            @test !haskey(f, "chain")
            @test read_attribute(f, "schema_version") == 2
        end

        m = ML.read_metadata(path)
        @test m["nsamples"] == 10
        @test m["parameters"]["J"] == 0.4
        @test m["sector"]["nup"] == 3
        @test m["observables"] == ["energy", "entropy"]
        @test m["couplings"] == ["t", "J", "t_prime"]
        @test length(ML.sha256_file(path)) == 64
        @test !ML.is_lfs_pointer(path)

        # other temperatures and parameter sets on the same lattice share the file
        write_ensemble(root, with(e; beta=3.0); tag="a")
        write_ensemble(root, with(e; parameters=merge(e.parameters, Dict("J" => 0.5))); tag="a")
        @test length(filter(endswith(".toml"), readdir(joinpath(root, "tJ", "square.L4.W2.cyl")))) == 1
        # a different lattice under the same name is refused
        other = with(e; lattice=square_lattice_toml(4, 2; yperiodic=false, bonds=TJ_BONDS))
        @test_throws ErrorException write_ensemble(root, other; tag="a")
        # ... and accepted under its own name
        rel4 = write_ensemble(root, with(other; lattice_name="square.L4.W2.open"); tag="a")
        latf4 = joinpath(root, "tJ", "square.L4.W2.open", "square.L4.W2.open.toml")
        @test isfile(latf4)
        @test read_ensemble(joinpath(root, rel4)).lattice_name == "square.L4.W2.open"
        # missing or modified lattice file is detected on read
        mv(latf4, latf4 * ".bak")
        @test_throws ErrorException read_ensemble(joinpath(root, rel4))
        write(latf4, e.lattice)
        @test_throws ErrorException read_ensemble(joinpath(root, rel4))
        @test_throws ErrorException ML.read_metadata(joinpath(root, rel4))
    end
end

@testset "states handoff" begin
    e = synthetic_tj()
    s = ML.state_indices(e, 1)
    @test s isa Vector{Int}
    @test all(1 .<= s .<= 3)
    @test state_labels(e, 1) == e.local_states[s]
    @test count(==("Up"), state_labels(e, 1)) == 3

    ss = initial_states(e, 4; rng=MersenneTwister(0))
    @test length(ss) == 4
    @test all(length(v) == nsites(e) for v in ss)
    @test length(initial_states(e)) == nsamples(e)
    @test initial_states(e)[1] == s
    @test length(initial_states(e; thin=2)) == 5
    @test length(initial_states(e, 3; thin=2)) == 3
    @test_throws ArgumentError initial_states(e, 100)
    @test_throws ArgumentError initial_states(e; basis="X")
end

@testset "index, query, local fetch" begin
    mktempdir() do root
        e1 = synthetic_tj(beta=2.0)
        e2 = synthetic_tj(beta=4.0, seed=2)
        e3 = synthetic_tj(Lx=6, beta=2.0, seed=3)
        for e in (e1, e2, e3)
            write_ensemble(root, e; tag="a")
        end
        mkpath(joinpath(root, ".git")); touch(joinpath(root, ".git", "x.h5"))

        idx = build_index(root)
        @test isfile(joinpath(root, "index.toml"))
        @test length(idx["ensembles"]) == 3
        parsed = TOML.parsefile(joinpath(root, "index.toml"))
        @test length(parsed["ensembles"]) == 3
        first_entry = parsed["ensembles"][1]
        @test first_entry["couplings"] == ["t", "J", "t_prime"]
        @test first_entry["lattice"] == joinpath("tJ", "square.L4.W2.cyl", "square.L4.W2.cyl.toml")
        @test first_entry["lattice_name"] == "square.L4.W2.cyl"
        @test length(first_entry["lattice_sha256"]) == 64
        @test length(unique(x["lattice"] for x in parsed["ensembles"])) == 2

        @test length(ensembles(parsed; model="tJ")) == 3
        @test length(ensembles(parsed; beta=2.0)) == 2
        @test length(ensembles(parsed; nsites=12)) == 1
        @test length(ensembles(parsed; lattice_name="square.L6.W2.cyl")) == 1
        @test length(ensembles(parsed; t=3.0, nup=3)) == 3
        @test length(ensembles(parsed; J=0.5)) == 0
        @test length(ensembles(parsed; code="test")) == 3

        withenv("METTSLIBRARY_PATH" => root, "METTSLIBRARY_CACHE" => joinpath(root, "cache")) do
            @test ML.data_root() == root
            li = load_index()
            @test length(li["ensembles"]) == 3
            entry = only(ensembles(li; nsites=12))
            @test ML.fetch_file(entry["path"]; sha256=entry["sha256"]) == joinpath(root, entry["path"])
            e = load(entry)
            @test nsites(e) == 12
            @test nsamples(e) == 10
            @test e.lattice_name == "square.L6.W2.cyl"
        end

        # fetching into a cache reproduces the tree so the lattice reference works
        mktempdir() do cache
            withenv("METTSLIBRARY_PATH" => nothing, "METTSLIBRARY_CACHE" => cache) do
                clear_remotes!()
                add_remote!("file://" * root)
                li = load_index(source=root)
                entry = only(ensembles(li; nsites=12))
                e = load(entry)
                @test nsamples(e) == 10
                @test isfile(joinpath(cache, entry["lattice"]))
                clear_remotes!()
            end
        end
    end
end

@testset "converter: legacy C++ dump" begin
    mktempdir() do dir
        N, M = 6, 5
        ps = Matrix{Float64}(Float64[0 1 2 0 1 2; 1 1 2 0 0 2; 2 1 0 1 0 2; 0 0 1 2 1 2; 1 2 0 2 1 0]')  # (N, M)
        dump = joinpath(dir, "run.dump.h5")
        h5open(dump, "w") do f
            f["ProductState"] = ps
            f["H"] = reshape(collect(1.0:M), 1, M)
            f["Entropy"] = reshape(fill(0.5, M), 1, M)
            f["Sz"] = zeros(N, M)               # per-site, not per-step scalar: ignored
        end
        lat = joinpath(dir, "chain.lat")
        write(lat, "[Coordinates]\n" * join(string.(0:N-1), "\n") *
                   "\n[Interactions]\n" * join(["HOP T $i $(i+1)" for i in 0:N-2], "\n") * "\n")

        kw = (model="tJ", site_type="tJ", temperature=0.5, basis="Z",
              parameters=Dict("T" => 1.0, "J" => 0.3), sector=Dict("nup" => 2, "ndn" => 2))
        e = from_legacy_cpp(dump; lattice=lat, kw...)
        @test e.lattice_name == "chain"                     # from the file name
        @test nsites(e) == N
        @test nsamples(e) == M
        @test e.beta == 2.0
        @test e.states[:, 1] == UInt8[0, 1, 2, 0, 1, 2]
        @test e.observables["energy"] == collect(1.0:M)
        @test e.observables["entropy"] == fill(0.5, M)
        @test !haskey(e.observables, "sz")
        @test size(lattice(e).coordinates) == (1, N)
        @test lattice_couplings(lattice(e)) == ["T"]
        @test e.provenance["source_format"] == "legacy_cpp_dump_h5"
        @test e.provenance["source_lattice"] == abspath(lat)
        @test e.collapse_bases == ["Z"]

        # coupling in the lattice without a parameter value is rejected
        @test_throws ArgumentError from_legacy_cpp(dump; lattice=lat, kw..., parameters=Dict("J" => 0.3))
        # lattice may also be TOML text, then a name is required
        @test_throws ArgumentError from_legacy_cpp(dump; lattice=read_lattice(lat), kw...)
        e2 = from_legacy_cpp(dump; lattice=read_lattice(lat), lattice_name="chain6", kw...)
        @test e2.states == e.states
        @test e2.lattice_name == "chain6"
        @test !haskey(e2.provenance, "source_lattice")

        # UpDn is not allowed for tJ, fine for Electron
        h5open(dump, "r+") do f
            d = f["ProductState"]; v = read(d); v[1, 1] = 3.0; write(d, v)
        end
        @test_throws ErrorException from_legacy_cpp(dump; lattice=lat, kw...)
        e3 = from_legacy_cpp(dump; lattice=lat, model="Hubbard", site_type="Electron", beta=1.0,
                             parameters=Dict("T" => 1.0))
        @test e3.states[1, 1] == 3
        @test startswith(ML.relpath_for(e3; tag="x"), joinpath("Hubbard", "chain"))
    end
end

@testset "converter: METTS.jl samples.txt" begin
    mktempdir() do base
        dir = joinpath(base, "metts.cylinder.tj.tjp", "outfiles.metts", "L.3.W.2", "J.0.4000", "t.3",
                       "t_prime.-0.3000", "filling.0.83333", "T.0.20000", "D.1000",
                       "tau.0.10.cutoff.1.0e-08.seed.7")
        mkpath(dir)
        write(joinpath(dir, "samples.txt"),
              "1: [1, 2, 3, 2, 3, 2]\n2: [2, 2, 3, 1, 3, 2]\n\n4: [3, 2, 2, 3, 1, 2]\n")
        p = ML.parse_ttj_path(dir)
        @test p.L == 3 && p.W == 2 && p.J == 0.4 && p.t == 3.0 && p.t_prime == -0.3
        @test p.T == 0.2 && p.D == 1000 && p.tau == 0.1 && p.cutoff == 1e-8 && p.seed == 7

        e = from_ttj_run(dir)
        @test nsites(e) == 6 && nsamples(e) == 3
        @test e.step == Int32[1, 2, 4]
        @test e.states[:, 1] == UInt8[0, 1, 2, 1, 2, 1]
        @test e.collapse_bases == ["X"]
        @test e.beta ≈ 5.0
        @test e.sector["n"] == 5
        @test e.lattice_name == "square.L3.W2.cyl"
        @test e.algorithm["seed"] == 7
        @test state_labels(e, 1) == ["Emp", "Up", "Dn", "Up", "Dn", "Up"]
        @test sort(lattice_couplings(lattice(e))) == ["J", "t", "t_prime"]
        @test lattice(e).coordinates[:, 2] == [0.0, 1.0]
    end
end

end
