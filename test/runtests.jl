using METTSLibrary
using Test
using HDF5
using TOML
using Random
using Sockets
using JSON

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
        model="tJ", project="testproject", site_type="tJ",
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
    @test rel == joinpath("tJ", "testproject", "square.L4.W2.cyl", "J=0.4_t=3.0_t_prime=-0.3",
                          "ndn=2_nup=3", "T=00000.500000_beta=00002.000000", "chain01.h5")
    # a project is part of an ensemble's identity: the same physics run by two
    # people is two datasets, and they must not land on one path
    @test ML.relpath_for(with(e; project="other"); tag="a") !=
          ML.relpath_for(with(e; project="mine");  tag="a")
    @test_throws ArgumentError ML.relpath_for(with(e; project=""); tag="a")
    # a temperature whose inverse is not exactly representable must still give a
    # clean directory: beta = 1/0.0375 = 26.666..., and 1/beta round-trips to
    # 0.037500000000000006
    @test ML._beta_dir(with(e; beta = 1 / 0.0375)) == "T=00000.037500_beta=00026.666667"
    # fixed width, so a listing sorts numerically
    @test length(ML._tfield(0.0125)) == length(ML._tfield(80.0)) == 12
    @test sort([ML._tfield(x) for x in (80.0, 2.0, 9.0, 26.666667, 0.0125)]) ==
          [ML._tfield(x) for x in (0.0125, 2.0, 9.0, 26.666667, 80.0)]
    @test_throws ArgumentError ML.relpath_for(e; tag="")
    @test_throws ArgumentError ML.relpath_for(e; tag="a/b")

    # the tag is name=value like the rest of the path, built from whatever
    # algorithm parameters the run recorded, in TAG_FIELDS order
    alg = Dict{String,Any}("code" => "metts", "seed" => 3, "maxdim" => 1000,
                           "tau" => 0.1, "cutoff" => 1e-10, "nwarm" => 20)
    ea = with(e; algorithm=alg, collapse_bases=["X"])
    @test default_tag(ea) == "basis=X_maxdim=1000_tau=0.1_cutoff=1.0e-10_seed=3"
    # fields absent from algorithm are skipped, not written as empty
    @test default_tag(with(ea; algorithm=Dict{String,Any}("seed" => 3))) == "basis=X_seed=3"
    # unlisted keys stay out of the filename even though they are in the file
    @test !occursin("nwarm", default_tag(ea)) && haskey(ea.algorithm, "nwarm")
    # the field set is not fixed: another time evolution names its own
    @test default_tag(with(ea; algorithm=merge(alg, Dict{String,Any}("order" => 4)));
                      fields=["basis", "order", "seed"]) == "basis=X_order=4_seed=3"
    # a mixed-basis file is labelled by every basis it holds
    @test startswith(default_tag(with(ea; collapse_bases=["Z", "X"])), "basis=ZX_")
    # nothing to name is an error, not a nameless file
    @test_throws ArgumentError default_tag(with(e; algorithm=Dict{String,Any}()); fields=["seed"])
    # write_ensemble and relpath_for default to it
    @test ML.relpath_for(ea) == ML.relpath_for(ea; tag=default_tag(ea))
    @test ML.unflatten_path(ML.flatten_path(rel)) == rel
    # model, project and lattice are found by position, not stored in the file:
    # <model>/<project>/<lattice>/<params>/<sector>/<T>/<tag>.h5
    @test ML.lattice_dir("/lib/tJ/proj/sq/par/sec/T/b.h5") == "/lib/tJ/proj/sq"
    @test ML.layout_names("/lib/tJ/proj/sq/par/sec/T/b.h5") == ("tJ", "proj", "sq")
    # a slash in any of the three would silently change what the file claims
    for bad in (with(e; model="a/b"), with(e; project="a b"), with(e; lattice_name="a/b"))
        @test_throws ArgumentError ML.relpath_for(bad; tag="x")
        @test_throws ArgumentError validate(bad)
    end
    # the name=value codec round-trips, including keys with _ and negative values
    @test ML._kv_parse(ML._kv_dir(Dict("t_prime" => -0.3, "J" => 0.4))) ==
          ["J" => "0.4", "t_prime" => "-0.3"]
    @test ML._kv_parse("default") == Pair{String,String}[]
    # and they round-trip with what relpath_for builds
    @test ML.layout_names(joinpath("/lib", ML.relpath_for(e; tag="x"))) ==
          (e.model, e.project, e.lattice_name)
end

@testset "write, read, lattice sharing" begin
    e = synthetic_tj()
    mktempdir() do root
        rel = write_ensemble(root, e; tag="a")
        @test rel == ML.relpath_for(e; tag="a")
        path = joinpath(root, rel)
        @test isfile(path)
        @test_throws ErrorException write_ensemble(root, e; tag="a")     # append-only

        # the lattice lives in the lattice directory, inside the project
        latf = joinpath(root, "tJ", "testproject", "square.L4.W2.cyl", "square.L4.W2.cyl.toml")
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
            # one file is one run, so the seed lives in /algorithm, not a column
            @test !haskey(f, "chain")
            # names the path already carries are not duplicated inside the file
            ks = keys(attributes(f))
            for k in ("model", "project", "lattice_name", "lattice_sha256")
                @test !(k in ks)
            end
            # what IS stored is what the bytes cannot be read without
            for k in ("site_type", "local_states", "beta", "temperature", "collapse_bases")
                @test k in ks
            end
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
        @test length(filter(endswith(".toml"),
                            readdir(joinpath(root, "tJ", "testproject", "square.L4.W2.cyl")))) == 1
        # a different lattice under the same name is refused
        other = with(e; lattice=square_lattice_toml(4, 2; yperiodic=false, bonds=TJ_BONDS))
        @test_throws ErrorException write_ensemble(root, other; tag="a")
        # ... and accepted under its own name
        rel4 = write_ensemble(root, with(other; lattice_name="square.L4.W2.open"); tag="a")
        latf4 = joinpath(root, "tJ", "testproject", "square.L4.W2.open", "square.L4.W2.open.toml")
        @test isfile(latf4)
        @test read_ensemble(joinpath(root, rel4)).lattice_name == "square.L4.W2.open"
        # a missing lattice file is detected on read
        mv(latf4, latf4 * ".bak")
        @test_throws ErrorException read_ensemble(joinpath(root, rel4))
        @test_throws ErrorException ML.read_metadata(joinpath(root, rel4))
        # so is a second one: a lattice directory holds exactly one lattice
        write(latf4, e.lattice)
        write(joinpath(dirname(latf4), "extra.toml"), e.lattice)
        @test_throws ErrorException read_ensemble(joinpath(root, rel4))
        rm(joinpath(dirname(latf4), "extra.toml"))
        # a lattice that does not fit the data is caught by validate rather than
        # by a checksum: wrong site count here
        write(latf4, square_lattice_toml(6, 2; yperiodic=true, bonds=TJ_BONDS))
        @test_throws ArgumentError read_ensemble(joinpath(root, rel4))
        # A same-shape substitution is NOT detected, and that is the deliberate
        # cost of locating the lattice by position instead of pinning it with a
        # hash: this file has the same 8 sites and the same coupling names, only
        # different bonds. Paying it is what makes renaming a lattice a mv
        # rather than a rewrite of every file that references it.
        write(latf4, e.lattice)
        @test read_ensemble(joinpath(root, rel4)) isa Ensemble
        rm(latf4); mv(latf4 * ".bak", latf4)

        # Renaming a model, project or lattice is a mv: nothing inside a file
        # names them, so no file is rewritten and nothing can go stale.
        mv(joinpath(root, "tJ", "testproject"), joinpath(root, "tJ", "renamed.project"))
        moved = joinpath(root, "tJ", "renamed.project", splitpath(rel)[3:end]...)
        @test isfile(moved)
        r = read_ensemble(moved)
        @test r.project == "renamed.project"
        @test r.model == "tJ" && r.lattice_name == e.lattice_name
        @test r.states == e.states && r.beta == e.beta
        idx2 = build_index(root; write=false)
        @test all(x["project"] == "renamed.project" for x in idx2["ensembles"])
        mv(joinpath(root, "tJ", "renamed.project"), joinpath(root, "tJ", "testproject"))
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

@testset "serve: catalogue routes" begin
    mktempdir() do root
        write_ensemble(root, synthetic_tj(beta=2.0); tag="a")
        write_ensemble(root, synthetic_tj(beta=4.0, seed=2); tag="b")
        build_index(root)

        port = 8000 + (getpid() % 1000)
        task = @async METTSLibrary.serve(root; port=port, verbose=false)
        sleep(2)

        function get_(p)
            s = Sockets.connect("127.0.0.1", port)
            write(s, "GET $p HTTP/1.1\r\nHost: localhost\r\n\r\n")
            raw = read(s, String); close(s)
            i = findfirst("\r\n\r\n", raw)
            (head = raw[1:i[1]-1], body = raw[i[end]+1:end])
        end

        page = get_("/")
        @test occursin("200 OK", page.head)
        @test occursin("<title>METTS library</title>", page.body)

        s = JSON.parse(get_("/api/summary").body)
        @test s["nfiles"] == 2
        @test length(s["projects"]) == 1
        proj = s["projects"][1]["project"]
        lat  = first(keys(s["projects"][1]["lattices"]))

        l = JSON.parse(get_("/api/lattice?project=$proj&lattice=$lat").body)
        @test l["nsites"] > 0
        @test length(l["coordinates"]) == l["nsites"]
        # every bond must point at a site that exists, or the drawing is wrong
        @test all(all(0 .<= b["sites"] .< l["nsites"]) for b in l["bonds"])

        ens = JSON.parse(get_("/api/ensembles?project=$proj&lattice=$lat").body)
        @test !isempty(ens["runs"])
        @test issorted([r["temperature"] for r in ens["runs"]])

        d = JSON.parse(get_("/api/series?path=$(replace(ens["runs"][1]["path"], "=" => "%3D"))").body)
        @test d["nsamples"] > 0
        for (_, v) in d["observables"]
            @test length(v) == d["nsamples"]          # one value per sample
        end

        # a path from the browser is resolved through the index, never joined
        # onto the root, so it cannot escape the library
        for bad in ("../../../etc/passwd", "/etc/passwd")
            @test haskey(JSON.parse(get_("/api/series?path=$(replace(bad, "/" => "%2F"))").body),
                         "error")
        end
        @test occursin("404", get_("/api/nope").head)

        schedule(task, InterruptException(); error=true)
    end
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
        @test length(idx["ensembles"]) == 3

        # the root index is a manifest; the entries live per project
        manifest = TOML.parsefile(joinpath(root, "index.toml"))
        @test !haskey(manifest, "ensembles")
        @test only(manifest["projects"])["project"] == "testproject"
        @test only(manifest["projects"])["nfiles"] == 3
        @test isfile(joinpath(root, "tJ", "testproject", "index.toml"))
        pidx = TOML.parsefile(joinpath(root, "tJ", "testproject", "index.toml"))
        @test length(pidx["ensembles"]) == 3
        # each lattice is described once, not once per run that shares it
        @test length(pidx["lattices"]) == 2
        # and what the path already carries is not repeated in an entry
        for k in ("model", "project", "lattice_name", "parameters", "sector",
                  "temperature", "nsites", "couplings", "lattice", "provenance")
            @test !haskey(pidx["ensembles"][1], k)
        end

        # load_index expands it back to the flat shape callers expect
        parsed = load_index(source=root)
        @test length(parsed["ensembles"]) == 3
        first_entry = parsed["ensembles"][1]
        @test first_entry["couplings"] == ["t", "J", "t_prime"]
        @test first_entry["lattice"] ==
              joinpath("tJ", "testproject", "square.L4.W2.cyl", "square.L4.W2.cyl.toml")
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
            # a genuinely complex per-step scalar, as the Hubbard runs have in
            # `Polarization`: kept whole as two real observables
            f["Polarization"] = reshape(ComplexF64[0.5 + 0.25im, -0.5 + 0.0im,
                                                   0.0 - 0.75im, 1.0 + 1.0im,
                                                   -0.25 - 0.5im], 1, M)
        end
        lat = joinpath(dir, "chain.lat")
        write(lat, "[Coordinates]\n" * join(string.(0:N-1), "\n") *
                   "\n[Interactions]\n" * join(["HOP T $i $(i+1)" for i in 0:N-2], "\n") * "\n")

        kw = (model="tJ", project="testproject", site_type="tJ", temperature=0.5, basis="Z",
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
        # complex per-step scalars are split rather than reduced: neither the
        # modulus nor the real part alone would let the other be recovered
        @test !haskey(e.observables, "polarization")
        @test e.observables["polarization_re"] == [0.5, -0.5, 0.0, 1.0, -0.25]
        @test e.observables["polarization_im"] == [0.25, 0.0, -0.75, 1.0, -0.5]
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
        e3 = from_legacy_cpp(dump; lattice=lat, model="Hubbard", project="testproject",
                             site_type="Electron", beta=1.0,
                             parameters=Dict("T" => 1.0))
        @test e3.states[1, 1] == 3
        @test startswith(ML.relpath_for(e3; tag="x"), joinpath("Hubbard", "testproject", "chain"))
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

        # Real samples.txt files are written from a Dict and so arrive in hash
        # order, not chain order -- a kagome run reads 5, 56, 35, 55, ... while
        # covering 1..N exactly. The converter sorts by label, which is what
        # makes "sample position == step" true afterwards; without it the chain
        # would be scrambled and, with no step field, unrecoverable.
        let d = joinpath(dir, "scrambled")
            mkpath(d)
            write(joinpath(d, "samples.txt"),
                  "3: [3, 3, 3, 3, 3, 3]\n1: [1, 1, 1, 1, 1, 1]\n2: [2, 2, 2, 2, 2, 2]\n")
            s = from_samples_txt(joinpath(d, "samples.txt");
                                 lattice=square_lattice_toml(3, 2; yperiodic=true,
                                                             bonds=[("t", "HOP", :nn)]),
                                 lattice_name="square.L3.W2.cyl", model="tJ",
                                 project="testproject", site_type="tJ", temperature=0.5,
                                 parameters=Dict("t" => 1.0))
            # stored in chain order 1, 2, 3 -- not the file's 3, 1, 2
            @test s.states[1, :] == UInt8[0, 1, 2]
        end

        # the labels here are 1, 2, 4 -- not 1..N. They are not stored, but the
        # sample order is kept and the gap is reported rather than swallowed.
        e = @test_logs (:warn,) match_mode=:any from_ttj_run(dir; project="testproject")
        @test nsites(e) == 6 && nsamples(e) == 3
        @test e.states[:, 1] == UInt8[0, 1, 2, 1, 2, 1]
        @test e.states[:, 3] == UInt8[2, 1, 1, 2, 0, 1]   # the "4:" sample, third in order
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
