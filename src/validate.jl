"""
    validate(e::Ensemble) -> e

Check internal consistency. Throws an `ArgumentError` listing every problem
found. Checks:
- the lattice text parses and has as many sites as `states`
- every interaction refers to existing sites
- every coupling named in the lattice has a value in `parameters`
- `model`, `project` and `lattice_name` are usable as directory names
- label tables are non-empty and short enough for UInt8 indices
- every stored state and basis index refers to an existing label
- basis and observables have matching sample counts
- if a sector is given, every Z-basis sample has the declared particle numbers
"""
function validate(e::Ensemble)
    errs = String[]
    N, M = size(e.states)

    lat = nothing
    try
        lat = parse_lattice(e.lattice)
    catch err
        push!(errs, "lattice does not parse: " * sprint(showerror, err))
    end
    if lat !== nothing
        nsites(lat) == N ||
            push!(errs, "lattice has $(nsites(lat)) sites, states has $N")
        for (k, (cpl, typ, sites)) in enumerate(lat.interactions)
            all(0 .<= sites .< N) ||
                push!(errs, "interaction $k ($cpl $typ $sites) refers to sites outside 0:$(N-1)")
        end
        missing = [c for c in lattice_couplings(lat) if !haskey(e.parameters, c)]
        isempty(missing) ||
            push!(errs, "lattice couplings without a value in parameters: $(join(missing, ", "))")
    end
    # model, project and lattice_name each become one path component and are
    # read back out of the path, so a slash or a space in any of them would
    # change what the written file claims to be.
    for (what, s) in ("model" => e.model, "project" => e.project,
                      "lattice_name" => e.lattice_name)
        (isempty(s) || occursin(r"[/\\\s]", s)) &&
            push!(errs, "$what must be non-empty and contain no slashes or whitespace, got '$s'")
    end

    isempty(e.local_states) && push!(errs, "local_states is empty")
    length(e.local_states) > 255 && push!(errs, "more than 255 local states")
    isempty(e.collapse_bases) && push!(errs, "collapse_bases is empty")
    if !isempty(e.local_states) && M > 0
        mx = maximum(e.states)
        mx >= length(e.local_states) &&
            push!(errs, "state index $mx out of range for $(length(e.local_states)) local states")
    end
    if !isempty(e.collapse_bases) && M > 0
        mb = maximum(e.basis)
        mb >= length(e.collapse_bases) &&
            push!(errs, "basis index $mb out of range for $(length(e.collapse_bases)) bases")
    end

    length(e.basis) == M || push!(errs, "basis has length $(length(e.basis)), expected $M")
    for (k, v) in e.observables
        size(v, ndims(v)) == M ||
            push!(errs, "observable '$k' has last dimension $(size(v, ndims(v))), expected $M")
    end
    e.beta > 0 || push!(errs, "beta must be positive, got $(e.beta)")

    if !isempty(e.sector) && M > 0 && isempty(errs)
        _check_sector!(errs, e)
    end

    isempty(errs) || throw(ArgumentError("invalid Ensemble:\n  " * join(errs, "\n  ")))
    return e
end

# Sample particle numbers must match the declared sector.
#
# The two halves of this check have different scope, and conflating them hides
# real corruption. **Total charge `n` is conserved by every collapse basis** --
# an X collapse rotates spin, it does not move electrons -- so `n` is checkable
# in ANY basis. `nup`, `ndn` and `sz2` are conserved only by a Z collapse, so
# they stay restricted to Z samples.
#
# Checking `n` everywhere is what catches a dump's unwritten `ProductState`
# rows: the legacy code writes that dataset extensibly, and a row allocated but
# never filled reads back as all-empty, i.e. n=0. Eight such rows turned up in
# hubbard.triangular.metts, in an all-X project where the old Z-only check
# returned immediately and verified nothing at all.
function _check_sector!(errs, e::Ensemble)
    nup_of = [l in ("Up", "UpDn") ? 1 : 0 for l in e.local_states]
    ndn_of = [l in ("Dn", "UpDn") ? 1 : 0 for l in e.local_states]
    zidx = findfirst(==("Z"), e.collapse_bases)
    wantn = get(e.sector, "n", nothing)
    badn = 0; badz = 0
    for j in 1:nsamples(e)
        col = @view e.states[:, j]
        nup = sum(nup_of[s + 1] for s in col)
        ndn = sum(ndn_of[s + 1] for s in col)
        wantn === nothing || nup + ndn == wantn || (badn += 1)
        if zidx !== nothing && e.basis[j] == zidx - 1
            exp = Dict("nup" => nup, "ndn" => ndn, "n" => nup + ndn, "sz2" => nup - ndn)
            all(get(exp, k, v) == v for (k, v) in e.sector) || (badz += 1)   # unknown keys are not checked
        end
    end
    badn > 0 && push!(errs, "$badn sample(s) violate the declared particle number n=$wantn")
    badz > 0 && push!(errs, "$badz Z-basis sample(s) violate the declared sector $(e.sector)")
end
