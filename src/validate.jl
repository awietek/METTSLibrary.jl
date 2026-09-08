"""
    validate(e::Ensemble) -> e

Check internal consistency. Throws an `ArgumentError` listing every problem
found. Checks:
- the lattice text parses and has as many sites as `states`
- every interaction refers to existing sites
- every coupling named in the lattice has a value in `parameters`
- `lattice_name` is usable as a directory name
- label tables are non-empty and short enough for UInt8 indices
- every stored state and basis index refers to an existing label
- basis, step and observables have matching sample counts
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
    (isempty(e.lattice_name) || occursin(r"[/\\\s]", e.lattice_name)) &&
        push!(errs, "lattice_name must be non-empty and contain no slashes or whitespace, got '$(e.lattice_name)'")

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
    length(e.step)  == M || push!(errs, "step has length $(length(e.step)), expected $M")
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

# Particle numbers of every Z-basis sample must match the declared sector.
function _check_sector!(errs, e::Ensemble)
    zidx = findfirst(==("Z"), e.collapse_bases)
    zidx === nothing && return
    nup_of = [l in ("Up", "UpDn") ? 1 : 0 for l in e.local_states]
    ndn_of = [l in ("Dn", "UpDn") ? 1 : 0 for l in e.local_states]
    expected(nup, ndn) = Dict("nup" => nup, "ndn" => ndn, "n" => nup + ndn, "sz2" => nup - ndn)
    bad = 0
    for j in 1:nsamples(e)
        e.basis[j] == zidx - 1 || continue
        col = @view e.states[:, j]
        exp = expected(sum(nup_of[s + 1] for s in col), sum(ndn_of[s + 1] for s in col))
        all(get(exp, k, v) == v for (k, v) in e.sector) || (bad += 1)   # unknown keys are not checked
    end
    bad > 0 && push!(errs, "$bad Z-basis sample(s) violate the declared sector $(e.sector)")
end
