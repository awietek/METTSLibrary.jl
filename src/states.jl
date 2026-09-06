# Sample j as 1-based ITensors state indices.
state_indices(e::Ensemble, j::Integer) = Int.(@view(e.states[:, j])) .+ 1

"""
    state_labels(e, j) -> Vector{String}

Sample `j` as ITensors state names, e.g. `["Emp", "Up", "Dn", ...]`, directly
usable as `MPS(sites, state_labels(e, j))`.
"""
state_labels(e::Ensemble, j::Integer) = e.local_states[state_indices(e, j)]

"""
    initial_states(e, n=nothing; basis=nothing, thin=1, rng=Random.default_rng()) -> Vector{Vector{Int}}

Return `n` product states from `e` as 1-based ITensors index vectors, the
same form `METTS.random_product_state` returns, ready to seed METTS chains.
Samples are drawn without replacement from the eligible ones, so `n` may not
exceed their number. With `n=nothing` all eligible samples are returned in
file order.

- `basis`: restrict to samples collapsed in this basis label, e.g. "Z".
- `thin`:  use only every `thin`-th step along each chain to reduce autocorrelation.
"""
function initial_states(e::Ensemble, n::Union{Nothing,Integer}=nothing;
                        basis=nothing, thin::Integer=1, rng::AbstractRNG=Random.default_rng())
    thin >= 1 || throw(ArgumentError("thin must be >= 1"))
    bcode = basis === nothing ? nothing : begin
        i = findfirst(==(basis), e.collapse_bases)
        i === nothing && throw(ArgumentError("basis '$basis' not among $(e.collapse_bases)"))
        UInt8(i - 1)
    end
    eligible = [j for j in 1:nsamples(e)
                if (bcode === nothing || e.basis[j] == bcode) && (thin == 1 || e.step[j] % thin == 0)]
    n === nothing && return [state_indices(e, j) for j in eligible]
    n <= length(eligible) ||
        throw(ArgumentError("requested $n states but only $(length(eligible)) eligible samples"))
    chosen = eligible[sort!(randperm(rng, length(eligible))[1:n])]
    return [state_indices(e, j) for j in chosen]
end
