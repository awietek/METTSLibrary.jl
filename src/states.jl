# Sample j as 1-based ITensors state indices.
state_indices(e::Ensemble, j::Integer) = Int.(@view(e.states[:, j])) .+ 1

"""
    state_labels(e, j) -> Vector{String}

Sample `j` as ITensors state names, e.g. `["Emp", "Up", "Dn", ...]`, directly
usable as `MPS(sites, state_labels(e, j))`.
"""
state_labels(e::Ensemble, j::Integer) = e.local_states[state_indices(e, j)]
