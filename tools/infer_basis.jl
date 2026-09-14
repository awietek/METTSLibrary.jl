#  Infer the collapse basis from the states themselves.
#  X-basis collapse does not conserve total Sz; Z-basis does.
#  Legacy C++ encoding: 0 Emp, 1 Up, 2 Dn, 3 UpDn  =>  2*Sz = n(Up) - n(Dn)
using HDF5

const KMAX = 400            # samples read per file; enough to see Sz move

function szstats(path::String)
    try
        return h5open(path, "r") do f
            haskey(f, "ProductState") || return nothing
            d = f["ProductState"]
            sz = size(d)                       # (nsites, nsteps) in Julia order
            length(sz) == 2 || return nothing
            nsites, nsteps = sz
            nsteps == 0 && return nothing
            k = min(nsteps, KMAX)
            a = d[:, 1:k]                      # hyperslab: only what we need
            tot = Vector{Int}(undef, k)
            for j in 1:k
                up = 0; dn = 0
                for i in 1:nsites
                    v = Int(round(a[i, j]))
                    v == 1 && (up += 1)
                    v == 2 && (dn += 1)
                end
                tot[j] = up - dn
            end
            vals = Dict{Int,Int}()
            for t in tot; vals[t] = get(vals, t, 0) + 1; end
            modecount = maximum(values(vals))
            return (nsites, nsteps, k, length(vals), minimum(tot), maximum(tot), modecount)
        end
    catch
        return nothing
    end
end

listfile = ARGS[1]; outfile = ARGS[2]
open(outfile, "w") do io
    for line in eachline(listfile)
        isempty(line) && continue
        p = String(line)
        tag = occursin("updates.x", p) ? "x" : occursin("updates.z", p) ? "z" : "unknown"
        r = szstats(p)
        if r === nothing
            println(io, join([p, tag, "ERR", 0, 0, 0, 0, 0], "\t"))
        else
            (nsites, nsteps, k, ndist, mn, mx, mc) = r
            println(io, join([p, tag, "OK", nsites, nsteps, k, ndist, mn, mx, mc], "\t"))
        end
    end
end
