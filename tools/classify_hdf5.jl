using HDF5

# Classify one HDF5 file: does it hold collapsed product states, and if so how
# many sites/steps and how many per-step scalar observables.
# Columns: path, status, bytes, nsites, nsteps, nscalars, ndatasets, names
function classify(path::String)
    bytes = try filesize(path) catch; 0 end
    try
        return h5open(path, "r") do f
            ks = collect(keys(f))
            nd = length(ks)
            if "ProductState" in ks
                d = f["ProductState"]
                sz = size(d)                      # Julia order: (nsites, nsteps)
                si = length(sz) >= 1 ? sz[1] : 0
                st = length(sz) >= 2 ? sz[2] : 0
                nsc = 0
                for k in ks
                    k == "ProductState" && continue
                    try
                        o = f[k]
                        o isa HDF5.Dataset || continue
                        s = size(o)
                        length(s) == 2 && s[1] == 1 && s[2] == st && (nsc += 1)
                    catch
                    end
                end
                return (path, "STATES", bytes, si, st, nsc, nd, join(ks, ","))
            else
                return (path, "NOSTATES", bytes, 0, 0, 0, nd, join(ks, ","))
            end
        end
    catch err
        return (path, "ERROR", bytes, 0, 0, 0, 0, string(typeof(err)))
    end
end

listfile = ARGS[1]
outfile  = ARGS[2]
open(outfile, "w") do io
    for (i, line) in enumerate(eachline(listfile))
        isempty(line) && continue
        r = classify(String(line))
        println(io, join(r, "\t"))
        i % 2000 == 0 && flush(io)
    end
end
