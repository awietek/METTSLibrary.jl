# API

## Ensemble

```@docs
Ensemble
nsites
nsamples
lattice
validate
```

## Reading and writing

```@docs
write_ensemble
read_ensemble
```

## Lattice file functions

```@docs
Lattice
parse_lattice
read_lattice
lattice_couplings
square_lattice_toml
```

## Index and fetching

```@docs
build_index
load_index
ensembles
load
add_remote!
clear_remotes!
```

## Product states

```@docs
initial_states
state_labels
```

## Converters

```@docs
from_legacy_cpp
from_samples_txt
from_ttj_run
```

## Zenodo

```@docs
publish_zenodo
```

## Internals

Not exported; reachable as `METTSLibrary.<name>`.

```@docs
METTSLibrary.relpath_for
METTSLibrary.lattice_path
METTSLibrary.lat_to_toml
METTSLibrary.fetch_file
METTSLibrary.data_root
METTSLibrary.cache_dir
METTSLibrary.hf_token
METTSLibrary.parse_ttj_path
METTSLibrary.zenodo_files
```
