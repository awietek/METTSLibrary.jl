# Sharing

## Within the group

People on the institute cluster set `METTSLIBRARY_PATH` to the shared clone
and need nothing else: no account, no token, no network.

## Collaborators elsewhere

Add them as collaborators on the private Hugging Face dataset. They create a
read token, run `huggingface-cli login` or set `HF_TOKEN`, and the package
downloads what they load. They see the live library, including everything
added after you granted access.

## A frozen subset through Zenodo

For samples that should be citable, or for a collaborator who should not
need any relationship to your infrastructure, publish a Zenodo record. It
holds the chosen ensemble files, the lattice files they reference, and a
matching `index.toml`.

```julia
root  = ENV["METTSLIBRARY_PATH"]
paths = [x["path"] for x in ensembles(load_index(source = root); model = "tJ", Ly = 4, beta = 4.0)]

rec = publish_zenodo(paths, root;
        title       = "METTS product states, t-t'-J cylinder W4 L32, beta 4",
        description = "Collapsed product states and per-sample energies from ...",
        creators    = [Dict("name" => "Wietek, Alexander", "affiliation" => "MPI-PKS")],
        access      = "restricted",        # or "open"
        sandbox     = true)                # test on sandbox.zenodo.org first
rec["id"], rec["doi"]
```

A personal token goes in `ZENODO_TOKEN`, or `ZENODO_SANDBOX_TOKEN` for the
sandbox. A published record cannot be deleted and holds at most 100 files
and 50 GB. Zenodo stores files flat, so paths are flattened into names with
`__` as separator; the package undoes this on load.

The collaborator then does:

```julia
using METTSLibrary
add_remote!(zenodo = 1234567)          # plus token = "..." for a restricted record
e = load(ensembles(load_index(); beta = 4.0)[1])
```

For a restricted record they need a Zenodo account that you have granted
access, and their token in `ZENODO_TOKEN`.

## Other HTTPS sources

Any web server that exposes the library tree works as a source:

```julia
add_remote!("https://example.org/mettslibrary")        # base URL of the tree
add_remote!("https://huggingface.co/datasets/someone/fork/resolve/main"; token = "hf_...")
```

Remotes are tried in the order added, before the default Hugging Face
repository. `clear_remotes!()` returns to the default.
