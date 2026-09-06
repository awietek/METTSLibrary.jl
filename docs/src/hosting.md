# Hosting and copies

This page is for the maintainer. It records how the data repository is set
up and why, so that it can be moved when the time comes.

## Three copies

| copy | role |
|---|---|
| Hugging Face dataset `awietek/mettslibrary`, private | the Git server and off-site safe. Receives pushes; nothing reads from it day to day. |
| clone on the institute cluster | the primary. Samples are produced and ingested here; the group reads from it through `METTSLIBRARY_PATH`. |
| clone on the laptop | convenience for analysis and travel. |

All three are full Git clones with Git LFS, so each holds every file and the
complete history. Syncing is `git pull` and `git push`. There is no tool
between them and no partial copies to keep track of.

## Setting up a clone

```bash
brew install git-lfs          # or conda install -c conda-forge git-lfs on a cluster
git lfs install
git clone https://huggingface.co/datasets/awietek/mettslibrary
```

Pushing needs a Hugging Face token with write access, used as the password
with your username, or an SSH key registered with Hugging Face. On the
cluster keep the token in a file readable only by you. Make the clone
directory group-readable but not group-writable, so nobody edits the tree
by hand.

## Rules that keep it healthy

- Push after every ingest. The window in which new samples exist in one place
  only should be minutes, not weeks.
- Never modify or delete an ensemble file. Git LFS and Hugging Face store
  every version of a changed large file in full.
- Commit batches, not single files. Hugging Face degrades after a few
  thousand commits.
- The laptop's backup covers the laptop clone, which makes it a fourth,
  unregistered copy.

## Public snapshots

At each publication, publish the relevant subset to Zenodo (see
[Sharing](@ref)). That is the permanent, citable copy and it survives any of
the three working copies disappearing. Only what is published there is
guaranteed to outlive the grant, the institute affiliation and the hosting
company.

## Moving to another host

Hugging Face is the only service-specific piece, and it appears in exactly
one place in the package: the default remote URL. To move:

1. On a full clone, `git lfs fetch --all` to be sure every LFS object is local.
2. Create the empty repository on the new Git host with LFS (GitLab, GIN,
   Codeberg, a self-hosted Gitea or Forgejo).
3. `git remote set-url origin <new>`, `git push --all origin`, `git lfs push --all origin`.
4. Clone from the new host into a scratch directory and check a few files are
   real content, not pointers.
5. Change the default remote in `src/remotes.jl` and release the package.

Group members notice nothing, since they read from the cluster path. If the
new home is not a Git host, such as a bucket or a disk, the cluster clone
becomes the Git server, the laptop pulls from it over SSH, and the off-site
copy is an `rsync` of the whole clone including `.git`.

## Leaving the institute

The cluster clone is one copy of three. Copy it, or simply let the laptop and
Hugging Face copies stand, set up a clone at the new place, and point
`METTSLIBRARY_PATH` there. Nothing in the files or the package refers to the
institute.
