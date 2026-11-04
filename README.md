# AssemblyReview

An interactive reviewer for structural variants that sets a trio's short reads, long reads and assemblies side by
side, in the reference's coordinates and in the assembly's own. It takes neither the call nor the assembly as the
truth: it tests both on the raw reads and puts their disagreements first.

Status: planning (2026-11-04). Nothing is built yet.

The first case is HG01943's 1.09-Mb de novo tandem duplication on chr2, from the 1000 Genomes trios. Short and long
reads support the call, while HPRC's assembly of the child places the extra copy wrongly. A reviewer has to be able
to show both sides of that.

| | |
|---|---|
| [docs/brief.md](docs/brief.md) | What the reviewer must do: the questions, the views, the data to prepare, the architecture and the acceptance tests. |
| [cases/hg01943/](cases/hg01943/README.md) | The first case: findings, coordinates, its data and their public sources, the case page and its numbers. |
| [CLAUDE.md](CLAUDE.md) | Working rules for this repo. |

## Data

Each case keeps all the data it needs in its own `data/` folder. That data is region slices of public releases (1000
Genomes, HPRC and the 1000 Genomes ONT project) and the files made from them. The case's `data/MANIFEST.tsv` lists
every file with its size, md5, public source and the command that made it.
- Files under GitHub's 100-MB limit are committed. A larger one would be listed in the manifest with its public source
  and named in `.gitignore`.
- No file holds a local path: every header is cleaned, and every file is scanned.

All of the data is open access; see each source's own terms of use.
