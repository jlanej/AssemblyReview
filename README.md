# AssemblyReview

An interactive reviewer for structural variants that sets a trio's short reads, long reads and assemblies side by
side, in the reference's coordinates and in the assembly's own. It takes neither the call nor the assembly as the
truth: it tests both on the raw reads and puts their disagreements first.

The first case is HG01943's 1.09-Mb de novo tandem duplication on chr2, from the 1000 Genomes trios. Short and long
reads support the call, while HPRC's assembly of the child places the extra copy wrongly. The page shows both sides.

**The page:** `index.html` at the root of this repository, built for GitHub Pages. Once Pages is enabled for the
repository (Settings, Pages, deploy from the `main` branch, folder `/`), it is served at
`https://jlanej.github.io/AssemblyReview/`. Nothing else is needed: the page is static, its tables are in the
repository, and the read panels fetch the BAM and FASTA slices by byte range from the same site.

To run it locally:

```bash
python3 serve.py
```

then open `http://localhost:8000/`. (Python's own `http.server` ignores range requests, which the read panels need;
`serve.py` answers them.) Opened straight from the file system, every panel but the read panels still works.

## What the page shows

A short paper: an abstract; Table 1, each line of evidence by data type; Figure 1, the locus on one zoomable GRCh38
axis (copy number from short reads and from HiFi by parental origin, allele fractions, the assembly aligned to
GRCh38 as contig position against reference position, the parental origin of the reads each contig segment receives,
each segment's identity to hap2 by direct alignment, and the call-set records); Figure 2, the junction at base
resolution, every junction-spanning HiFi molecule with its parent-informative alleles, and dot plots of any one of
them against reference and assembly; Figure 3, the structure of the locus, each hap1 segment aligned directly to
hap2, and the proband's reads realigned to the assembly, with Table 2 of the contig segments; Figure 4, the maternal
haplotype carried by each of the two maternal copies and the reads that span the transitions; three interactive
igv.js views (GRCh38, the assembly, the junction as called); and methods with every threshold.
[docs/reviewer.md](docs/reviewer.md) describes the design.

| | |
|---|---|
| [index.html](index.html), [reviewer/](reviewer/) | The page and its modules: no build step, igv.js 3.8.10 from jsdelivr for the read panels, everything else drawn from the tables. |
| [tools/prepare_case.py](tools/prepare_case.py) | Makes a case's tables, tracks and data bundle from its `data/` folder, and checks every number against the case's `numbers.json`. |
| [cases/hg01943/](cases/hg01943/README.md) | The first case: findings, coordinates, its data and their public sources, the case page and its numbers. |
| [docs/reviewer.md](docs/reviewer.md) | What was built, how, and where it departs from the brief. |
| [docs/brief.md](docs/brief.md) | The original brief for the reviewer. |
| [CLAUDE.md](CLAUDE.md) | Working rules for this repo. |

## Data

Each case keeps all the data it needs in its own `data/` folder. That data is region slices of public releases (1000
Genomes, HPRC and the 1000 Genomes ONT project) and the files made from them, including the reviewer's own tables in
`data/reviewer/`. The case's `data/MANIFEST.tsv` lists every file with its size, md5, public source and the command
that made it.
- Files under GitHub's 100-MB limit are committed. A larger one would be listed in the manifest with its public source
  and named in `.gitignore`.
- No file holds a local path: every header is cleaned, and every file is scanned.

All of the data is open access; see each source's own terms of use.

## Rebuild a case

With samtools 1.22, bcftools 1.22, minimap2 2.30 and Python 3.14 (standard library only) on the path:

```bash
python3 -I tools/prepare_case.py cases/hg01943 --align
```

writes `cases/hg01943/data/reviewer/` and the junction-as-called alignment, and prints every number that differs from
`numbers.json` with the reason. `checks.tsv` keeps the full comparison.
