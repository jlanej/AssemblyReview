# The reviewer: what was built

A static page (`index.html`, `reviewer/*.js`) that draws every panel from tables made by `tools/prepare_case.py`
out of the case's `data/` folder. It runs on GitHub Pages or any server that answers byte-range requests
(`serve.py`); opened from the file system, everything but the read panels still works. One exact igv.js release
(3.8.10, jsdelivr) is the only dependency, and only the read panels use it.

## The design

The case has five sides of evidence, and the page's core puts them on one axis so that they can be weighed at once:

| Comparison | What the board shows | Where to dig |
|---|---|---|
| Short reads (the guide) | the trio's depth in copies, bin by bin; the pairs and split reads at the junction | the read panels on GRCh38 |
| Long reads vs GRCh38 | HiFi depth by parent; the 31 junction molecules; the reads straight through each end | the junction-read table, dot plots, the read panels |
| Assembly vs GRCh38 | HPRC's own alignment: contig coordinate against GRCh38, each piece a diagonal, each gap a jump back | the picture, the pieces table |
| Assembly vs long reads | whose reads each piece takes once realigned to it; the junction reads splitting at M's ends; no read across a gap | depth by parent on HPRC's sequence; the assembly browser |
| Assembly vs assembly | each piece of the paternal scaffold, and U, aligned directly to the maternal chromosome: B and D identical to M | the hap1-against-hap2 dot plot |

Above the board, five cards give each comparison's verdict and numbers. Everything is hoverable, the board zooms
(drag) and resets (double-click), and a click sends the read panels to that place. Pieces are clickable everywhere
and open in the assembly browser.

Then, in order: the picture of what happened (the child's two chromosomes between HPRC's two, with bands from each
piece to the copy it holds); the read panels on GRCh38 (the trio's short reads coloured by pair orientation, the
child's HiFi reads coloured by parent, the mother's nanopore reads; the junction reads outlined); the 31 junction
reads with k-mer dot plots on demand (against GRCh38's two ends, GRCh38 joined as called, M's two ends, M joined as
called, and the paternal scaffold at gaps 1 and 4); the assembly tested on its own reads (depth by parent on HPRC's
sequence, the direct hap1-against-hap2 alignment, the pieces, the joins, every pair of pieces, HPRC's QC calls, and
igv.js on HPRC's local sequence and on the junction as called); where the copy was made (the mother's phase, each
copy's haplotype site by site, the reads that switch); every call set; every claim with its recount; the verdict row;
methods and provenance with every threshold.

Colours are the same everywhere: blue the father's, amber the mother's, pink the junction and whatever is missing or
differs, green and violet the mother's two haplotypes, grey unphased, black a gap.

## The data preparation

`tools/prepare_case.py cases/<case>` reads `case.json` (coordinates, file paths, thresholds) and `data/`, and writes
`data/reviewer/`:
- the tables (`reads.tsv`, `ont_reads.tsv`, `sites.tsv`, `mother_phase.tsv`, `depth_grch38.tsv`,
  `depth_assembly.tsv`, `pieces.tsv`, `joins.tsv`, `piece_pairs.tsv`, `qc_placed.tsv`, `calls.tsv`, `sr_pairs.tsv`,
  `self_matches.tsv`), the direct assembly alignment (`asm_vs_asm.paf`), the igv.js tracks (`tracks/*.bed`), the
  page's bundle (`bundle.js`, `sequences.js`) and `checks.tsv`;
- with `--align`, the junction as called (`realign/junction_as_called.fa`) and the HiFi reads aligned to it;
- a bgzipped GRCh38 chr2 with N outside the window (`ref/grch38.chr2.fa.gz`), so that the read panels need no
  reference from anywhere else.

Every threshold is in `case.json` and listed on the page. The script reproduces the case's `numbers.json` where the
original workup's definition is known (444 of 478 numbers, including every junction, phasing, piece, join, hole and
stretch count); the rest differ in the third decimal or in a definition the original did not publish (raw depth,
the phased-depth rule, the switch-read classes, the mother's last phase block), and each is listed with its reason in
`checks.tsv` and on the page.

## Where it departs from the brief

- **Reads are coloured by parent through a colour function keyed on the read name** (igv.js takes one), so the BAMs
  stay exactly as published: no tagged copies.
- **GRCh38 for the read panels is a padded chr2 from the repository**, not igv's hosted genome: no outside
  dependency, and it works offline with `serve.py`. A button adds RefSeq genes from UCSC's mirror when wanted.
- **The assembly against itself** is a sixth view the brief did not ask for: each gap-free piece aligned directly to
  the maternal chromosome, which shows B and D identical to M without any reference in between.
- **The verdict** is a row the reviewer copies or downloads, not a file the page writes: the page is static.
- Not made: the trio's short reads and the mother's nanopore reads on HPRC's sequence (two more alignments, about 90
  MB of BAM), and the breakpoint figures as separate panels, which the read panels cover.

## Adding a case

Make `cases/<case>/case.json` and `cases/<case>/data/` on the pattern of HG01943 (the same file roles; the piece
letters and local offsets in `assembly`), run `tools/prepare_case.py`, list the new files in the manifest, and open
`index.html?case=<case>`.
