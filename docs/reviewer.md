# The reviewer: what was built

A static page (`index.html`, `reviewer/*.js`) in the form of a short paper: a title, an abstract, one evidence
table, four numbered figures with lettered panels and captions, three supplementary interactive views, methods, and a
data note. Every figure is drawn from tables that `tools/prepare_case.py` rebuilds from the case's `data/` folder.
It runs on GitHub Pages or any server that answers byte-range requests (`serve.py`); opened from the file system,
everything but the interactive views still works. One pinned igv.js release (3.8.10, jsdelivr) is the only
dependency, used only by the interactive views.

## The page

| Element | What it shows |
|---|---|
| Table 1 | Each line of evidence by data type (short reads of the trio, HiFi of the proband, ONT of the mother, the HPRC assembly): copy number, junction, parental origin of the extra copy, zygosity, reads realigned to the assembly, other call sets. |
| Figure 1 | The locus on one GRCh38 axis, zoomable: (a) short-read copy number of the trio with the junction pairs; (b) HiFi copy number by parental origin with the junction-spanning reads; (c) maternal allele fraction at parent-opposite SNVs; (d) maternal-specific allele fraction at mother-only SNVs; (e) the assembly aligned to GRCh38 as contig position against reference position; (f) parental origin of the reads each contig segment receives; (g) identity of each hap1 segment and U to hap2 by direct alignment; (h) call-set records. |
| Figure 2 | The junction: (a) the sequence with its microhomology; (b) all junction-spanning HiFi reads aligned at the junction with their parent-informative alleles; (c) dot plots of one read (selectable) against GRCh38, GRCh38 joined as called, hap2, and hap2 joined as called. |
| Figure 3, Table 2 | The assembly: (a) the structure of the locus to scale, reads against assembly; (b) each gap-free hap1 segment and U aligned directly to hap2; (c) the proband's reads realigned to the assembly, depth by parental origin, gaps, NucFlag. Table 2 lists each contig segment with its alleles, reads, identity to hap2 and join support. |
| Figure 4 | The two maternal copies: (a) copies carrying the maternal-specific allele along the segment; (b) the maternal haplotype each copy carries within the mother's ONT phase block; (c, d) HiFi reads spanning each haplotype transition. |
| Interactive views | igv.js on GRCh38 (the trio at both breakpoints), on HPRC's sequence (the realigned reads) and on the junction as called. Reads coloured by parental origin through a colour function keyed on the read name, so the BAMs stay as published. |
| Methods | Nine numbered paragraphs with every threshold, and the recount against the case record. |

Colours are the same throughout: paternal blue, maternal orange, junction and absent sequence magenta, maternal
haplotype 1 teal and haplotype 2 purple, unphased grey, gaps black.

## The data preparation

`tools/prepare_case.py cases/<case>` reads `case.json` (coordinates, file paths, thresholds) and `data/`, and writes
`data/reviewer/`:
- the tables (`reads.tsv`, `ont_reads.tsv`, `sites.tsv`, `mother_phase.tsv`, `depth_grch38.tsv`,
  `depth_assembly.tsv`, `pieces.tsv`, `joins.tsv`, `piece_pairs.tsv`, `qc_placed.tsv`, `calls.tsv`, `sr_pairs.tsv`,
  `self_matches.tsv`), the direct assembly alignment (`asm_vs_asm.paf`), the igv.js tracks (`tracks/*.bed`), the
  page's bundle (`bundle.js`, `sequences.js`) and `checks.tsv`;
- with `--align`, the junction as called (`realign/junction_as_called.fa`) and the HiFi reads aligned to it;
- a bgzipped GRCh38 chr2 with N outside the window (`ref/grch38.chr2.fa.gz`), so that the interactive views need
  no reference from anywhere else.

Every threshold is in `case.json` and listed in the page's Methods. The script reproduces the case's `numbers.json`
where the original workup's definition is known (444 of 478 numbers, including every junction, phasing, segment,
join, hole and interval count); the rest differ in the third decimal or in a definition the original did not publish
(raw depth, the phased-depth rule, the switch-read classes, the mother's last phase block), and each is listed with
its reason in `checks.tsv`.

## Where it departs from the brief

The brief (`docs/brief.md`) described a reviewer with many views and a verdict file. The page was redesigned as a
figure-led scientific presentation instead:
- the sides of the evidence are weighed at once (Table 1 and Figure 1) rather than visited view by view;
- the assembly is compared with itself by direct alignment (Figure 3b), which the brief did not ask for and which
  shows B and D identical to hap2 without any reference in between;
- reads are coloured by parent through a colour function, not through tagged BAM copies;
- GRCh38 for the interactive views is a padded chr2 from the repository, not a hosted genome;
- the verdict row, the claims scorecard and the breakpoint figures as separate panels were dropped; the recount
  against `numbers.json` lives in `checks.tsv` and the Methods;
- not made: the trio's short reads and the mother's nanopore reads realigned to HPRC's sequence.

## Adding a case

Make `cases/<case>/case.json` and `cases/<case>/data/` on the pattern of HG01943 (the same file roles; the segment
letters and local offsets in `assembly`), run `tools/prepare_case.py`, list the new files in the manifest, and open
`index.html?case=<case>`.
