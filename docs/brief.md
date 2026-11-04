# The reviewer: brief

What to build first: one interactive page for the region of the first case, HG01943's chr2 duplication
([cases/hg01943/](../cases/hg01943/README.md)). It sets the trio's short reads, the long reads and the child's
assembly side by side. The design should carry over to later cases.

Read the case first: this brief builds on its findings (section 1), coordinates (2), data map (3) and traps (4).
Figure numbers below are the case page's (`cases/hg01943/reference/hg01943_chr2_dup.html`). File paths such as
`aln/...` or `realign/...` are inside `cases/hg01943/data/`.

## 1. Goal

One page where a reviewer can test every claim of the case's findings (`cases/hg01943/README.md`, section 1) on the data itself, in GRCh38's coordinates and the
assembly's own, with each read's parent shown everywhere, and then record a verdict. It takes neither the call nor
the assembly as the truth: it shows both models side by side and puts their disagreements first.

## 2. What a reviewer must be able to check

| Question | Evidence | View |
|---|---|---|
| Is there a junction, and in whom? | The trio's short reads at both breakpoints: split reads, and pairs in a duplication's orientation. The child's 31 HiFi junction reads. The mother's nanopore reads. | breakpoints; junction reads |
| How many copies, and whose? | The trio's short-read depth; HiFi depth by parent; allele shares at parent-opposite and mother-only sites. | overview |
| In every cell? | Junction against straight-through reads at the two ends: 31 / 64 / 69, against 1 : 2 : 2. | breakpoints |
| Unique, not a paralog? | Where the 240 pieces map (`aln/segment_pieces.paf`). | overview |
| What does the assembly hold, and does it stand up? | Pieces A-E, M and U; the joins and gaps; NucFlag and Flagger; the reads realigned to it, by parent; the reads at each join. | assembly |
| Where was the copy made? | The mother's heterozygous sites and phase blocks; each copy's alleles over the stretch; the HiFi reads that switch. | stretch |
| What do the other call sets say? | Yale's 16 sets, the review's 6, and the 1KGP record. | calls |

## 3. Views

1. **Overview** (GRCh38, WIN):
   - the segment and the junction;
   - the trio's short-read depth, in copies, and HiFi depth by parent;
   - the parent-opposite and mother-only sites, with the child's shares;
   - where each assembly piece aligns, and the missing copy;
   - the call sets, and a gene track if wanted.

   A click moves the read views to that place.
2. **Breakpoints side by side**: the segment's start and end, 20 kb either side.
   - The trio's short reads, pairs linked, split reads shown.
   - The child's HiFi reads, grouped and coloured by parent.
   - The mother's nanopore reads.
   - The junction reads marked.
3. **Junction reads**: the 31 as a table with each read's name, SMRT cell, instrument, length, two sides, MAPQ,
   parent and the sites behind that parent. Each read opens its alignments and a dot plot against:
   - GRCh38's two ends;
   - the assembly (M and the paternal scaffold);
   - the junction as called.
4. **The assembly**: hap1, hap2 and U in their own coordinates, linked to GRCh38.
   - The pieces as ribbons to GRCh38 (figures 1 and 7).
   - The gaps, and the NucFlag and Flagger calls.
   - The reads realigned to it, coloured by parent, with MAPQ 0 shown, not hidden.
   - Depth by parent (figure 10), and the reads at each join (figure 12).

   Selecting a read in any view finds its alignments in the others: a read has the same name in every BAM.
5. **The stretch** (about chr2:47.55-47.75 Mb): the mother's phase blocks, each copy's allele at her heterozygous
   sites (figure 13), and the reads that switch haplotype inside one molecule.
6. **Calls**: every record over the segment from every set, with the trio's genotypes where a set holds them, and what
   the record reflects: the junction, nothing, or the scaffold.
7. **Claims**: each claim of the case's findings with:
   - its value from numbers.json;
   - the threshold behind it (the page's section 10);
   - a link that opens the view filtered to the reads or sites behind it.

   Where the reviewer can recount a value from the data it has loaded, it shows both and flags any mismatch.
8. **Verdict**: a row in svreview's `verdicts.tsv`, with the columns candidate, trio, chrom, start, end, sv_type,
   verdict, note, reviewer and reviewed_at, so that the screen's review can read it back.

## 4. Data to prepare

The reads, alignments and calls are already files here (the case's section 3). What the original workup computed
in memory, the reviewer needs as tables. This repo's own data preparation should make them from the files in
`cases/hg01943/data/`, and they must reproduce `cases/hg01943/numbers.json`. The workup's thresholds are in the
case page's section 10; its script is not public. The tables:
- `reads.tsv`, one row per HiFi read:
  - its name, SMRT cell, length and MAPQ;
  - its parent (mother, father or unphased) and the sites behind it;
  - at the start, the end, and each join's two points: whether it is a junction read, straight through, or stops;
  - on the assembly: its piece, its MAPQ, and whether it splits at M's segment ends.
- `sites.tsv`, one row per SNV of the trio in WIN:
  - its position and kind (parent-opposite, mother-only, father-only, both heterozygous);
  - each member's genotype and allele depths, and the child's HiFi counts;
  - whether it lies inside the segment and the stretch, and whether the two maternal copies agree there.
- `mother_phase.tsv`: her nanopore phase blocks, and the haplotype of each of her heterozygous sites.
- `depth.tsv`: on GRCh38, short-read depth by member and HiFi depth by parent, in 10-kb bins; on the assembly, depth
  by parent in 5-kb bins.
- `assembly.bed`, in GRCh38 and in each local system: the pieces, the gaps, the missing copy, and the NucFlag and
  Flagger records.
- `calls.tsv`: one row per record over LWIN, per set.
- Copies of `hifi/HG01943.bam` and `realign/hifi_on_hprc.bam` tagged by parent: `HP:i:1` the father's and `HP:i:2`
  the mother's (HPRC's own hap1/hap2 numbering), unphased reads untagged. Any viewer can then colour and group reads
  by parent.
- `numbers.json`, as the claims' values (here, `cases/hg01943/numbers.json`).

New alignments worth making, seconds to minutes each on the region:
- **The junction as called.** GRCh38 from 30 kb before the segment's end to that end, then from its start to 30 kb
  after, with the HiFi reads aligned to it. The junction reads should align end to end; figure 4b shows the same in
  a dot plot.
- **The trio's short reads on HPRC's sequence.** The child's split reads and pairs should land at M's segment ends,
  and no pair should span a gap. The parents' reads show the same places without the junction.
- **The mother's nanopore reads on HPRC's sequence**, to show her molecules against the pieces she gave.

## 5. How to build it

Recommended:
- **Code** in this repo: an `index.html` and a few JavaScript modules. **Data**: the case's `data/` folder as it
  is, plus the reviewer's own tables in `cases/hg01943/data/reviewer/`, listed in the manifest and committed like
  the rest (each file under 100 MB).
- **Read panels**: igv.js, one exact release pinned, loaded from jsdelivr or vendored as its one file (no npm).
  - One browser instance per reference on the same page: GRCh38 (the analysis set by its public URL, or chr2
    rebuilt from `data/ref/grch38.fa` with N outside the window), `hprc_local.fa`, and the junction as called.
  - The short reads are already BAM here, with clean headers, which svreview's browser also needs.
  - Colour and group reads by the HP tag, and use igv.js's multi-locus view for the two breakpoints.
- **Custom panels**, in SVG or canvas: the overview, ribbons, depth, alleles, the stretch, the claims and the dot
  plots. Port build_case.py's figure code, and port its k-mer dot plot (16-mers, runs of 40 bp or more) to
  JavaScript so that any read can be plotted on demand.
- **Serve** it with a server that answers range requests (Python's `http.server` does not), such as the `serve`
  command of [svreview](https://github.com/jlanej/svreview-).

The alternative is to build long-read and assembly lanes into svreview itself, which already has its own BAM reader,
the review-set format and verdicts. That is the better long-term home, since every card whose child has an HPRC
assembly could use it, but it is a much larger change. Either way, this case and its numbers.json are the acceptance
test.

## 6. Acceptance: what it must reproduce

All from `cases/hg01943/numbers.json`:
- **Junction reads:** 31; 29 the mother's, none the father's, 2 unphased; both sides at MAPQ 60.
- **Straight through:** 64 at the start (33 the mother's, 31 the father's) and 69 at the end.
- **The mother's nanopore reads:** no junction read; 10 straight through each end.
- **Short reads:** the child has 6 split reads and 31 pairs, the parents none; 3.02 copies against the parents.
- **Phased HiFi:** 1.88 maternal copies, 0.89 paternal.
- **On the assembly:**
  - B holds 578 of the mother's reads and none of the father's;
  - no read crosses any join;
  - the junction reads split at M's segment ends: 29 at both, 30 at each.
- **The missing copy:** 186,866 bp.
- **The stretch:** the copies differ at 111 of 117 sites.

## 7. Keeping it clean

The data are public, and the files here carry no local paths: their headers are clean (see
`cases/hg01943/data/README.md`). Anything added must stay that way: a clean header, a scan of its raw and
decompressed bytes, and a manifest row naming its public source.

## Appendix: how the screen's review reads an assembly

The trio screen's review gives each card the line "an assembled haplotype holds it".
Its rules, adopted 2026-11-03, take neither the assembly nor the call as the truth:
1. A haplotype of the child holds the event: ok.
2. Every haplotype of the child runs straight through a DEL or INV: fails. The exception: if the child's own long
   reads carry the junction (two or more), it is one to follow up.
3. A haplotype runs straight through a DUP, a BND, an INS or an untyped card, or cannot be read across a DUP or an
   event of 100 kb or more: borderline, one to follow up, never counted against the call.
4. Both of the child's haplotypes hold it: borderline, one to follow up.
5. A parent's haplotype holds it, by its junction (or its insertion placed at the partner): fails, as inherited. An
   insertion near a lone breakpoint only: borderline.

A join across a gap is never a junction. HG01943 is the case rule 3 exists for. Its maternal haplotype runs straight
through and its paternal one cannot be read, while every test on the reads supports the call.
