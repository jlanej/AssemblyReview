# HG01943's 1.09-Mb duplication on chr2

The first case for this repo. A trio screen of the 1000 Genomes high-coverage genomes (kmer_denovo_filter in discovery
mode) calls a de novo tandem duplication of chr2:46,970,725-48,057,041 in HG01943 (father HG01941, mother HG01942).
HPRC's release 2 assembly of the child disagrees with the call. So every claim of both sides was tested on the raw
reads, and neither was taken as the truth. The call holds; the assembly is wrong in structure.

| In this folder | |
|---|---|
| this file | The findings (section 1), the coordinates in GRCh38 and in the assembly (2), the data (3), how the evidence was built and checked (4), and the open questions (5). |
| [data/](data/README.md) | Every file the case needs: region slices of the public releases and the files made from them. [data/MANIFEST.tsv](data/MANIFEST.tsv) gives each one's size, md5, public source and the command that made it. |
| [numbers.json](numbers.json) | Every number the case page states: the reviewer's acceptance values. |
| [case.json](case.json) | The case for the reviewer: coordinates, file roles, the assembly's contigs and local offsets, every threshold. `tools/prepare_case.py` reads it. |
| [reference/hg01943_chr2_dup.html](reference/hg01943_chr2_dup.html) | The original case page, self-contained; it opens offline. Its section 10 lists every source, every threshold, and a command to check each key claim. The interactive reviewer is the repository's `index.html`. |
| [reference/figures/](reference/figures/) | The page's 13 figures as SVG. |
| [lab_notes.md](lab_notes.md) | The case's lab-notebook entry (2026-11-03). |

What to build from it: [docs/brief.md](../../docs/brief.md).

## 1. Findings

The child carries three copies of the segment: one on the father's chromosome and two in tandem on the mother's.
Neither parent carries the duplication, and it is in every cell. It was made in the mother's germline, from both of
her chromosomes. HPRC's assembly has the right sequence in the wrong structure. It places two pieces of the second
maternal copy inside the paternal chromosome, at four gaps that no read crosses, and loses 187 kb of that copy,
the junction included.

### 1.1 The call holds

| Test | Result | numbers.json |
|---|---|---|
| The junction in short reads | The child has 6 split reads and 31 pairs (MAPQ 20 or more, median 60). Neither parent has any. | `split`, `pairs`, `parents_sr` |
| An independent call | 1KGP's GATK-SV set has HGSV_24204, a 1,086,317-bp DUP, PASS. The child is 0/1, both parents 0/0, and no other of its 3,202 genomes carries it. | `carriers` |
| Unique sequence | All 240 five-kb pieces of the segment, its flanks and its breakpoints map first to themselves at MAPQ 60. The 9 that match elsewhere are interspersed repeats. Neither breakpoint is in a segmental duplication; there is one 2,996-bp, 95.2% match inside the segment, 143 kb from its end. | `self_pieces`, `self_first`, `self_bp` |
| Copies, short reads | 3.02 inside the segment against the parents, bin by bin (father 2.00, mother 2.01). | `sr_rel` (x2) |
| Copies, phased HiFi | 1.88 maternal and 0.89 paternal, over the 29 bins where reads can be phased. | `hifi_mat`, `hifi_pat` |
| Allele balance | At the 221 parent-opposite sites inside, the child's maternal allele is 0.65 in short reads and 0.68 in HiFi. At the 163 outside it is 0.50. | `oh_in`, `oh_out` |
| HiFi junction reads | 31 reads, both sides at MAPQ 60, from 2 instruments (Sequel II 15, Revio 16) and 5 SMRT cells. All 29 that can be phased are the mother's. | `jr`, `jr_mat`, `jr_runs` |
| In every cell, not mosaic | There are 31 junction reads against 64 straight through the start and 69 through the end, so 0.19 of the reads cross the junction. A heterozygous tandem duplication in every cell gives 1/5 (1 junction : 2 : 2). | `through_s`, `through_e` |
| The mother's nanopore reads | 28 reach a breakpoint; none carries the junction. | `ont` |
| Yale's HiFi call sets | 6 of the 10 read-based callers write the junction to within 3 bp: Sniffles2, DeBreak, DELLY, SVIM, pbsv and Sawfish. | page section 6 |

The junction has 1 bp of microhomology (`mh`).

### 1.2 The assembly is wrong in structure, not in sequence

HPRC release 2 assembled HG01943 with trio-phased hifiasm (0.19.9): hap1 is the father's, hap2 the mother's.
Over the locus:
- **M**, the maternal chromosome (CM089362.1), holds the segment once, straight through.
- **U**, an unplaced maternal contig (JBHIIL010000075.1), holds 230 kb from inside the segment.
- The paternal chromosome (CM089354.1) is a scaffold of five pieces, **A-E**, that reads parts of the segment twice.
  - A, C and E carry the father's alleles: A has 352 of his own alleles and none of hers.
  - B and D are the mother's. B has 279 of her own alleles and none of his; D has 3 and none.
  - B and D are base-identical to M, over 368,852 and 300,439 bp.

What is wrong:
- **The joins.** The four joins between the pieces are the paternal chromosome's only gaps, 100 N each.
  - No HiFi read crosses any of them, while 79-99 run straight through each join's points on GRCh38.
  - NucFlag calls every join a MISJOIN and calls C a COLLAPSE. Flagger marks only the N bases.
- **The missing copy.** The second maternal copy survives only as fragments (B, U and D). 186,866 bp of it is
  missing, the junction included.
- **False duplications.** The paternal scaffold doubles two stretches where its pieces overlap: A and C over
  20,508 bp, and C and E over 17,767 bp.

Tested on the assembly itself: every primary HiFi read of the window was realigned to HPRC's own sequence, keeping
the parent its alleles gave it on GRCh38.
- The 31 junction reads still split, at M's two segment ends. 29 split at both ends; the other 2 put a short side on
  one of the father's pieces at MAPQ 0. None runs through.
- B takes 578 of the mother's reads and none of the father's, all at MAPQ 0, because B is identical to M.
- No read crosses a gap.

The call sets built from the assembly inherit the scaffold:
- PAV writes four paternal insertions (369, 21, 21 and 300 kb);
- cuteSV-asm writes four paternal duplications;
- SVision-pro-asm writes breakends at the joins;
- dipcall and svim-asm write nothing.

### 1.3 Where the copy was made

- The child's two maternal copies agree at 118 of the mother's heterozygous sites before chr2:47,590,351 and at 133
  after chr2:47,704,313. Between those points they differ at 111 of 117. There the child's share of her allele is
  1/3, which takes three copies.
- The child inherited one maternal chromosome, and a copy made after conception would repeat it. So the copy was
  made in the mother's germline, from both her chromosomes.
- Her nanopore reads phase her heterozygous sites into 4 blocks, with no conflicting links. Read against the block
  across the stretch (chr2:47,557,874-47,730,934):
  - before the stretch, the copies carry her haplotype 2 (15 sites; in the patchy kilobases, 7 read both haplotype 1
    and 8 one of each);
  - inside it, one copy carries each haplotype (109 of 115 sites);
  - after it, both carry haplotype 1 (22 of 22).
- So each copy changes from her haplotype 2 to haplotype 1, at a different place: one near chr2:47,590,351, the
  other between chr2:47,704,313 and 47,713,612. HiFi reads show both switches inside single molecules.
- Still open: two exchanges about 120 kb apart in one meiosis would be unusual, and one of them may be a template
  switch made while the copy was made. Which copy comes first, the reads cannot say.

### 1.4 For the review

The review's rule 3 kept this true call standing. A haplotype that runs straight through, or cannot be read, at a
duplication or at an event of 100 kb or more is borderline and flagged to follow up, never counted against the call.

The page proposes four changes; none is applied:
- read the pieces of a haplotype that reads the segment twice by the trio's alleles, and count its copies against
  the call;
- give an assembly-based call at a gap no weight;
- add phased long-read depth to the long-read line for large duplications;
- mark this card "followed up", with a link, instead of leaving it on the list.

## 2. Coordinates

GRCh38, 1-based and inclusive, unless a table says otherwise.

| | |
|---|---|
| The segment | chr2:46,970,725-48,057,041 (1,086,317 bp) |
| The junction | the segment's last base, chr2:48,057,041, joined to its first, chr2:46,970,725 |
| The screen's card | chr2:46970724-46970852. Its other end, chr2:48056899-48057041, shares the bundle. Stage 4e: tier 1, rank 187; the child has 6 split reads and 33 discordant pairs, depth ratio 1.55. |
| Short-read and assembly window (`WIN`) | chr2:46,400,001-48,650,000 |
| Long-read and call-set window (`LWIN`) | chr2:46,850,001-48,180,000 |
| The stretch where the copies differ | chr2:47,590,351-47,704,313 |

### 2.1 The assembly's pieces

These come from HPRC's own Winnowmap alignment of each haplotype to GRCh38 (`pieces` in numbers.json). The last
column counts each parent's own alleles that the piece carries inside the segment.

| Piece | Haplotype | Contig | On the contig | On GRCh38 (chr2) | Parent | Mother's / father's own alleles |
|---|---|---|---|---|---|---|
| A | hap1, paternal | CM089354.1 | 44,643,737-47,396,431 | 44,588,993-47,340,307 | father | 0 / 352 |
| B | hap1 | CM089354.1 | 47,396,532-47,765,590 | 46,971,835-47,340,686 | mother | 279 / 0 |
| C | hap1 | CM089354.1 | 47,765,691-48,221,141 | 47,319,800-47,774,341 | father | 0 / 122 |
| D | hap1 | CM089354.1 | 48,221,250-48,522,660 | 47,755,434-48,055,872 | mother | 3 / 0 |
| E | hap1 | CM089354.1 | 48,522,761-49,322,940 | 47,756,575-48,555,799 | father | 0 / 7 |
| M | hap2, maternal | CM089362.1 | 44,613,400-48,584,116 | 44,588,649-48,555,799 | mother | 422 / 0 |
| U | hap2, unplaced | JBHIIL010000075.1 | 1-230,135 (whole) | 47,426,809-47,656,968 | mother | 27 / 0 |

### 2.2 The joins, the missing copy and HPRC's QC calls

| Join | Pieces | Gap on CM089354.1 | GRCh38 points: left piece's end / right piece's start | HiFi reads straight through each point | Reads crossing it |
|---|---|---|---|---|---|
| 1 | A to B | 47,396,432-47,396,531 | 47,340,307 / 46,971,835 | 80 / 94 | 0 |
| 2 | B to C | 47,765,591-47,765,690 | 47,340,686 / 47,319,800 | 79 / 84 | 0 |
| 3 | C to D | 48,221,142-48,221,241 | 47,774,341 / 47,755,434 | 80 / 81 | 0 |
| 4 | D to E | 48,522,661-48,522,760 | 48,055,872 / 47,756,575 | 99 / 81 | 0 |

At join 3, 108 bases are left unaligned: the 100 N and 8 bases beside them that no alignment record holds.

The second copy is missing over chr2:46,970,725-46,971,834, 47,340,687-47,426,808, 47,656,969-47,755,433 and
48,055,873-48,057,041: 186,866 bp in all (`holes`).

NucFlag's calls, placed on GRCh38 (`flags_placed`):
- MISJOIN at the ends of A (chr2:47,324,847-47,340,307), B (47,322,847-47,340,686), C (47,765,076-47,774,341) and
  D (48,038,699-48,055,872), at 47,744,602-47,744,619 inside C, and at both ends of U;
- COLLAPSE in C (chr2:47,686,528-47,711,595). Where the second copy is missing over 47,656,969-47,755,433, its reads
  land on C (53 of 55 informative sites), and this collapse lies inside that stretch.

### 2.3 The case files' own coordinates

The region files name each sequence by its region, so positions in them are offsets:

| File | Sequence | Is | Conversion |
|---|---|---|---|
| `ref/grch38.fa` | `chr2:46400001-48650000` | GRCh38 over WIN | GRCh38 = local + 46,400,000 |
| `asm/pat.region.fa` | `HG01943#1#CM089354.1:46454978-49418103` | the paternal chromosome | contig = local + 46,454,977 |
| `asm/mat.region.fa` | `HG01943#2#CM089362.1:46426428-48679377` | M | contig = local + 46,426,427 |
| | `HG01943#2#JBHIIL010000075.1:1-231135` | U | contig = local. The name runs past U's end: U is 230,135 bp. |
| `realign/hprc_local.fa`, `realign/hifi_on_hprc.bam` | `hap1` | CM089354.1:46,805,838-49,047,195 | contig = local + 46,805,837 |
| | `hap2` | CM089362.1:46,777,285-48,308,368 | contig = local + 46,777,284 |
| | `U` | JBHIIL010000075.1, whole | contig = local |

On `hprc_local.fa`:
- M's copy of the segment runs from hap2 220,695 to 1,308,103 (CM089362.1:46,997,979-48,085,387);
- the gaps' N runs are at hap1 590,595-590,694, 959,754-959,853, 1,415,305-1,415,404 and 1,716,824-1,716,923.

HPRC names its contigs in PanSN form (`HG01943#1#CM089354.1`) in its FASTAs, its alignment records and its QC tracks.
numbers.json and the page use the bare accession.

## 3. Data

Everything the case needs is in [data/](data/README.md): 102 files, 294 MB, each under GitHub's 100-MB limit. Two are
over GitHub's 50-MB warning: `realign/hifi_on_hprc.bam` (78.9 MB) and `hifi/HG01943.bam` (67.6 MB). Every file is a
region slice of a public release, or was made from such slices; [data/MANIFEST.tsv](data/MANIFEST.tsv) gives each
file's size, md5, public source and the command that made it. All of the data is open access, from 1000 Genomes,
HPRC and the 1000 Genomes ONT project. Paths below are inside `data/`.

### 3.1 What is in data/

| Path | Holds | From |
|---|---|---|
| `sr/{child,father,mother}.bam` | the trio's short reads over WIN, decoded from the CRAM slices so that no reference is needed to read them | the 1000 Genomes high-coverage CRAMs (NYGC), AWS mirror |
| `hifi/HG01943.bam` | the child's HiFi reads over LWIN | HPRC's reads as Yale aligned them to GRCh38_no_alt |
| `ont/mother.bam` | the mother's nanopore reads over LWIN | 1000 Genomes ONT (Vienna) |
| `asm/{pat,mat}.winnowmap.bam` | HPRC's own alignment of each haplotype to GRCh38, over WIN | HPRC release 2 |
| `asm/{pat,mat}.region.fa` | HPRC's published sequence of each contig aligned in WIN, 1 kb either side | HPRC release 2 FASTAs |
| `qc/{hap1,hap2}.{gaps,nucflag,flagger_hifi,flagger_ont}.bed` | HPRC's gaps, NucFlag and Flagger rows on the child's contigs over the locus | HPRC release 2 annotation and QC |
| `calls/yale.<caller>.vcf` | Yale's 16 call sets of the child over LWIN: 10 from its HiFi reads, 5 from its assembly, 1 merged | HPRC release 2, Yale |
| `calls/<i>.<set>.vcf` + `.json` | the review's six published sets over LWIN, with the trio's columns where a set holds them: 0 1KGP, 1 HPRC, 2 and 3 HGSVC3 (insdel, inv), 4 ONT SVIM-asm, 5 1KGP-asm. Each `.json` names the URL and the columns. | the sets' releases |
| `calls/1KGP.event_all_samples.vcf` | the 1KGP SV set's record of the event, with every genome's column | 1KGP GATK-SV |
| `ref/grch38.fa` | GRCh38 over WIN | the 1000 Genomes GRCh38 analysis set |
| `aln/trio.sr.bcf` | the trio's SNVs over WIN | made here: bcftools mpileup and call on `sr/` |
| `aln/{sr,hifi}.bins.bed` | the 10-kb depth bins over WIN and LWIN | made here |
| `aln/segment_pieces.paf` | the 240 pieces against all of GRCh38 | made here: minimap2 -x asm20 |
| `realign/hprc_local.fa` | HPRC's sequence over the locus as one FASTA: hap1, hap2 and U | made here from `asm/*.region.fa` |
| `realign/hifi_on_hprc.bam` | the 6,765 primary HiFi reads of LWIN, realigned to `hprc_local.fa` | made here: minimap2 -ax map-hifi --eqx -Y |
| `realign/junction_as_called.fa`, `realign/hifi_on_junction.bam` | the junction as called (30 kb up to the segment's end joined to its first 30 kb), and the HiFi reads aligned to it: all 31 junction reads cross the join unbroken | made here by `tools/prepare_case.py --align` |
| `ref/grch38.chr2.fa.gz` | GRCh38 chr2 for the page's read panels: the window in place, N elsewhere, bgzipped | made here from `ref/grch38.fa` |
| `reviewer/` | the reviewer's tables, the direct piece-to-M alignment, the igv.js tracks, the page's data bundle and `checks.tsv`, every number of numbers.json against its recount | made here by `tools/prepare_case.py` |

With them are their indexes (`.bai`, `.csi`, `.fai`). The headers are clean: no `@PG` lines, no `UR` tags, and no
header line naming a file path ([data/README.md](data/README.md) lists what was removed). The read records are exactly
as published: each BAM holds as many records as its source.

### 3.2 The public sources

| What | Where |
|---|---|
| HPRC release 2, the paternal haplotype (hap1) | `https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01943/assemblies/release2/HG01943_pat_hprc_r2_v1.0.1.fa.gz` |
| its alignment to GRCh38 | `.../release2/alignments/winnowmap/GRCh38_no_alt/HG01943.GRCh38_no_alt.hap1.sorted.bam` |
| the maternal haplotype (hap2) and its alignment | `.../release2/HG01943_mat_hprc_r2_v1.0.1.fa.gz` and `.../HG01943.GRCh38_no_alt.hap2.sorted.bam` |
| the gaps | `.../release2/annotation/chrom_assignment/HG01943_{pat,mat}_hprc_r2_v1.0.1.gaps.bed` |
| NucFlag | `https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA--HPRC_R2_QC_NUCFLAG/HG01943/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/HG01943_{pat,mat}_hprc_r2_v1.0.1.nucflag.bed` |
| Flagger, from HiFi and from ONT reads | `https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/HG01943/hprc_r2/assembly_qc/hmm_flagger/v1.2.0_{hifi,ont}/...`. The full file names are in `data/qc/*.bed` and the manifest. |
| the child's HiFi reads | `https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/samples/HG01943/alignments/GRCh38_no_alt/hifi_reads/HG01943.GRCh38_no_alt.sorted.bam` |
| Yale's call sets of the child | `.../samples/HG01943/variant_calls/GRCh38_no_alt/<caller>/HG01943.GRCh38_no_alt.<caller>.vcf.gz`; for svdss the last part is `svdss.genotyped`, for merged `autosomes.merged` |
| the mother's nanopore reads | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG01942.hg38.cram` |
| the trio's short reads | child `https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988955/HG01943.final.cram`; father `.../1000G_2504_high_coverage/data/ERR3242006/HG01941.final.cram`; mother `.../data/ERR3242007/HG01942.final.cram` |
| the 1KGP SV set | `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz` |
| HPRC release 2's Minigraph-Cactus VCF | `https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.wave.vcf.gz` |
| HGSVC3's PAV calls: insertions and deletions, and inversions | `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_sv_{insdel,inv}_sym_HGSVC2024v1.0.vcf.gz` |
| the 1000 Genomes ONT project's SVIM-asm calls | `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.1/svim-asm-hg38/svim.asm.hg38.bcf` |
| 1KGP_PacBio_WGS's integrated HiFi call set ("1KGP-asm") | `https://s3.amazonaws.com/1000g-ont/1KGP_PacBio_WGS/SV_Callset/Integrated_Callset/GRCh38/disco_truvari_collapsed.deduped.vcf.gz` |
| the GRCh38 analysis set | `https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa` |

### 3.3 Long reads: who has what

Of the five long-read sets the review reads (the 1000 Genomes ONT project, HPRC's HiFi, 1KGP-PacBio, 1kGP-LRSC ONT
and HGSVC1 CLR):
- the child has HPRC's HiFi reads only;
- the mother has the 1000 Genomes ONT project's nanopore reads only;
- the father has none.

By its file name, HPRC's Flagger ONT track was made from the child's own R9.4.1 nanopore reads, aligned to its assembly. The
2026-10-30 survey found HPRC's nanopore reads published as raw reads only, and it was not checked whether those
alignments are published. Aligned to GRCh38, those reads would give the child a second long-read technology.

## 4. How the evidence was built and checked

- **The slices.** Each file in `data/` was read by region from its public source (HTTP range reads), by the command
  in its manifest row, and the derived files were made from those slices by the commands listed there. To fetch a
  file again, run its row's command against its source, then clean its header the same way.
- **The numbers.** Every number in [numbers.json](numbers.json) and every figure of the case page was computed from
  these same files by the original workup's script, `build_case.py`, which is not public. Its thresholds are listed in
  the page's section 10. On 2026-11-04 a rebuild from a copy of the files reproduced the page and its 13 figures byte
  for byte, apart from the build date, and numbers.json came out equal as data.
- **Tools:** samtools 1.22.1, bcftools 1.22, minimap2 2.30-r1287 and Python 3.14.
- **Check a claim yourself.** The page's section 10 gives commands against the public files, and they run as well on
  the files here. For example, the 31 junction reads:

  ```bash
  cd cases/hg01943/data
  samtools view -q 20 -F 0x704 hifi/HG01943.bam chr2:46970725-46970725 | awk '$4 == 46970725 && /SA:Z:chr2,4805?[0-9]/' | wc -l
  ```

- **The reviewer's recount.** `tools/prepare_case.py` recomputes every number here from these files and writes
  `data/reviewer/checks.tsv`: 444 of 478 equal numbers.json; the rest differ in the third decimal or in a definition
  the original did not publish, each with its reason.

Traps:
- The region files' positions are offsets (section 2.3).
- B and D are identical to M, so on the assembly every read on B or D maps at MAPQ 0 (622 of 622 on B). A MAPQ
  filter hides the mother's second copy.
- The short reads are BAM here. The original CRAMs decode only against the 1000 Genomes analysis set; IGV's hosted
  hg38 fails their MD5 check.
- When fetching again, EBI's range reads sometimes hang: retry.

## 5. Open questions

- Which maternal copy comes first? And are the two haplotype switches two meiotic exchanges, or one exchange and a
  template switch?
- How common is HPRC's scaffolding error at large tandem duplications? Rule 3's follow-up list holds the other cards
  where an assembly runs straight through, or cannot be read, at a large event.
- The father has no long reads, and the child has HiFi only. The child's nanopore reads at HPRC (section 3.4) are
  the nearest second technology.
