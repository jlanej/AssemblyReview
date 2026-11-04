# HG01943: the data

Everything the case needs: 66 files, 270 MB, each under GitHub's 100-MB limit. Every file is a region slice of a
public release, or was made from such slices. [MANIFEST.tsv](MANIFEST.tsv) lists each file with its size, md5, what it
holds, its public source and the command that made it.

## Layout

There are two windows, in GRCh38 (1-based):
- **WIN**, chr2:46,400,001-48,650,000: the short reads and the assembly's alignments;
- **LWIN**, chr2:46,850,001-48,180,000: the long reads and the call sets.

| Folder | Holds |
|---|---|
| `sr/` | the trio's short reads over WIN, as BAM: `child.bam`, `father.bam`, `mother.bam` |
| `hifi/` | the child's HiFi reads over LWIN |
| `ont/` | the mother's nanopore reads over LWIN |
| `asm/` | HPRC's alignment of each of the child's haplotypes to GRCh38 over WIN, and the haplotypes' published sequence there |
| `qc/` | HPRC's gaps, NucFlag and Flagger rows on the child's contigs over the locus, as BED |
| `calls/` | Yale's 16 call sets of the child, the review's 6 published sets and the 1KGP record of the event, as VCF; one `.json` per published set names what was read |
| `ref/` | GRCh38 over WIN. Its one sequence is named `chr2:46400001-48650000`, so positions in it are offsets: GRCh38 = position + 46,400,000. |
| `aln/` | the trio's SNVs (BCF), the depth bins (BED) and the self-scan of the segment's pieces (PAF) |
| `realign/` | HPRC's sequence over the locus as one FASTA (hap1, hap2, U), and the child's HiFi reads realigned to it |

The coordinates of the event, the assembly's pieces and joins, and every file's own coordinate system are in the
case's README, section 2.

## What differs from the sources

The read records are exactly as published: each BAM holds as many records as its source (checked 2026-11-04). What
changed:
- **Headers.**
  - Every BAM has a new header with no `@PG` lines and no `UR` tags (`samtools reheader -P`), and its records were
    copied as they are.
  - The VCFs lost the header lines that name a file path: `##reference` lines pointing at a producer's file, and
    HPRC's `##bcftools_*Command` lines naming its workflow's files. That is 14 lines in 10 files.
  - The trio-sites BCF keeps its records under a header without the lines naming local files, plus one
    `##assemblyreview_source` line giving the command.
- **The short reads** are the CRAM slices decoded to BAM, so reading them needs no reference.
- **HPRC's QC tracks** are their rows on the child's contigs, as BED. Each file's first line names its source.

No file holds a local path: a scan of every file's raw and decompressed bytes found none (2026-11-04).

## Check them

```bash
cd cases/hg01943/data
python3 -c 'import csv, hashlib; rows = list(csv.DictReader(open("MANIFEST.tsv"), delimiter="\t")); bad = [r["path"] for r in rows if hashlib.md5(open(r["path"], "rb").read()).hexdigest() != r["md5"]]; print(len(rows), "files:", "all match" if not bad else bad)'
```

## Fetch them again

- **Each file:** run its manifest row's command against its source, then clean the header as above.
- **The CRAMs** (the trio's short reads and the mother's nanopore reads) need GRCh38 to be read. Use the analysis set
  from its URL in the manifest, or htslib's lookup of the reference by MD5 at EBI's CRAM reference registry.
- **EBI's range reads** sometimes hang: retry.

Only region slices are kept here because the full public files are large: the GRCh38 analysis set alone is 3.3 GB,
and each genome's read set is many gigabytes.
