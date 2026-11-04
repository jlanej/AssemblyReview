# AssemblyReview: working rules

What this repo is: [README.md](README.md). The first case and its data: [cases/hg01943/README.md](cases/hg01943/README.md).
What to build: [docs/brief.md](docs/brief.md).

**The stance.** Neither the SV call nor the assembly is the truth. The reviewer tests both on the raw reads, shows
both, and puts their disagreements first.

**A public repo.**
- No local paths, user names or machine names in any committed file. Code takes paths as arguments, or relative to
  the repo.
- Alignment and variant files carry paths in their headers: `@SQ` UR, `@PG` CL, the CRAM file id, and VCF
  `##reference` and command lines. Clean every header before a file is committed:
  - `samtools reheader -P` with the `@PG` lines and `UR` tags removed, and never `samtools reheader -i`;
  - `bcftools reheader` for VCF and BCF.

  Then scan the raw and the decompressed bytes.

**Data.**
- Each case's data lives in `cases/<case>/data/`. Every file is listed in its `MANIFEST.tsv` (path, bytes, md5, what,
  public source, the command that made it). Update the manifest with every change.
- Commit files under 100 MB (GitHub's limit; it warns over 50 MB). A larger file is listed in the manifest with its
  public source and named in `.gitignore`.
- Fetch only from public sources, by region where possible, and list each source in the manifest.
- Read downloaded files with `python3 -I`, from their own folder.

**Git.**
- Commit only when asked, and name explicit paths.
- Push only when asked: a push publishes.
- At the end of a piece of work, propose a commit message and a one-line standup entry.

**Tools.**
- samtools 1.22.1, bcftools 1.22, minimap2 2.30 and Python 3.14 (standard library first).
- Shareable pages are static HTML files in the repo.

Machine-specific rules, if any, live in `CLAUDE.local.md`, which is not tracked.
