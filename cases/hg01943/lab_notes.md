# HG01943: the lab-notebook entry

The case's entry in the workup's lab notebook, as of 2026-11-03, lightly edited for this repo: paths point into
`data/`, and sizes are as of that day. `build_case.py` is the workup's script, which is not public.

## HG01943's duplication against its assembly: which side is wrong (2026-11-03)

* The first card on the follow-up list: HG01943's 1.09-Mb DUP (chr2:46,970,725-48,057,041, standing). The review's
  line: the maternal haplotype runs straight through; the paternal one cannot tell, as it reads parts of the segment twice
  across four gaps. Neither side was taken as the truth: each claim was tested on the raw
  reads. The workup's script, build_case.py, fetched everything by region and built
  hg01943_chr2_dup.html (in `reference/`). Its data is in `data/` (167 MB then),
  numbers in its numbers.json.
* The call holds on every test:
  * Short reads: 6 split reads and 31 pairs in the child (pairs at MAPQ 60; the split reads at 0-60), none in either
    parent. 1KGP's GATK-SV set calls HGSV_24204 independently, child 0/1, parents 0/0, no other carrier of 3,202. The
    junction has 1 bp of microhomology.
  * Not a mapping artifact: 240 of 240 5-kb pieces of the segment and its flanks map first to themselves at MAPQ 60
    (minimap2 against all of GRCh38). The only matches over 1 kb at 90% are interspersed repeats, plus one 3-kb, 95% pair
    inside the segment. Neither breakpoint lies in a segmental duplication.
  * Copies: 3.02 against the parents bin by bin, short reads; HiFi 2.79. Phased HiFi reads give 1.88 maternal copies to
    0.89 paternal (29 bins). At the 221 parent-opposite sites inside, the child's maternal allele is 0.65 short / 0.68
    HiFi; at the 163 outside, 0.50.
  * HiFi junction reads: 31, from two instruments (Sequel II 15, Revio 16) and 5 SMRT cells, both sides at MAPQ 60.
    29 of the 29 phasable carry the mother's alleles. Through the start: 33 maternal and 31 paternal. 31 against 64 and
    69 is constitutional, not mosaic.
  * The mother's nanopore reads: 28 reach a breakpoint and none carries the junction. The father has no long reads.
  * Yale's HPRC R2 call sets of the child: 6 of the 10 read-based ones write the same junction to within 3 bp
    (Sniffles2, debreak, delly, svim, pbsv, sawfish).
* The assembly fails in structure, not in sequence:
  * By the trio's alleles, the paternal scaffold's pieces A, C and E are the father's (A: 352 of his own alleles, none
    of hers). B and D are the mother's: B has 279 of her own alleles and none of his; D has 3 of hers and none of his,
    while E over the same 300 kb has 7 of his.
  * The four joins (the paternal chr2's only gaps, from HPRC's chromosome assignment) are crossed by 0 HiFi reads, with
    79-99 running straight through each point. NucFlag calls each a MISJOIN (and a COLLAPSE in C); Flagger marks only the
    gap bases (NNN).
  * Copies: the assembly holds 3 in 87 of the segment's bins, 2 in 19 and 4 in 3, where the reads hold 3 throughout.
    The second maternal copy is three fragments: B, the unplaced JBHIIL010000075.1 (230 kb) and D. 186,866 bp of it is
    missing: chr2:47,340,687-47,426,808, chr2:47,656,969-47,755,433, ~1.1 kb either side of the junction, and the
    junction itself.
  * The assembly-based call sets inherit the scaffold. PAV writes four paternal insertions (369, 21, 21 and 300 kb,
    1|0); cuteSV-asm four paternal DUPs; SVision-pro-asm breakends at the joins. dipcall and svim-asm write nothing.
* Where the copy was made:
  * The two maternal copies agree at 118 of the mother's het sites before chr2:47,590,351 and at 133 after 47,704,313.
    Between those points they differ at 111 of 117 (the child's share of her own allele 1/3, which needs three copies).
    9 more sites differ in patches in the 95 kb before.
  * The child holds one maternal chromosome. A copy made after conception would repeat it, so the copy was made in the
    mother's germline, from both her chromosomes.
  * Her nanopore reads phase her het sites in 4 blocks, with 0 conflicting links. In the one spanning the stretch
    (47,557,874-47,730,934, 177 sites), both copies carry her haplotype 2 before it, one of each inside it, and both
    haplotype 1 after it (22 of 22).
  * So each copy switches from her haplotype 2 to 1, at a different place. HiFi reads show both switches inside
    molecules: at the right end, 10 go 2 to 1 and 8 stay on 1; at the left, 21 stay on 2 and 21 go 2 to 1.
  * Two exchanges ~120 kb apart in one meiosis are unusual; a template switch while the copy was made may be one of
    them. Which copy is first the reads cannot say.
* For the review, proposed on the page and not applied:
  * rule 3 held: kept;
  * the traversal could read a twice-read haplotype's pieces by the trio's alleles and count the copies against the
    call;
  * an assembly call at a gap counts for nothing;
  * phased long-read depth could join the long-read line for large DUPs;
  * this card could be marked "followed up" with a link, not left on the list.
* The page, tuned the same day for expert review. It opens with a one-picture copy map and five findings. Every number,
  term and figure mark carries a hover note, and a click pins it in a side Inspector (1,908 notes). It adds
  junction-read, join and haplotype-switch figures, six "Show the data" tables, a reviewer's FAQ, and every source,
  threshold and a command to check each claim. The checks added for it:
  * The pieces against each other, every SNV and indel from their CIGARs:
    * B and D are identical to M (0 differences over 368,852 and 300,439 bp);
    * A differs from B at 563 places over 368 kb;
    * E differs from D at 21 places over 299 kb;
    * U differs from M at 130 places, 114 of them in the stretch.
  * The paternal scaffold holds two false duplications where its pieces overlap: A and C over 20,508 bp (6 indels), C
    and E over 17,767 bp (identical). These are the 4-copy bins; Flagger does not flag them.
  * Where the second copy is missing, its reads must land on the best-matching sequence:
    * over chr2:47,340,687-47,426,808 that is M (100 of 100 informative sites);
    * over 47,656,969-47,755,433 it is C (53 of 55, M 28 of 55), where NucFlag calls C collapsed (47,686,528-47,711,595);
    * neither tool flags M.
  * The HiFi phasing checks itself:
    * reads phased to the mother by the parent-opposite sites carry 2,829 mother-only alleles and 10 father-only;
    * reads phased to the father carry 4,336 father-only and 0 mother-only.
  * The only contigs aligned over chr2:46.4-48.65 Mb are CM089354.1, CM089362.1 and JBHIIL010000075.1, so no aligned
    contig holds the junction. The HiFi oversamples the father a little in the left flank (the mother's allele 0.43
    against 0.50 in the short reads).
* The page, reordered later the same day. It now runs short reads (section 1), long reads (2), then the assembly (3).
  A new figure 1 replaces the copy map, which did not read well. It draws the child's two chromosomes as the reads
  hold them, HPRC's two chromosomes and U, and a band from each piece to the copy it holds, numbered 1-5 for what went
  wrong. Figure 2 is new too: the trio's short-read depth with the child's pairs. The claims scorecard moved to section 7.
* The reads given back to the assembly (build_case.py's reads_on_assembly; data in `data/realign/`):
  * HPRC's sequence over the locus comes from the published FASTAs by region read (fetch: asm/pat.region.fa,
    asm/mat.region.fa). It is identical to the bases the Winnowmap records hold, except 8 bases beside gap 3 that no
    record holds.
  * All 6,765 primary HiFi reads of the window were aligned to hap1 46,805,838-49,047,195, hap2
    46,777,285-48,308,368 and U at once (minimap2 -ax map-hifi --eqx -Y, about 10 s).
  * The junction reads all split on the assembly: 29 end at M's segment end and start again at its start (local hap2
    1,308,103 / 220,695). The other 2 put their short side (556 and 1,657 bp) on A or E at MAPQ 0. None runs through.
  * Primary reads per piece, by GRCh38 phase:
    * B: 578 the mother's, 0 the father's, all at MAPQ 0;
    * D: 2 the mother's, 0 the father's, 502 unphased (few informative sites; E is nearly M there);
    * A, C, E: the father's reads only; M and U: the mother's only.
  * One haplotype's depth is 37x (left flank: 41.2 paternal, 33.6 maternal). Depth within 2 kb of each gap's edges
    is 0.01-0.09 of it. The reads straight through each join's points on GRCh38 (the same 80/94, 79/84, 80/81, 99/81)
    align elsewhere unbroken: through A's end, C 23 and M 57; through B's start, A 30 and M 64. None crosses a gap.
  * Caught while building: the gap-edge depth first read 5-kb bins, and a 2-kb window inside one bin came back empty
    (a false 0.00). It is now measured at base resolution.
* Dot plots (exact 16-mer runs, 40 bp or longer drawn, no aligner):
  * Figure 4: the mother's junction read m84081_231124_031433_s4/133041701/ccs (22,246 bp, junction at read base
    10,510) against GRCh38 at the two ends (bases 1-10,510 up to E, 10,513-22,246 from S), and against GRCh38 joined
    as called (one diagonal).
  * Figure 9: the same read against M (the same break), the paternal scaffold (bases 1-9,343 up to D's last base,
    11,625-22,246 from B's first: 2.3 kb across the junction in no piece of that copy), and M joined as called (one
    diagonal).
  * Figure 11, at join 1: a father's read through A's end runs along A to the gap (bases 1-11,310), and its whole
    length lies on C. A mother's read through B's start picks up only at B's first base, and its whole length lies on M.
  * The window hover notes chain runs along the diagonal: right by start, left by end. A looser first version swept in
    stray repeat matches.
