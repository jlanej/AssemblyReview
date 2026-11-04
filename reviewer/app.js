/* AssemblyReview: the page. Loads the case's bundle and builds the figures, tables and text from it. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});
  const { h, fmt } = AR;
  const caseId = new URLSearchParams(location.search).get('case') || 'hg01943';

  async function boot() {
    try {
      await AR.loadScript(`cases/${caseId}/data/reviewer/bundle.js`);
      AR.loadScript(`cases/${caseId}/data/reviewer/sequences.js`).then(() => { AR.SEQ = window.ASSEMBLYREVIEW_SEQ; document.dispatchEvent(new Event('ar:sequences')); }).catch(() => {});
    } catch (e) {
      document.getElementById('main').textContent = `The case bundle could not be loaded (cases/${caseId}/data/reviewer/bundle.js). Run tools/prepare_case.py first.`;
      return;
    }
    build(window.ASSEMBLYREVIEW_CASE);
  }

  function build(B) {
    AR.B = B;
    const C = B.case, [S, E] = C.segment, L = E - S + 1, chrom = C.chrom, P = AR.pal;
    const reads = AR.rows(B.reads);
    const chk = new Map(B.checks.map((c) => [c.key, c])); const got = (k) => (chk.get(k) || {}).got;
    const jr = reads.filter((r) => r.junction);
    const n = {
      jr: jr.length, jrMat: AR.count(jr, (r) => r.jr_parent === 'mother'), jrPat: AR.count(jr, (r) => r.jr_parent === 'father'), jrUn: AR.count(jr, (r) => r.jr_parent === 'unphased'),
      thrS: AR.count(reads, (r) => r.through_S), thrE: AR.count(reads, (r) => r.through_E),
      jrBoth: AR.count(jr, (r) => r.asm_split_end && r.asm_split_start), jrThrough: AR.count(jr, (r) => r.asm_through_M),
      cells: new Set(jr.map((r) => r.cell)).size, instruments: new Set(jr.map((r) => r.instrument)).size,
      cn: got('sr_rel') * 2, cnF: got('sr_father') * 2, cnM: got('sr_mother') * 2, cnH: got('hifi_in'), cnHm: got('hifi_mat'), cnHp: got('hifi_pat'),
      afIn: got('oh_in.sr'), afOut: got('oh_out.sr'), afInH: got('oh_in.hifi'), afOutH: got('oh_out.hifi'), nIn: got('oh_in.sr_n'), nOut: got('oh_out.sr_n'),
      holes: B.holes.reduce((s, [a, b]) => s + b - a + 1, 0), differ: got('stretch_differ'), stretchN: got('stretch_n'), before: got('agree_before'), after: got('agree_after'),
      copies: got('copies') || {},
    };
    const perPiece = {}; for (const r of reads) { if (!r.asm_piece) continue; const c = (perPiece[r.asm_piece] = perPiece[r.asm_piece] || { reads: 0, mother: 0, father: 0, unphased: 0, q0: 0 }); c.reads++; c[r.parent]++; if (r.asm_mapq === 0) c.q0++; }
    const vsM = {}; for (const a of B.asm_vs_asm || []) vsM[a.piece] = a;
    const ident = (l, d = 3) => (vsM[l] ? (vsM[l].identity >= 0.99999 ? '100.000%' : fmt.pct(vsM[l].identity, d)) : '–');
    const pB = B.pieces.find((p) => p.piece === 'B'), pD = B.pieces.find((p) => p.piece === 'D');
    const crosses = B.joins.reduce((s, j) => s + j.realign.cross, 0);
    const yale = B.yale_summary, nReadCallers = Object.values(yale).filter((y) => y.made_from === 'HiFi reads').length, nCallIt = Object.values(yale).filter((y) => y.verdict === 'calls it').length;
    const main = document.getElementById('main'); AR.clear(main);
    const sec = (id, ...kids) => { const s = h('section', { id }); for (const k of kids) if (k) s.appendChild(k); main.appendChild(s); return s; };
    const caption = (label, title, text) => h('p', { class: 'caption' }, h('b', {}, `${label} | ${title}`), ' ', text);
    const figbox = (id) => h('div', { class: 'figure', id });
    const stretchLen = B.stretch ? B.stretch[1] - B.stretch[0] + 1 : 0;

    // ---------------------------------------------------------------------------- title and abstract
    document.title = `${C.short} · AssemblyReview`;
    document.getElementById('title').textContent = `A ${fmt.mb(L)} de novo tandem duplication on chromosome 2 in HG01943, and its misassembly in the HPRC release 2 trio-phased assembly`;
    document.getElementById('byline').textContent = `Proband ${C.trio.child}; father ${C.trio.father}; mother ${C.trio.mother} (1000 Genomes). Duplicated segment ${fmt.range(chrom, S, E)} (GRCh38), ${fmt.bp(L)}. Data: 1000 Genomes high-coverage short reads (trio), HPRC HiFi reads (proband), 1000 Genomes ONT reads (mother), HPRC release 2 assembly of the proband (hifiasm, trio mode).`;
    sec('abstract', h('p', { class: 'abstract' },
      `A trio screen of the 1000 Genomes high-coverage genomes reports a de novo tandem duplication of ${fmt.range(chrom, S, E)} in ${C.trio.child}. The HPRC release 2 assembly of the same individual contradicts it: the maternal haplotype (hap2) carries the segment once, and the paternal scaffold (hap1) carries parts of it twice, across four gaps. Both were tested against the raw reads. Read depth gives ${fmt.f(n.cn)} copies in the proband against 2 in each parent; ${n.jr} HiFi molecules span the junction, ${n.jrMat} of the ${n.jrMat + n.jrPat} phasable ones maternal; the maternal allele fraction rises from ${fmt.f(n.afOut)} outside the segment to ${fmt.f(n.afIn)} within it; junction-spanning and continuous reads at the breakpoints are in the ratio of a heterozygous constitutional event. Realigned to the assembly, no read spans any scaffold gap, ${n.jrBoth} of the ${n.jr} junction reads split at both ends of hap2's single copy, and the scaffold segments B and D carry only maternal alleles and maternal reads and are ${ident('B')} and ${ident('D')} identical to hap2: the second maternal copy, broken into three pieces, was placed in the paternal haplotype, with ${fmt.kb(n.holes)} of it, the junction included, in no contig. The two maternal copies differ over ${fmt.kb(stretchLen)}, where one carries the mother's other haplotype, so the duplication arose in the maternal germline from both maternal homologues. The assembly holds the sequence of all three copies in the wrong structure; the call stands.`));

    // ---------------------------------------------------------------------------- Table 1
    const sr = B.sr;
    const t1 = [
      ['Copy number of the segment', `proband ${fmt.f(n.cn)} relative to the parents, bin by bin; father ${fmt.f(n.cnF)}, mother ${fmt.f(n.cnM)}`, `${fmt.f(n.cnH)} overall; by the origin of the reads, maternal ${fmt.f(n.cnHm)}, paternal ${fmt.f(n.cnHp)}`, null, `${Object.entries(n.copies).sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k} copies in ${v} bins`).join(', ')} of 109; the second maternal copy in three contig segments (B, U, D)`],
      ['Junction, segment end to start', `${sr.child.split} split reads and ${sr.child.pairs} discordant pairs in the proband; 0 and 0 in each parent`, `${n.jr} junction-spanning reads, MAPQ 60 on both sides; ${n.instruments} instruments, ${n.cells} SMRT cells`, `0 of ${B.ont.reach} reads reaching a breakpoint; ${B.ont.through_s} and ${B.ont.through_e} continuous across the start and end`, `absent: ${fmt.bp(n.holes)} of the second copy, junction included, in no contig`],
      ['Parental origin of the extra copy', `maternal allele fraction ${fmt.f(n.afIn)} within the segment (${n.nIn} SNVs), ${fmt.f(n.afOut)} outside (${n.nOut})`, `${fmt.f(n.afInH)} within, ${fmt.f(n.afOutH)} outside; ${n.jrMat} of ${n.jrMat + n.jrPat} phasable junction reads maternal, ${n.jrPat} paternal`, null, `hap2 (maternal) one copy; in hap1 (paternal), segments B and D carry ${pB.mom_own} and ${pD.mom_own} maternal-specific alleles and ${pB.dad_own} and ${pD.dad_own} paternal-specific`],
      ['Zygosity', null, `${n.jr} junction-spanning against ${n.thrS} and ${n.thrE} continuous reads at the start and end: heterozygous, constitutional`, null, null],
      ['Reads realigned to the assembly', null, `${crosses} reads span any of the ${B.joins.length} gaps; ${n.jrBoth} of ${n.jr} junction reads split at both ends of hap2's copy, ${n.jrThrough} continuous; B receives ${perPiece.B.mother} maternal reads and ${perPiece.B.father} paternal, all at MAPQ 0`, null, `${B.joins.length} gaps of 100 N in hap1, NucFlag MISJOIN at each; B and D ${ident('B')} and ${ident('D')} identical to hap2 by direct alignment`],
      ['Other call sets', `1KGP GATK-SV: DUP of ${fmt.bp(L)}, PASS; proband 0/1, parents 0/0; the sole carrier among ${fmt.int(B.published[0].columns)} genomes`, `${nCallIt} of ${nReadCallers} HiFi callers report the junction (Sniffles2, pbsv, Sawfish, DeBreak, DELLY, SVIM)`, null, `assembly-based callers report the scaffold: PAV four insertions, cuteSV-asm four duplications, on hap1`],
    ];
    sec('table1', h('h2', {}, 'Table 1 | Evidence by data type'),
      h('div', { class: 'tbl' }, h('table', { class: 'evidence' }, h('thead', {}, h('tr', {}, h('th', {}, ''), h('th', {}, 'Short reads, trio (Illumina, ~35×)'), h('th', {}, 'HiFi, proband (~37× per haplotype)'), h('th', {}, 'ONT, mother (~15×)'), h('th', {}, 'HPRC release 2 assembly, proband'))),
        h('tbody', {}, ...t1.map((r) => h('tr', {}, h('th', {}, r[0]), ...r.slice(1).map((c) => h('td', { class: c == null ? 'na' : null }, c == null ? '—' : c))))))),
      h('p', { class: 'caption' }, 'Copy number from read depth in 10-kb bins with the flanks taken as two copies. Allele fractions are medians over SNVs at which the parents are homozygous for opposite alleles. Every number is recomputed from the files in this repository (Methods 9).'));

    // ---------------------------------------------------------------------------- Figure 1
    const f1 = figbox('fig1');
    const zoom = h('p', { class: 'controls' }, 'View: ');
    f1.appendChild(zoom);
    const fig1 = AR.fig.locus(f1, B, { onLocus: (pos) => gotoGrch38([fmt.locus(chrom, pos - 12000, pos + 12000)]), onPiece: (p) => gotoPiece(p) });
    const zb = (label, a, b) => zoom.appendChild(h('button', { onclick: () => fig1.zoom(a, b) }, label));
    zb('whole window', C.windows.WIN[0], C.windows.WIN[1]); zb('segment ±100 kb', S - 100000, E + 100000); zb('start ±50 kb', S - 50000, S + 50000); zb('end ±50 kb', E - 50000, E + 50000);
    if (B.stretch) zb('discordant interval', B.stretch[0] - 60000, B.stretch[1] + 60000);
    zoom.appendChild(h('span', { class: 'hint' }, ' drag to zoom, double-click to reset, click a position to open it in the interactive view'));
    AR.fig.legend(f1, [['proband', P.ink], ['paternal', P.pat], ['maternal', P.mat], ['junction; copies held by the assembly; absent sequence; discordant SNVs', P.jx], ['scaffold gap', P.gap], ['duplicated segment', P.band], ['discordant interval', P.band2]]);
    sec('figure1', h('h2', {}, 'Figure 1 | Evidence along the locus'), f1,
      caption('Figure 1', 'The duplication in the reads and in the assembly, on one GRCh38 axis.', `(a) Short-read depth of the trio in copies: the proband reads ${fmt.f(n.cn)} copies within the segment (shaded) and 2 outside, each parent 2 throughout; the magenta step is the number of copies of each bin present in the assembly. The arc marks the ${sr.child.pairs} discordant pairs and ${sr.child.split} split reads joining the segment end to its start in the proband; the parents have none. (b) HiFi depth of the proband in all, and by the parental origin of the reads (phased at parent-informative SNVs): about two maternal copies to one paternal. The arc marks the ${n.jr} junction-spanning reads; ${n.thrS} and ${n.thrE} reads are continuous across the start and the end. (c) Fraction of the maternal allele in the proband's short reads at SNVs where the parents are homozygous for opposite alleles: ${fmt.f(n.afIn)} within the segment, ${fmt.f(n.afOut)} outside (dashed line, 2/3). (d) At SNVs heterozygous only in the mother, the fraction of her second allele in the proband's HiFi reads is 0, 1/3 or 2/3 and never 1/2: three copies, the two maternal ones differing (magenta, 1/3) over the discordant interval (purple band). (e) HPRC's alignment of its assembly to GRCh38, drawn as contig position against reference position. hap2 (M) crosses the segment once, continuously. hap1 is five gap-free segments; at each gap (black) the next segment restarts upstream on GRCh38, so B and D, maternal by their alleles, repeat sequence already covered by A and C. U, a maternal contig, is unplaced. (f) The proband's HiFi reads realigned to the assembly: the parental origin of the reads each segment receives. (g) Each hap1 segment and U aligned directly to hap2: B and D are ${ident('B')} and ${ident('D')} identical to hap2; A, C and E differ from it as the paternal haplotype does from the maternal. (h) Structural-variant records over the segment from 22 call sets: the junction (magenta) from read-based callers and the 1KGP set; scaffold artefacts (blue) from assembly-based callers.`));

    // ---------------------------------------------------------------------------- Figure 2
    const f2 = figbox('fig2');
    const seqBox = h('div', { class: 'seqpanel' }); f2.appendChild(seqBox);
    const readsBox = h('div'); f2.appendChild(readsBox);
    AR.fig.junctionReads(readsBox, B, { onRead: (r) => showDotplots(r) });
    const dotBox = h('div', { class: 'dotpanel' }); f2.appendChild(dotBox);
    const dotHead = h('p', { class: 'controls' }, h('b', {}, 'c'), ' One junction-spanning read against the reference and the assembly, as dot plots. Read: ');
    const sel = h('select', { onchange: () => { const r = jr.find((q) => q.name === sel.value); if (r) showDotplots(r, false); } });
    for (const r of jr.slice().sort((a, b) => a.name.localeCompare(b.name))) sel.appendChild(h('option', { value: r.name }, `${r.name} (${r.jr_parent === 'unphased' ? 'unphased' : r.jr_parent === 'mother' ? 'maternal' : 'paternal'}, ${fmt.bp(r.length)})`));
    dotHead.appendChild(sel); dotBox.appendChild(dotHead);
    const dotGrid = h('div', { class: 'dotgrid' }); dotBox.appendChild(dotGrid);
    sec('figure2', h('h2', {}, 'Figure 2 | The breakpoint junction'), f2,
      caption('Figure 2', 'The junction at base resolution, in every spanning molecule, and in one read against reference and assembly.', `(a) GRCh38 at the segment end and start, and the junction as the reads read it: the last base of the segment (${fmt.pos(chrom, E)}) is joined to its first (${fmt.pos(chrom, S)}) with ${B.mh} bp of microhomology (highlighted), too little for recombination between repeats. (b) The ${n.jr} HiFi reads whose alignments end at the segment end and resume at its start on the same strand, each with 100 bp or more clipped at the junction and both alignments at MAPQ 60; ${n.instruments} instruments, ${n.cells} SMRT cells, ${new Set(jr.map((r) => r.side1_start)).size} distinct start positions, so distinct molecules. Ticks mark parent-informative alleles, orange maternal and blue paternal: within 30 kb of the breakpoints ${n.jrMat} reads carry maternal alleles only and ${n.jrPat} paternal; ${n.jrUn} reach no informative site. (c) Exact 16-mer matches of 40 bp or longer, no aligner, between the chosen read (x) and: GRCh38 around the segment start and end, stacked (the diagonal reaches the end and resumes at the start); GRCh38 joined as called (one diagonal); hap2 around the two ends of its single copy (the same break); hap2 joined as called. Click a read in (b) or choose one from the list.`));

    // ---------------------------------------------------------------------------- Figure 3 and Table 2
    const f3 = figbox('fig3');
    AR.fig.structure(f3, B, { onPiece: (p) => gotoPiece(p) });
    f3.appendChild(h('p', { class: 'caption small' }, `1 The two maternal copies collapsed into one, M. 2 Segments of the other maternal copy, B and D, placed in the paternal scaffold; U unplaced. 3 Four joins across 100-N gaps, none supported by a read. 4 ${fmt.kb(n.holes)} of the second copy, and the junction, absent. 5 The paternal segments overlap: A and C over 20.5 kb, C and E over 17.8 kb.`));
    AR.fig.asmDotplot(f3, B);
    AR.fig.assemblyDepth(f3, B, { onLocus: (seq, a, b) => gotoAssembly(`${seq}:${Math.max(1, a - 5000)}-${b + 5000}`) });
    AR.fig.legend(f3, [['maternal origin (alleles, reads)', P.mat], ['paternal origin', P.pat], ['unphased reads', P.grey], ['MAPQ 0 depth', P.ink, 'ln'], ['junction; absent sequence', P.jx], ['gap (100 N)', P.gap], ['NucFlag MISJOIN', '#f59e0b'], ['NucFlag COLLAPSE', '#b71c1c']]);
    sec('figure3', h('h2', {}, 'Figure 3 | The assembly against its own reads and against itself'), f3,
      caption('Figure 3', 'Right sequence, wrong structure.', `(a) The proband's haplotypes as the reads define them (middle rows) and HPRC's assembly of them (top and bottom), to scale; bands join each contig segment to the copy whose sequence it carries. (b) Each gap-free segment of hap1, and U, aligned directly to hap2 with minimap2 (asm20): B and D fall on hap2 at ${ident('B')} and ${ident('D')} identity over ${fmt.kb(vsM.B ? vsM.B.blen : 0)} and ${fmt.kb(vsM.D ? vsM.D.blen : 0)}; A, C and E at ${ident('A', 2)}, ${ident('C', 2)} and ${ident('E', 2)}; U at ${ident('U', 2)}, differing from hap2 mainly within the discordant interval. (c) All ${fmt.int(reads.length)} primary HiFi reads of the window realigned to hap1, hap2 and U together, each read keeping the parental origin its alleles gave it on GRCh38. B and D receive maternal and unphased reads only, every one at MAPQ 0 because an identical placement exists on hap2; depth falls below 0.1 of haploid within 2 kb of every gap and no read spans one; where no contig carries the second copy, hap2 receives two copies' worth of reads (shaded); the ${n.jr} junction-spanning reads split at both ends of hap2's single copy. NucFlag marks every join as a misjoin and part of C as collapsed; Flagger marks only the N bases.`));
    const t2 = B.pieces.map((p) => { const c = perPiece[p.piece] || {}; const ji = B.joins.findIndex((j) => j.left === p.piece); const j = ji >= 0 ? B.joins[ji] : null;
      return { seg: p.piece, contig: `${p.hap} ${p.contig}`, grch38: fmt.range(chrom, p.rs, p.re), len: p.re - p.rs + 1, alleles: `${p.mom_own} / ${p.dad_own}`, origin: p.parent, reads: `${fmt.int(c.mother)} / ${fmt.int(c.father)} / ${fmt.int(c.unphased)}`, q0: c.reads ? fmt.pct(c.q0 / c.reads, 0) : '—', ident: ident(p.piece), join: j ? `gap ${ji + 1} to ${j.right}: ${j.carry} / ${j.realign.cross}` : '—' }; });
    sec('table2', h('h2', {}, 'Table 2 | The contig segments over the duplicated segment'),
      AR.table(t2, [{ key: 'seg', label: 'Segment' }, { key: 'contig', label: 'Haplotype, contig' }, { key: 'grch38', label: 'GRCh38 interval' }, { key: 'len', label: 'Length (bp)', num: true, fmt: (v) => fmt.int(v) }, { key: 'alleles', label: 'Maternal- / paternal-specific alleles' }, { key: 'origin', label: 'Origin by alleles' }, { key: 'reads', label: 'Realigned reads: maternal / paternal / unphased' }, { key: 'q0', label: 'MAPQ 0', num: true }, { key: 'ident', label: 'Identity to hap2', num: true }, { key: 'join', label: 'Join to the next segment: reads supporting it on GRCh38 / spanning the gap on the assembly' }]),
      h('p', { class: 'caption' }, 'Parent-specific alleles are counted at the trio\'s informative SNVs within the segment; a segment is assigned to a parent with 3 or more of that parent\'s alleles and 2% or fewer of the other\'s. Identity from the direct alignment of Figure 3b. A join is supported by a read whose alignment ends within 1 kb of its first point and resumes within 1 kb of the second, both clipped.'));

    // ---------------------------------------------------------------------------- Figure 4
    const f4 = figbox('fig4');
    AR.fig.composition(f4, B);
    AR.fig.legend(f4, [['maternal haplotype 1', P.h1], ['maternal haplotype 2', P.h2], ['the copies differ', P.jx], ['the copies agree', P.mat], ['paternal read', P.pat]]);
    const blk = B.blocks.find((b) => B.stretch && b.start <= B.stretch[0] && b.end >= B.stretch[1]);
    sec('figure4', h('h2', {}, 'Figure 4 | The two maternal copies carry different maternal haplotypes over one interval'), f4,
      caption('Figure 4', 'Where the copy was made.', `(a) At SNVs heterozygous in the mother and homozygous in the father, the fraction of her second allele in the proband's HiFi reads counts the copies carrying it: 0, 1 or 2 of 3. The two maternal copies agree at ${n.before} sites before ${fmt.pos(chrom, B.stretch ? B.stretch[0] : 0)} and at ${n.after} after ${fmt.pos(chrom, B.stretch ? B.stretch[1] : 0)}, and differ at ${n.differ} of ${n.stretchN} between (purple): one copy carries the mother's other haplotype there. A copy made after fertilisation would reproduce the single maternal homologue the proband inherited, so the duplication arose in the maternal germline, from both maternal homologues. (b) The mother's heterozygous sites phased from her ONT reads (${blk ? blk.n : 0} sites in the block spanning the interval, ${B.ont.conflicts} conflicting links): before the interval both copies carry haplotype 2, within it one of each, after it both haplotype 1; each copy switches haplotype once, at a different position. (c, d) HiFi reads spanning each transition, coloured at her phased sites by the haplotype they carry: the switch is seen within single molecules at both ends (${B.switch.left['2>1'] || 0} reads at the proximal transition, ${B.switch.right['2>1'] || 0} at the distal), beside reads that stay on one haplotype and paternal reads. Two exchanges 123 kb apart in one meiosis would be unusual; one may be a template switch during the duplication itself. Which copy is proximal, the reads cannot say.`));

    // ---------------------------------------------------------------------------- interactive views
    const views = sec('views', h('h2', {}, 'Interactive views'), h('p', {}, 'Read-level views (igv.js), loaded when scrolled into view. Reads are coloured by the parental origin of their alleles, orange maternal, blue paternal, grey unphased; junction-spanning reads are outlined in magenta. The views need the page served over HTTP (GitHub Pages or serve.py); the figures do not.'));
    const view = (id, title, btns, creator, note) => {
      const box = h('div', { class: 'view' }); const ctl = h('p', { class: 'controls' }, h('b', {}, title), ' ');
      for (const [label, fn] of btns) ctl.appendChild(h('button', { onclick: fn }, label));
      const holder = h('div', { class: 'igv-holder', id }); holder.appendChild(h('p', { class: 'hint' }, note));
      box.appendChild(ctl); box.appendChild(holder); views.appendChild(box);
      AR.browsers.pending[id.replace('igv-', '')] = () => creator(holder);
      return holder;
    };
    const v1 = view('igv-grch38', 'S1 · GRCh38, the trio at both breakpoints.', [
      ['both breakpoints', () => gotoGrch38()], ['segment start', () => gotoGrch38([fmt.locus(chrom, S - 3000, S + 3000)])], ['segment end', () => gotoGrch38([fmt.locus(chrom, E - 3000, E + 3000)])],
      ...(B.bounds && B.bounds[0] ? [['proximal transition', () => gotoGrch38([fmt.locus(chrom, B.bounds[0] - 6000, B.bounds[1] + 6000)])], ['distal transition', () => gotoGrch38([fmt.locus(chrom, B.bounds[2] - 6000, B.bounds[3] + 6000)])]] : []),
      ...B.joins.map((j, i) => [`gap ${i + 1} points`, () => gotoGrch38([fmt.locus(chrom, j.p1 - 5000, j.p1 + 5000), fmt.locus(chrom, j.p2 - 5000, j.p2 + 5000)])]),
      ['RefSeq genes', async () => { await gotoGrch38(); AR.browsers.addGenes(); }]], (holder) => AR.browsers.grch38(holder, B), 'Tracks: segment, junction-spanning reads, contig segments, absent sequence, NucFlag, informative SNVs, SV records; proband HiFi; trio short reads coloured by pair orientation; mother ONT. Reference: GRCh38 chr2 from this repository.');
    const v2 = view('igv-assembly', 'S2 · The assembly, with the proband HiFi reads realigned to hap1, hap2 and U.', [
      ...B.gaps_local.map((g, i) => [`gap ${i + 1}`, () => gotoAssembly(`${g[0]}:${g[1] - 12000}-${g[2] + 12000}`)]),
      ...(B.m_segment ? [['hap2: copy start', () => gotoAssembly(`hap2:${B.m_segment[0] - 8000}-${B.m_segment[0] + 8000}`)], ['hap2: copy end', () => gotoAssembly(`hap2:${B.m_segment[1] - 8000}-${B.m_segment[1] + 8000}`)]] : []),
      ...Object.entries(B.piece_local).map(([l, [s, a, b]]) => [`segment ${l}`, () => gotoAssembly(`${s}:${Math.max(1, a - 2000)}-${b + 2000}`)])], (holder) => AR.browsers.assembly(holder, B),
      `Coordinates are local to the sequence slice: hap1 = ${C.assembly.sequences.hap1.contig} − ${fmt.int(C.assembly.sequences.hap1.local_offset)}; hap2 = ${C.assembly.sequences.hap2.contig} − ${fmt.int(C.assembly.sequences.hap2.local_offset)}; U whole. MAPQ 0 reads are drawn faint, not hidden.`);
    const v3 = view('igv-junction', 'S3 · The junction as called: GRCh38 segment end (30 kb) joined to segment start (30 kb), with the proband HiFi reads.', [['the join point', () => AR.browsers.goto('junction', 'junction_chr2_48027042-48057041_46970725-47000724:22000-38000')]], (holder) => AR.browsers.junction(holder, B),
      `${got('junction_as_called.jr_unbroken') ?? '?'} of the ${n.jr} junction-spanning reads cross the join point (position 30,000) unbroken; reads continuous in GRCh38 are clipped there.`);
    lazy(v1, () => gotoGrch38()); lazy(v2, () => gotoAssembly()); lazy(v3, () => AR.browsers.goto('junction', 'junction_chr2_48027042-48057041_46970725-47000724:22000-38000'));

    // ---------------------------------------------------------------------------- methods
    const th = C.thresholds;
    const bad = B.checks.filter((c) => !c.ok);
    const cons = B.numbers.consistency || {};
    sec('methods', h('h2', {}, 'Methods'), h('ol', { class: 'methods' },
      h('li', {}, h('b', {}, 'Data. '), `Region slices of public releases: the trio's 1000 Genomes high-coverage CRAMs (${fmt.range(chrom, C.windows.WIN[0], C.windows.WIN[1])}), the proband's HPRC HiFi reads aligned to GRCh38 by Yale (${fmt.range(chrom, C.windows.LWIN[0], C.windows.LWIN[1])}), the mother's 1000 Genomes ONT reads, the HPRC release 2 assembly with its Winnowmap alignments, gaps, NucFlag and Flagger tracks, 16 Yale call sets of the proband and 6 published sets. Every file, its md5, source and command are in the manifest; headers hold no local paths; read records are as published.`),
      h('li', {}, h('b', {}, 'Informative SNVs. '), `bcftools mpileup and call on the trio's short reads; biallelic SNVs with QUAL ≥ ${th.site_qual}, each parent's depth within 0.5–2× its median, a homozygous call with at most 1 read (or 5%) of the other allele, a heterozygous call with ≥ 20% of reads for each allele. Classes: parents homozygous for opposite alleles (${(B.numbers.sites || {}).OH}), mother-only heterozygous (${(B.numbers.sites || {}).MH}), father-only (${(B.numbers.sites || {}).FH}), both heterozygous (${(B.numbers.sites || {}).BH}).`),
      h('li', {}, h('b', {}, 'Parental origin of reads. '), `From the alleles of a read's primary alignment at MAPQ ≥ ${th.mapq}: ${th.phase_opposite_min} or more parent-opposite sites with ${Math.round(th.phase_opposite_agree * 5)} in 5 agreeing; otherwise ${th.phase_own_min} or more of one parent's specific alleles and none of the other's. Junction-spanning reads are phased from both alignments within ${fmt.kb(th.junction_phase_window)} of a breakpoint. The assignment checks itself: reads phased maternal at parent-opposite sites carry ${fmt.int(cons['mat|mother-only alleles'])} maternal-specific and ${cons['mat|father-only alleles']} paternal-specific alleles; reads phased paternal carry ${cons['pat|mother-only alleles']} and ${fmt.int(cons['pat|father-only alleles'])}.`),
      h('li', {}, h('b', {}, 'Junction and continuity. '), `A read spans the junction if one alignment ends within ${th.junction_slack} bp of the segment end with ≥ ${th.junction_clip} bp clipped and another on the same strand begins within ${th.junction_slack} bp of the segment start, clipped; the same rule tests each scaffold join on GRCh38 (within ${fmt.int(th.gap_slack)} bp) and the ends of hap2's copy on the assembly. A read is continuous across a point if one alignment covers ${fmt.int(th.through_flank)} bp on either side. Short reads: split alignments with ≥ 10 bp clipped at each breakpoint; pairs whose forward mate ends within 1 kb before the end and reverse mate starts within 1 kb after the start, pair MAPQ ≥ ${th.mapq}.`),
      h('li', {}, h('b', {}, 'Copy number. '), `samtools bedcov at MAPQ ≥ ${th.mapq} in 10-kb bins; each sample's flank (bins more than 10 kb outside the segment) is taken as 2 copies, and the proband's bins are divided by the parents' mean. HiFi copies by origin split a bin's depth by the depth of reads phased at parent-opposite sites, over bins more than 10 kb inside the segment with ≥ 10 such reads of each parent. At mother-only SNVs with ≥ ${th.site_min_hifi} HiFi reads, the second allele's fraction is read as 0 (< 1/6), 1/3 (≤ 0.42), 1/2 (≤ 0.58, excluded) or 2/3 (≤ 0.83) copies; the discordant interval is the longest run of 1/3 sites, merged across at most two others.`),
      h('li', {}, h('b', {}, 'Realignment to the assembly. '), `All primary HiFi reads of the window aligned with minimap2 (map-hifi, --eqx -Y) to the published hap1 and hap2 sequence over the locus and to U, together; depth in 5-kb bins in multiples of the haploid depth of the left flank (${fmt.f(B.flank_depth.hap1, 1)}× on hap1, ${fmt.f(B.flank_depth.hap2, 1)}× on hap2). Assembly against itself: each gap-free hap1 segment and U aligned to hap2 with minimap2 asm20. Segment differences on GRCh38 from the Winnowmap records' CIGAR and sequence.`),
      h('li', {}, h('b', {}, 'Maternal phase. '), `Each pair of the mother's heterozygous sites on one of her ONT reads (MAPQ ≥ ${th.mapq}) votes; a link needs ≥ ${th.ont_link_min_reads} reads with ≥ ${Math.round(th.ont_link_agree * 100)}% agreement; links are joined strongest first (${B.ont.conflicts} conflicts). Haplotype 1 carries the reference allele at a block's first site.`),
      h('li', {}, h('b', {}, 'Dot plots. '), `Exact ${th.dotplot_k}-mers shared by the read (reference orientation) and each window, both strands, chained along diagonals; runs of ≥ ${th.dotplot_min_run} bp drawn.`),
      h('li', {}, h('b', {}, 'Reproducibility. '), `tools/prepare_case.py recomputes every table on this page from the files in the repository and compares each number with the case record (numbers.json): ${B.checks.length - bad.length} of ${B.checks.length} are equal; the ${bad.length} that differ do so in the third decimal of raw depth or under a definition the original workup did not publish, each listed with its reason in `, h('a', { href: `cases/${C.id}/data/reviewer/checks.tsv` }, 'checks.tsv'), '.')));
    sec('data', h('h2', {}, 'Data and code'), h('p', {}, 'Repository: ', h('a', { href: 'https://github.com/jlanej/AssemblyReview' }, 'github.com/jlanej/AssemblyReview'), `. The page is static (index.html, reviewer/); igv.js 3.8.10 is loaded from jsdelivr for the interactive views. Every figure is drawn from cases/${C.id}/data/reviewer/, rebuilt by `, h('code', {}, `python3 -I tools/prepare_case.py cases/${C.id} --align`), '. Sources and checksums: ', h('a', { href: `cases/${C.id}/data/MANIFEST.tsv` }, 'MANIFEST.tsv'), '. All data are open access (1000 Genomes, HPRC, the 1000 Genomes ONT project).'));

    const defaultRead = () => jr.find((r) => r.name === C.dotplot.read) || jr[0];
    document.addEventListener('ar:sequences', () => { renderSequence(seqBox, B); showDotplots(defaultRead(), false); });
    if (AR.SEQ) { renderSequence(seqBox, B); showDotplots(defaultRead(), false); }
    navSpy();

    // ---- helpers
    function gotoGrch38(loci) { return AR.browsers.goto('grch38', loci || [fmt.locus(chrom, S - 12000, S + 12000), fmt.locus(chrom, E - 12000, E + 12000)]); }
    function gotoAssembly(locus) { return AR.browsers.goto('assembly', locus || (B.gaps_local[0] ? `hap1:${B.gaps_local[0][1] - 12000}-${B.gaps_local[0][2] + 12000}` : 'hap1:1-50000')); }
    function gotoPiece(p) { const [sname, a, b] = B.piece_local[p.piece] || []; if (!sname) return; document.getElementById('igv-assembly').scrollIntoView({ behavior: 'smooth', block: 'center' }); gotoAssembly(`${sname}:${Math.max(1, a - 2000)}-${b + 2000}`); }
    function showDotplots(r, scroll = true) {
      if (!r || !AR.SEQ || !AR.SEQ.reads[r.name]) return;
      sel.value = r.name;
      AR.clear(dotGrid);
      const SEQ = AR.SEQ, read = SEQ.reads[r.name], w = SEQ.windows;
      const sub = (win, a, b) => win.seq.slice(a - win.coord[1], b - win.coord[1] + 1);
      const side1 = r.side1_end - r.side1_start + 1;
      const L1 = Math.min(30000, E - w.grch38_end.coord[1] + 1), L2 = Math.min(30000, w.grch38_start.seq.length - (S - w.grch38_start.coord[1]));
      const plots = [
        { title: 'GRCh38: segment start (lower) and segment end (upper)', windows: [Object.assign({}, w.grch38_start, { marks: [{ at: S - w.grch38_start.coord[1], label: 'segment start' }] }), Object.assign({}, w.grch38_end, { marks: [{ at: E - w.grch38_end.coord[1] + 1, label: 'segment end' }] })] },
        { title: 'GRCh38 joined as called: segment end, then segment start', windows: [{ name: `GRCh38 ${fmt.int(E - L1 + 1)}–${fmt.int(E)} + ${fmt.int(S)}–${fmt.int(S + L2 - 1)}`, coord: ['joined', 1], seq: sub(w.grch38_end, E - L1 + 1, E) + sub(w.grch38_start, S, S + L2 - 1), marks: [{ at: L1, label: 'junction' }] }] },
      ];
      if (w.M_start && w.M_end) { const [ms, me] = SEQ.m_segment; const m1 = Math.min(30000, me - w.M_end.coord[1] + 1), m2 = Math.min(30000, w.M_start.seq.length - (ms - w.M_start.coord[1]));
        plots.push({ title: 'hap2: copy start (lower) and copy end (upper)', windows: [Object.assign({}, w.M_start, { marks: [{ at: ms - w.M_start.coord[1], label: 'copy start' }] }), Object.assign({}, w.M_end, { marks: [{ at: me - w.M_end.coord[1] + 1, label: 'copy end' }] })] });
        plots.push({ title: 'hap2 joined as called', windows: [{ name: `hap2 ${fmt.int(me - m1 + 1)}–${fmt.int(me)} + ${fmt.int(ms)}–${fmt.int(ms + m2 - 1)}`, coord: ['joined', 1], seq: sub(w.M_end, me - m1 + 1, me) + sub(w.M_start, ms, ms + m2 - 1), marks: [{ at: m1, label: 'junction' }] }] }); }
      for (const p of plots) {
        const computed = AR.dotplot.compute(read, p.windows, C.thresholds.dotplot_k, C.thresholds.dotplot_min_run);
        const canvas = h('canvas'); const box = h('div', {}, h('p', { class: 'plottitle' }, p.title), canvas); dotGrid.appendChild(box);
        AR.dotplot.draw(canvas, read, computed, { width: 520, height: 400, xmarks: [{ at: side1, label: 'junction in the read' }], fwdColor: P.pat, revColor: P.mat, markColor: P.jx });
      }
      if (scroll) dotBox.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    }
  }

  function renderSequence(box, B) {
    const SEQ = AR.SEQ; if (!SEQ) return;
    const C = B.case, [S, E] = C.segment, chrom = C.chrom, w = SEQ.windows, mh = B.mh;
    const sub = (win, a, b) => win.seq.slice(a - win.coord[1], b - win.coord[1] + 1);
    const endCtx = sub(w.grch38_end, E - 39, E + 40), startCtx = sub(w.grch38_start, S - 40, S + 39), joined = sub(w.grch38_end, E - 39, E) + sub(w.grch38_start, S, S + 39);
    const mark = (s, from, to) => AR.esc(s.slice(0, from)) + `<span class="mh">${AR.esc(s.slice(from, to))}</span>` + AR.esc(s.slice(to));
    AR.clear(box);
    box.appendChild(h('p', { class: 'controls' }, h('b', {}, 'a'), ' The junction, base by base'));
    const lines = [
      `GRCh38 at the segment end    ${fmt.range(chrom, E - 39, E + 40)}   ` + mark(endCtx, 40 - mh, 40),
      `junction, as the reads read it                             ` + mark(joined, 40 - mh, 40),
      `GRCh38 at the segment start  ${fmt.range(chrom, S - 40, S + 39)}   ` + mark(startCtx, 40 - mh, 40),
    ];
    box.appendChild(h('div', { class: 'seqbox', html: lines.join('\n') }));
    box.appendChild(h('p', { class: 'hint' }, `The junction line is the segment's last 40 bases (to ${fmt.pos(chrom, E)}) followed by its first 40 (from ${fmt.pos(chrom, S)}); the ${mh} base shared by both sides is highlighted.`));
  }

  function lazy(el, fn) { if (!('IntersectionObserver' in window)) return; const io = new IntersectionObserver((es) => { for (const e of es) if (e.isIntersecting) { io.disconnect(); fn(); } }, { rootMargin: '200px' }); io.observe(el); }
  function navSpy() {
    const links = [...document.querySelectorAll('nav a')]; const secs = links.map((a) => document.getElementById(a.getAttribute('href').slice(1))).filter(Boolean);
    const io = new IntersectionObserver((es) => { for (const e of es) if (e.isIntersecting) links.forEach((a) => a.classList.toggle('active', a.getAttribute('href') === '#' + e.target.id)); }, { rootMargin: '-30% 0px -60% 0px' });
    secs.forEach((s) => io.observe(s));
  }
  boot();
})();
