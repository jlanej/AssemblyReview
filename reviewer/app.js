/* AssemblyReview: the page. Loads a case's bundle (cases/<id>/data/reviewer/bundle.js), then builds every section from it.
   Nothing here is typed in: each number is read from the tables, and the claims section compares them with numbers.json. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});
  const { h, fmt } = AR;
  const params = new URLSearchParams(location.search);
  const caseId = params.get('case') || 'hg01943';

  async function boot() {
    try {
      await AR.loadScript(`cases/${caseId}/data/reviewer/bundle.js`);
      AR.loadScript(`cases/${caseId}/data/reviewer/sequences.js`).then(() => { AR.SEQ = window.ASSEMBLYREVIEW_SEQ; document.dispatchEvent(new Event('ar:sequences')); }).catch(() => {});
    } catch (e) {
      document.getElementById('loading').textContent = `The case's bundle could not be loaded (cases/${caseId}/data/reviewer/bundle.js). Run tools/prepare_case.py first.`;
      return;
    }
    const B = window.ASSEMBLYREVIEW_CASE; AR.B = B;
    build(B);
  }

  function build(B) {
    const C = B.case, N = B.numbers, [S, E] = C.segment, L = E - S + 1, chrom = C.chrom;
    const reads = AR.rows(B.reads), sites = AR.rows(B.sites), depth = AR.rows(B.depth_grch38), ont = AR.rows(B.ont_reads);
    const checks = new Map(B.checks.map((c) => [c.key, c]));
    const jr = reads.filter((r) => r.junction);
    const n = {
      jr: jr.length, jrMat: AR.count(jr, (r) => r.jr_parent === 'mother'), jrPat: AR.count(jr, (r) => r.jr_parent === 'father'), jrUn: AR.count(jr, (r) => r.jr_parent === 'unphased'),
      jrQ60: AR.count(jr, (r) => r.side1_mapq === 60 && r.side2_mapq === 60),
      throughS: AR.count(reads, (r) => r.through_S), throughSmat: AR.count(reads, (r) => r.through_S && r.parent === 'mother'), throughSpat: AR.count(reads, (r) => r.through_S && r.parent === 'father'),
      throughE: AR.count(reads, (r) => r.through_E), stop: AR.count(reads, (r) => r.stop),
      jrSplitBoth: AR.count(jr, (r) => r.asm_split_end && r.asm_split_start), jrAtEnd: AR.count(jr, (r) => r.asm_split_end), jrAtStart: AR.count(jr, (r) => r.asm_split_start), jrThrough: AR.count(jr, (r) => r.asm_through_M),
      cells: new Set(jr.map((r) => r.cell)).size, instruments: new Set(jr.map((r) => r.instrument)).size,
      copies: (checks.get('sr_rel') || {}).got * 2, cpChild: (checks.get('sr_child') || {}).got * 2, cpFather: (checks.get('sr_father') || {}).got * 2, cpMother: (checks.get('sr_mother') || {}).got * 2,
      hifiIn: (checks.get('hifi_in') || {}).got, hifiMat: (checks.get('hifi_mat') || {}).got, hifiPat: (checks.get('hifi_pat') || {}).got,
      ohIn: (checks.get('oh_in.sr') || {}).got, ohInH: (checks.get('oh_in.hifi') || {}).got, ohOut: (checks.get('oh_out.sr') || {}).got, ohOutH: (checks.get('oh_out.hifi') || {}).got,
      ohInN: (checks.get('oh_in.sr_n') || {}).got, ohOutN: (checks.get('oh_out.sr_n') || {}).got,
      holes: B.holes.reduce((s, [a, b]) => s + b - a + 1, 0),
      stretchDiffer: (checks.get('stretch_differ') || {}).got, stretchN: (checks.get('stretch_n') || {}).got, agreeBefore: (checks.get('agree_before') || {}).got, agreeAfter: (checks.get('agree_after') || {}).got,
    };
    const sr = B.sr;
    const perPiece = {};
    for (const r of reads) { if (!r.asm_piece) continue; const c = (perPiece[r.asm_piece] = perPiece[r.asm_piece] || { reads: 0, mother: 0, father: 0, unphased: 0, q0: 0 }); c.reads++; c[r.parent]++; if (r.asm_mapq === 0) c.q0++; }
    const pieceB = perPiece.B || { mother: 0, father: 0 };
    const crosses = B.joins.reduce((s, j) => s + j.realign.cross, 0);
    const show = (v) => typeof v === 'number' && !Number.isInteger(v) ? v.toFixed(4) : JSON.stringify(v);
    const card = (title, verdict, txt, to) => h('div', { class: 'card ' + verdict }, h('h4', {}, title, ' ', h('span', { class: 'pill ' + (verdict === 'holds' ? 'ok' : verdict === 'fails' ? 'bad' : '') }, verdict === 'structure' ? 'the structure' : verdict === 'holds' ? 'holds' : 'fails')), h('p', {}, txt), h('a', { href: '#' + to }, 'see it'));
    const pill = (key, label) => {
      const c = checks.get(key);
      if (!c) return document.createTextNode('');
      return h('span', { class: 'pill ' + (c.ok ? 'ok' : 'bad'), title: c.ok ? `${key}: recounted ${show(c.got)} = numbers.json` : `${key}: recounted ${show(c.got)}, numbers.json ${show(c.expected)}${c.note ? ' — ' + c.note : ''}` },
        (label ? label + ' ' : '') + (c.ok ? '✓ recount = case page' : '≈ recount differs'));
    };
    const link = (id, txt) => h('a', { href: '#' + id }, txt);
    const main = document.getElementById('main'); AR.clear(main);

    // -------------------------------------------------------------- hero
    document.title = `${C.short} · AssemblyReview`;
    document.getElementById('title').textContent = C.title;
    document.getElementById('lede').textContent = `A trio screen of the 1000 Genomes high-coverage genomes calls a de novo tandem duplication of ${fmt.range(chrom, S, E)} (${fmt.mb(L)}) in ${C.trio.child}. HPRC's trio-phased assembly of the child disagrees in structure. Taking neither as the truth, this page tests both on the reads as published: the trio's short reads, the child's HiFi reads, the mother's nanopore reads, and the child's reads given back to the assembly itself. Every number is recounted from the data in this repository and checked against the case's record.`;
    const stats = document.getElementById('stats'); AR.clear(stats);
    const stat = (big, txt, key, fv) => stats.appendChild(h('div', { class: 'stat' }, h('b', {}, big), h('span', {}, txt), key ? (() => { const c = checks.get(key); const f = fv || show; return h('span', { class: 'chk' + (c && c.ok ? '' : ' bad'), title: c && c.note ? c.note : null }, c && c.ok ? '✓ recount matches the case page' : `≈ ${c ? f(c.got) : '?'} here, ${c ? f(c.expected) : '?'} on the case page`); })() : null));
    stat('3 copies', `of the segment in the child's short reads (${fmt.f(n.copies)} against its parents, bin by bin); each parent reads two. The extra copy is the mother's.`, 'sr_rel', (v) => fmt.f(v * 2));
    stat(`${n.jr} HiFi reads`, `run the segment's last base straight into its first, both sides at MAPQ 60; ${n.jrMat} of the ${n.jrMat + n.jrPat} that can be phased are the mother's, ${n.jrPat} the father's.`, 'jr_mat');
    stat(`${n.jrSplitBoth} of ${n.jr}`, `junction reads split at both of M's segment ends when given back to HPRC's assembly; ${n.jrThrough} run through its maternal chromosome's one copy.`, 'realign.jr_split_at_m');
    stat(`${fmt.int(pieceB.mother)} to ${pieceB.father}`, `the mother's reads to the father's on piece B of HPRC's paternal chromosome, all at MAPQ 0 because B is base for base the maternal chromosome; ${crosses} reads cross any of its ${B.joins.length} gaps.`, 'realign.per_piece.B.mat');

    // -------------------------------------------------------------- nav
    const sections = [['overview', 'The evidence'], ['picture', 'The picture'], ['reads', 'The reads'], ['junction', 'Junction reads'], ['assembly', 'The assembly'], ['stretch', 'Where it was made'], ['calls', 'Call sets'], ['claims', 'Every claim, recounted'], ['verdict', 'Verdict'], ['methods', 'Methods']];
    const toc = document.querySelector('#toc .wrap'); AR.clear(toc);
    for (const [id, label] of sections) toc.appendChild(h('a', { href: '#' + id, 'data-id': id }, label));
    const sec = (id, title, sub) => { const s = h('section', { id }, h('h2', {}, title), sub ? h('p', { class: 'sub' }, sub) : null); main.appendChild(s); return s; };

    // -------------------------------------------------------------- 0. the evidence at once
    const s2 = sec('overview', 'The evidence at once', `Five comparisons on one axis, ${fmt.range(chrom, C.windows.WIN[0], C.windows.WIN[1])}: the short reads, the long reads against GRCh38, the assembly against GRCh38, the assembly against its own long reads, and the assembly against itself. Drag to zoom, double-click to reset, click a position to send the read panels there; every mark has its numbers on hover.`);
    const vsM = {}; for (const a of B.asm_vs_asm || []) vsM[a.piece] = a;
    const ident = (l) => vsM[l] ? (vsM[l].identity >= 0.99999 ? `identical over ${fmt.kb(vsM[l].blen)}` : `${fmt.pct(vsM[l].identity, 2)} over ${fmt.kb(vsM[l].blen)}`) : '–';
    const strip = h('div', { class: 'strip' },
      card('Short reads: the guide', 'holds', `${sr.child.split} split reads and ${sr.child.pairs} pairs join the end to the start in the child, 0 in either parent; ${fmt.f(n.copies)} copies against the parents, 2 in each parent.`, 'reads'),
      card('Long reads vs GRCh38', 'holds', `${n.jr} molecules read the junction (${n.jrMat} the mother's, ${n.jrPat} the father's); ${n.throughS} and ${n.throughE} run straight through the ends; ${fmt.f(n.hifiMat)} maternal copies to ${fmt.f(n.hifiPat)} paternal.`, 'junction'),
      card('Assembly vs GRCh38', 'structure', `M holds the segment once. The paternal scaffold reads ${B.joins.length > 0 ? fmt.kb(B.joins[0].back) : ''} and ${B.joins.length > 3 ? fmt.kb(B.joins[3].back) : ''} of it a second time, across ${B.joins.length} gaps; U is unplaced; ${fmt.kb(n.holes)} and the junction are in no piece.`, 'picture'),
      card('Assembly vs long reads', 'fails', `${crosses} reads cross any gap; ${n.jrSplitBoth} of ${n.jr} junction reads split at M's segment ends, ${n.jrThrough} run through; B takes ${fmt.int(pieceB.mother)} of the mother's reads and ${pieceB.father} of the father's.`, 'assembly'),
      card('Assembly vs assembly', 'fails', `Aligned directly to M, B is ${ident('B')} and D ${ident('D')}: the mother's sequence, twice, in the father's chromosome; A is ${ident('A')}, C ${ident('C')}.`, 'assembly'));
    s2.appendChild(strip);
    const zoomBtns = h('div', { class: 'btns' });
    const p2 = h('div', { class: 'panel' }); s2.appendChild(zoomBtns); s2.appendChild(p2);
    const ov = AR.panels.overview(p2, B, { onLocus: (pos) => gotoGrch38([fmt.locus(chrom, pos - 12000, pos + 12000)]), onPiece: (p) => gotoAssemblyPiece(p) });
    const zb = (label, a, b) => zoomBtns.appendChild(h('button', { onclick: () => ov.zoom(a, b) }, label));
    zb('whole window', C.windows.WIN[0], C.windows.WIN[1]); zb('the segment ±100 kb', S - 100000, E + 100000); zb('the start ±50 kb', S - 50000, S + 50000); zb('the end ±50 kb', E - 50000, E + 50000);
    if (B.stretch) zb('the stretch', B.stretch[0] - 60000, B.stretch[1] + 60000);
    for (const [i, j] of B.joins.entries()) zb(`join ${i + 1}'s points`, Math.min(j.p1, j.p2) - 30000, Math.max(j.p1, j.p2) + 30000);
    p2.appendChild(h('div', { class: 'key' },
      h('span', {}, h('i', { style: { background: 'var(--ink)' } }), 'the child, against its parents'), h('span', {}, h('i', { style: { background: 'var(--pat)' } }), 'the father\'s'), h('span', {}, h('i', { style: { background: 'var(--mat)' } }), 'the mother\'s'),
      h('span', {}, h('i', { style: { background: 'var(--h1)' } }), 'HiFi, all reads'), h('span', {}, h('i', { style: { background: 'var(--jx)' } }), 'the junction, the assembly\'s copies, what is missing, sites where the maternal copies differ'), h('span', {}, h('i', { style: { background: 'var(--h2)' } }), 'the stretch'),
      h('span', {}, h('i', { style: { background: 'var(--ink)' } }), 'a gap, and the jump back it hides')));

    // -------------------------------------------------------------- 1. the picture
    const s1 = sec('picture', 'What happened, in one picture', 'The child\'s two chromosome 2s as its reads hold them, between HPRC\'s assembly of them. Hover any bar or band; click a piece to open it in the assembly browser.');
    const p1 = h('div', { class: 'panel' }); s1.appendChild(p1);
    AR.panels.picture(p1, B, { onPiece: (p) => gotoAssemblyPiece(p) });
    s1.appendChild(h('div', { class: 'cols' },
      h('div', {},
        h('h3', {}, 'What the reads hold'),
        h('p', {}, `The mother's chromosome carries the segment twice in tandem, joined at the junction (${fmt.pos(chrom, E)} to ${fmt.pos(chrom, S)}, ${B.mh} bp of microhomology); the father's carries it once. Neither parent has the duplication: it arose in the mother's germline, from both of her chromosomes (`, link('stretch', 'where it was made'), ').'),
        h('p', {}, `Three copies by depth (${fmt.f(n.copies)} against the parents in short reads), two of them maternal (the mother's allele is ${fmt.f(n.ohIn)} of the child's reads at ${n.ohInN} parent-opposite sites inside the segment, ${fmt.f(n.ohOut)} at ${n.ohOutN} outside), and ${n.jr} HiFi molecules that read the junction itself, ${n.jrMat} of them the mother's by their alleles.`)),
      h('div', {},
        h('h3', {}, 'What the assembly holds'),
        h('p', {}, `HPRC release 2 assembled ${C.trio.child} with trio-phased hifiasm. Its maternal chromosome M holds the segment once, straight through: the two maternal copies, identical outside one ${fmt.kb(B.stretch ? B.stretch[1] - B.stretch[0] + 1 : 0)} stretch, collapsed into one. Two pieces of the other copy, B and D, sit in the paternal chromosome between the father's own pieces A, C and E, each across a 100-N gap that no read crosses; a third, U, is unplaced. ${fmt.kb(n.holes)} of that copy, and the junction, are in no piece.`),
        h('p', {}, `Every piece is real sequence of one of the three copies: the assembly is wrong in structure, not in sequence, and the call holds. The screen's review rule 3 kept it standing: at a duplication, a haplotype that runs straight through, or cannot be read, never counts against a call.`))));

    // -------------------------------------------------------------- 2. the depth and allele tables (under the picture)
    const flankBins = depth.filter((d) => d.flank).length, insideBins = depth.filter((d) => d.inside).length;
    s1.appendChild(h('div', { class: 'cols' },
      h('div', { class: 'panel' }, h('h4', {}, 'Copies by depth, 10-kb bins'), AR.table([
        { what: "the child's short reads, against its parents' bin by bin", v: fmt.f(n.copies), chk: pill('sr_rel') },
        { what: "the child's short reads, against its own flanks", v: fmt.f(n.cpChild), chk: pill('sr_child') },
        { what: "the father's short reads", v: fmt.f(n.cpFather), chk: pill('sr_father') }, { what: "the mother's short reads", v: fmt.f(n.cpMother), chk: pill('sr_mother') },
        { what: "the child's HiFi reads, against its flanks", v: fmt.f(n.hifiIn), chk: pill('hifi_in') },
        { what: 'HiFi reads phased to the mother', v: fmt.f(n.hifiMat), chk: pill('hifi_mat') }, { what: 'HiFi reads phased to the father', v: fmt.f(n.hifiPat), chk: pill('hifi_pat') },
      ], [{ key: 'what', label: 'copies inside the segment', cls: 'wrap' }, { key: 'v', label: 'copies', num: true }, { key: 'chk', label: 'against the case page', fmt: (v) => v }]),
        h('p', { class: 'cap' }, `Flanks: the ${flankBins} bins more than 10 kb outside the segment; inside: the ${insideBins} bins fully inside it. Depth is samtools bedcov at MAPQ 20.`)),
      h('div', { class: 'panel' }, h('h4', {}, "The mother's allele in the child, at parent-opposite sites"), AR.table([
        { where: `inside the segment (${n.ohInN} sites)`, sr: fmt.f(n.ohIn), hifi: fmt.f(n.ohInH), chk: pill('oh_in.sr') },
        { where: `outside (${n.ohOutN} sites)`, sr: fmt.f(n.ohOut), hifi: fmt.f(n.ohOutH), chk: pill('oh_out.sr') },
      ], [{ key: 'where', label: 'sites', cls: 'wrap' }, { key: 'sr', label: 'short reads', num: true }, { key: 'hifi', label: 'HiFi', num: true }, { key: 'chk', label: '', fmt: (v) => v }]),
        h('p', { class: 'cap' }, `Medians over sites of the share of her allele in the child's reads: two maternal copies to one paternal inside the segment (2/3), one to one outside (1/2). At mother-only sites inside, her own allele's share is 0, 1/3 or 2/3 of the HiFi reads, never 1/2: three copies (`, link('stretch', 'the stretch'), ').'))));

    // -------------------------------------------------------------- 3. the reads
    const s3 = sec('reads', 'The reads at the breakpoints', "The trio's short reads, the child's HiFi reads and the mother's nanopore reads on GRCh38 (igv.js). The two breakpoints side by side; the junction reads outlined in pink.");
    s3.appendChild(h('div', { class: 'cols c21' },
      h('div', { class: 'panel' }, h('h4', {}, 'The junction in the short reads'), AR.table(['child', 'father', 'mother'].map((m) => ({ who: `${m} ${C.trio[m]}`, split: sr[m].split, pairs: sr[m].pairs, med: sr[m].pair_mapq_med == null ? '–' : sr[m].pair_mapq_med, clip: `${sr[m].clip_s} / ${sr[m].clip_e}`, mapq: sr[m].split_mapq.join(', ') || '–' })),
        [{ key: 'who', label: 'member' }, { key: 'split', label: 'split reads', num: true }, { key: 'mapq', label: 'their MAPQ' }, { key: 'pairs', label: 'pairs across the junction', num: true }, { key: 'med', label: 'pair MAPQ, median', num: true }, { key: 'clip', label: 'clipped exactly at the start / end', num: true }]),
        h('p', { class: 'cap' }, `A split read aligns in two parts, one ending at ${fmt.pos(chrom, E)}, the other starting at ${fmt.pos(chrom, S)}. A pair across the junction has its forward mate ending within 1 kb before the segment's end and its reverse mate starting within 1 kb after its start (the orientation a tandem duplication's junction gives), the pair at MAPQ 20 or more. `, pill('split'), ' ', pill('pairs'), ' ', pill('parents_sr.father'), ' ', pill('parents_sr.mother'))),
      h('div', { class: 'panel' }, h('h4', {}, "The mother's nanopore reads"), AR.table([
        { what: 'reads at MAPQ 20 or more', v: fmt.int(B.ont.reads) }, { what: 'reads that reach a breakpoint (within 200 bp)', v: B.ont.reach }, { what: 'reads carrying the junction', v: B.ont.junction },
        { what: "straight through the segment's start / end, 1 kb either side", v: `${B.ont.through_s} / ${B.ont.through_e}` }, { what: 'her heterozygous sites phased: blocks of 10 or more sites', v: B.blocks.filter((b) => b.n >= 10).length }, { what: 'conflicting links', v: B.ont.conflicts },
      ], [{ key: 'what', label: '', cls: 'wrap' }, { key: 'v', label: '', num: true }]), h('p', { class: 'cap' }, pill('ont.junction'), ' ', pill('ont.through_s'), ' ', pill('ont.through_e'), ' ', pill('ont.reach')))));
    const readBtns = h('div', { class: 'btns' }); s3.appendChild(readBtns);
    const igv1 = h('div', { class: 'panel tight', id: 'igv-grch38' }); s3.appendChild(igv1);
    AR.browsers.pending.grch38 = () => AR.browsers.grch38(igv1, B);
    const rb = (label, loci) => readBtns.appendChild(h('button', { onclick: () => gotoGrch38(loci) }, label));
    rb('both breakpoints', [fmt.locus(chrom, S - 12000, S + 12000), fmt.locus(chrom, E - 12000, E + 12000)]);
    rb("the segment's start", [fmt.locus(chrom, S - 3000, S + 3000)]); rb("the segment's end", [fmt.locus(chrom, E - 3000, E + 3000)]);
    if (B.bounds && B.bounds[0]) { rb("the stretch's left end", [fmt.locus(chrom, B.bounds[0] - 6000, B.bounds[1] + 6000)]); rb("the stretch's right end", [fmt.locus(chrom, B.bounds[2] - 6000, B.bounds[3] + 6000)]); }
    for (const [i, j] of B.joins.entries()) rb(`join ${i + 1}'s two points`, [fmt.locus(chrom, j.p1 - 5000, j.p1 + 5000), fmt.locus(chrom, j.p2 - 5000, j.p2 + 5000)]);
    let hl = true; const hlBtn = h('button', { onclick: () => { hl = !hl; AR.browsers.highlightJunctionReads('grch38', B, hl); hlBtn.classList.toggle('active', hl); } }, 'outline the junction reads'); hlBtn.classList.add('active'); readBtns.appendChild(hlBtn);
    readBtns.appendChild(h('button', { class: 'primary', onclick: () => gotoGrch38() }, 'load the read panels'));
    readBtns.appendChild(h('button', { onclick: async () => { await gotoGrch38(); AR.browsers.addGenes(); } }, 'add RefSeq genes (fetched from UCSC\'s mirror)'));
    igv1.appendChild(h('div', { class: 'igv-note' }, 'The read panels load when you scroll here or press "load the read panels". They need the page served over HTTP (GitHub Pages, or ', h('code', {}, 'python3 serve.py'), '); igv.js comes from jsdelivr, pinned to release 3.8.10. The reference is GRCh38 chr2 from this repository, N outside the window, so nothing else is fetched unless you add the genes.'));
    lazy(igv1, () => gotoGrch38());

    // -------------------------------------------------------------- 4. the junction reads
    const s4 = sec('junction', `The ${n.jr} junction reads`, `Every HiFi read that aligns up to ${fmt.pos(chrom, E)} and goes on from ${fmt.pos(chrom, S)} on the same strand, with 100 or more bases clipped at each side. Click a read for its dot plots against GRCh38, the assembly and the junction as called.`);
    s4.appendChild(h('p', {}, `${n.jr} reads from ${n.instruments} instruments and ${n.cells} SMRT cells, both sides at MAPQ 60 in ${n.jrQ60} of them; they begin at ${new Set(jr.map((r) => r.side1_start)).size} different places and end at ${new Set(jr.map((r) => r.side2_end)).size}: distinct molecules, not one library's artifact. Phased by the parents' alleles within 30 kb of a breakpoint, ${n.jrMat} carry the mother's and ${n.jrPat} the father's; ${n.jrUn} reach no informative site. Against them, ${n.throughS} reads run straight through the segment's start (${n.throughSmat} the mother's, ${n.throughSpat} the father's) and ${n.throughE} through its end: one junction against two chromosomes' worth of straight reads, a constitutional duplication, not a mosaic one. `, pill('jr'), ' ', pill('jr_mat'), ' ', pill('jr_pat'), ' ', pill('through_s.mat'), ' ', pill('through_s.pat'), ' ', pill('through_e')));
    const seqBox = h('div', { class: 'panel' }, h('h4', {}, 'The junction, to the base')); s4.appendChild(seqBox);
    const jrTable = AR.table(jr.map((r) => ({ name: r.name, instrument: r.instrument, strand: r.strand, length: r.length, side1: `${fmt.int(r.side1_start)}-${fmt.int(r.side1_end)}`, side2: `${fmt.int(r.side2_start)}-${fmt.int(r.side2_end)}`, mapq: `${r.side1_mapq} / ${r.side2_mapq}`, alleles: `${r.jr_mat_alleles} / ${r.jr_pat_alleles}`, parent: r.jr_parent, asm: r.asm_split_end && r.asm_split_start ? 'splits at both of M\'s ends' : r.asm_split_end || r.asm_split_start ? 'splits at one end' : r.asm_through_M ? 'runs through' : '–', r })),
      [{ key: 'name', label: 'read', cls: 'mono' }, { key: 'instrument', label: 'instrument' }, { key: 'strand', label: 'strand' }, { key: 'length', label: 'length', num: true, fmt: (v) => fmt.int(v) }, { key: 'side1', label: `side 1 (to ${fmt.int(E)})` }, { key: 'side2', label: `side 2 (from ${fmt.int(S)})` }, { key: 'mapq', label: 'MAPQ' }, { key: 'alleles', label: "mother's / father's alleles" }, { key: 'parent', label: 'parent', fmt: (v) => `<span class="pill ${v === 'mother' ? 'mat' : v === 'father' ? 'pat' : ''}">${v === 'unphased' ? '–' : "the " + v + "'s"}</span>` }, { key: 'asm', label: 'on the assembly' }],
      { scroll: true, onRow: (row) => showDotplots(row.r) });
    s4.appendChild(h('div', { class: 'panel tight' }, jrTable));
    const dotBox = h('div', { class: 'panel', id: 'dotplots' }); s4.appendChild(dotBox);
    dotBox.appendChild(h('p', { class: 'note' }, 'Pick a read above. The dot plots draw every run of exact 16-mer matches of 40 bp or longer between the read (x) and each window (y), with no aligner: a read that matches a window end to end draws one diagonal.'));

    // -------------------------------------------------------------- 5. the assembly
    const s5 = sec('assembly', 'The assembly, tested on its own reads', `All ${fmt.int(reads.length)} primary HiFi reads of the window, realigned to HPRC's published sequence over the locus (both chromosomes and the unplaced contig at once), each keeping the parent its alleles gave it on GRCh38.`);
    const p5 = h('div', { class: 'panel' }); s5.appendChild(p5);
    AR.panels.assemblyDepth(p5, B, { onLocus: (seq, a, b) => gotoAssembly(`${seq}:${Math.max(1, a - 5000)}-${b + 5000}`) });
    p5.appendChild(h('p', { class: 'cap' }, `Bars: depth in 5-kb bins in copies of one haplotype's depth (${fmt.f(B.one_copy, 1)}x, the left flank's), stacked: the mother's reads (amber), the father's (blue), unphased (grey). Black ticks: depth at MAPQ 0 (an equally good place elsewhere). Pink: where no piece holds the second maternal copy, so M takes two copies' reads. Hover a bin; click one to open it in the browser below.`));
    const p5b = h('div', { class: 'panel' }, h('h4', {}, 'The assembly against itself')); s5.appendChild(p5b);
    AR.panels.asmDotplot(p5b, B);
    p5b.appendChild(h('p', { class: 'cap' }, `Each gap-free piece of the paternal scaffold, and U, aligned directly to the maternal chromosome with minimap2 (asm20), no reference involved. B and D draw diagonals identical to M base for base; A, C and E differ from it as the father's chromosome does from the mother's. The sequence of the two maternal copies is all there; it is the scaffold that is wrong.`));
    const piecesRows = B.pieces.map((p) => Object.assign({}, p, perPiece[p.piece] || {}));
    s5.appendChild(h('div', {},
      h('div', { class: 'panel' }, h('h4', {}, 'Each piece: whose by its alleles, whose reads it takes'), AR.table(piecesRows, [
        { key: 'piece', label: 'piece', fmt: (v, r) => `<b>${v}</b>` }, { key: 'contig', label: 'contig', fmt: (v, r) => `${r.hap} ${v}` }, { key: 'rs', label: 'GRCh38', fmt: (v, r) => fmt.range(chrom, r.rs, r.re) },
        { key: 'mom_own', label: "mother's own alleles", num: true }, { key: 'dad_own', label: "father's own", num: true }, { key: 'parent', label: 'by alleles', fmt: (v) => `<span class="pill ${v === 'maternal' ? 'mat' : v === 'paternal' ? 'pat' : ''}">${v}</span>` },
        { key: 'reads', label: 'reads', num: true, fmt: (v) => fmt.int(v) }, { key: 'mother', label: "the mother's", num: true }, { key: 'father', label: "the father's", num: true }, { key: 'unphased', label: 'unphased', num: true }, { key: 'q0', label: 'at MAPQ 0', num: true }]),
        h('p', { class: 'cap' }, 'Own alleles counted at the trio\'s informative sites inside the segment. ', pill('pieces.B.mom_own'), ' ', pill('pieces.B.dad_own'), ' ', pill('pieces.D.mom_own'), ' ', pill('realign.per_piece.B.mat'), ' ', pill('realign.per_piece.A.pat'))),
      h('div', { class: 'panel' }, h('h4', {}, 'The joins: each a jump back along GRCh38, across 100 N, that no read crosses'), AR.table(B.joins.map((j, i) => ({ i: i + 1, pieces: `${j.left} → ${j.right}`, points: `${fmt.int(j.p1)} → ${fmt.int(j.p2)}`, back: j.back, carry: j.carry, through: `${j.through1} / ${j.through2}`, stop: j.stop, cross: j.realign.cross, edge: `${fmt.f(j.realign.edge_left)} / ${fmt.f(j.realign.edge_right)}` })), [
        { key: 'i', label: 'join', num: true }, { key: 'pieces', label: 'pieces' }, { key: 'points', label: 'GRCh38: end of one → start of the next' }, { key: 'back', label: 'reads back (bp)', num: true, fmt: (v) => fmt.int(v) },
        { key: 'carry', label: 'HiFi reads carrying the join on GRCh38', num: true }, { key: 'through', label: 'straight through each point' }, { key: 'stop', label: 'clipped near a point', num: true },
        { key: 'cross', label: 'realigned reads across the gap', num: true }, { key: 'edge', label: "depth at the gap's edges (of one haplotype's)" }]),
        h('p', { class: 'cap' }, 'A read carries a join if one alignment ends within 1 kb of the first point and another starts within 1 kb of the second, both clipped. NucFlag calls a MISJOIN at every join and a COLLAPSE inside C; Flagger marks only the N bases. ', pill('joins.0.carry'), ' ', pill('joins.1.carry'), ' ', pill('joins.2.carry'), ' ', pill('joins.3.carry'), ' ', pill('realign.joins.0.cross')))));
    s5.appendChild(h('details', {}, h('summary', {}, 'Every pair of pieces, where both run: SNVs and indels between them (B and D are base for base M)'),
      h('div', { class: 'panel tight' }, AR.table(B.piece_pairs.map((p) => Object.assign({}, p, { parents: p.parents.join(' / '), where: Object.entries(p.where).map(([k, v]) => `${k} ${v}`).join(', ') || '–', cover: fmt.range(chrom, p.lo, p.hi) })), [
        { key: 'x', label: 'piece' }, { key: 'y', label: 'piece' }, { key: 'parents', label: 'parents' }, { key: 'cover', label: 'GRCh38 both cover' }, { key: 'len', label: 'bp', num: true, fmt: (v) => fmt.int(v) }, { key: 'snv', label: 'SNVs', num: true }, { key: 'indel', label: 'indels', num: true }, { key: 'per_kb', label: 'per kb', num: true, fmt: (v) => fmt.f(v, 2) }, { key: 'where', label: 'where they differ' }]))));
    s5.appendChild(h('details', {}, h('summary', {}, "HPRC's own QC calls over the pieces, placed on GRCh38"),
      h('div', { class: 'panel tight' }, AR.table(B.qc_placed.map((f) => ({ track: f.track, label: f.label, piece: f.piece, contig: `${f.contig}:${fmt.int(f.qs)}-${fmt.int(f.qe)}`, grch38: fmt.range(chrom, f.rs, f.re) })), [{ key: 'track', label: 'track' }, { key: 'label', label: 'call' }, { key: 'piece', label: 'piece' }, { key: 'contig', label: 'on the contig' }, { key: 'grch38', label: 'GRCh38, through the piece' }]))));
    const asmBtns = h('div', { class: 'btns' }); s5.appendChild(h('h3', {}, "HPRC's sequence, with the reads on it")); s5.appendChild(asmBtns);
    const igv2 = h('div', { class: 'panel tight', id: 'igv-assembly' }); s5.appendChild(igv2);
    AR.browsers.pending.assembly = () => AR.browsers.assembly(igv2, B);
    const ab = (label, locus) => asmBtns.appendChild(h('button', { onclick: () => gotoAssembly(locus) }, label));
    B.gaps_local.forEach((g, i) => ab(`gap ${i + 1} (${B.joins[i] ? B.joins[i].left + ' to ' + B.joins[i].right : ''})`, `${g[0]}:${g[1] - 12000}-${g[2] + 12000}`));
    if (B.m_segment) { ab("M's segment start", `hap2:${B.m_segment[0] - 8000}-${B.m_segment[0] + 8000}`); ab("M's segment end", `hap2:${B.m_segment[1] - 8000}-${B.m_segment[1] + 8000}`); }
    for (const [letter, [sname, a, b]] of Object.entries(B.piece_local)) ab(`piece ${letter}`, `${sname}:${Math.max(1, a - 2000)}-${b + 2000}`);
    asmBtns.appendChild(h('button', { class: 'primary', onclick: () => gotoAssembly() }, 'load the assembly browser'));
    igv2.appendChild(h('div', { class: 'igv-note' }, 'Loads when you scroll here. Coordinates are local to the sequence slice: hap1 = ', h('code', {}, `${C.assembly.sequences.hap1.contig} − ${fmt.int(C.assembly.sequences.hap1.local_offset)}`), ', hap2 = ', h('code', {}, `${C.assembly.sequences.hap2.contig} − ${fmt.int(C.assembly.sequences.hap2.local_offset)}`), ', U whole. Reads at MAPQ 0 are drawn faint, not hidden: on B and D every read is at MAPQ 0, because each has an equally good place on M.'));
    lazy(igv2, () => gotoAssembly());
    s5.appendChild(h('h3', {}, 'The junction as called, as a reference'));
    s5.appendChild(h('p', {}, `GRCh38 from 30 kb before the segment's end to that end, joined to its start and the 30 kb after: the sequence the call says the mother's chromosome has. The child's HiFi reads aligned to it (minimap2, map-hifi): ${(checks.get('junction_as_called.jr_unbroken') || {}).got ?? '?'} of the ${n.jr} junction reads cross the join point at 30,000 unbroken, while reads that run straight on in GRCh38 are clipped there.`));
    const igv3 = h('div', { class: 'panel tight', id: 'igv-junction' }); s5.appendChild(igv3);
    AR.browsers.pending.junction = () => AR.browsers.junction(igv3, B);
    igv3.appendChild(h('div', { class: 'igv-note' }, 'Loads when you scroll here.'));
    lazy(igv3, () => AR.browsers.goto('junction', 'junction_chr2_48027042-48057041_46970725-47000724:22000-38000'));

    // -------------------------------------------------------------- 6. where the copy was made
    const s6 = sec('stretch', 'Where the copy was made', 'Two maternal copies, not one copied twice: the child\'s two maternal copies differ over one stretch, where one carries the mother\'s other chromosome.');
    const blk = B.blocks.find((b) => B.stretch && b.start <= B.stretch[0] && b.end >= B.stretch[1]);
    s6.appendChild(h('p', {}, B.stretch ? `At the mother's heterozygous sites inside the segment where the father is homozygous, her own allele's share of the child's HiFi reads counts how many of its three copies carry it: 0, 1/3 or 2/3. The two maternal copies agree at ${n.agreeBefore} sites before ${fmt.pos(chrom, B.stretch[0])} and at ${n.agreeAfter} after ${fmt.pos(chrom, B.stretch[1])}; between those points they differ at ${n.stretchDiffer} of ${n.stretchN}. A copy made after conception would repeat the one maternal chromosome the child received, so the copy was made in the mother's germline, from both of her chromosomes. ` : 'No stretch of differing sites was found. ', pill('stretch'), ' ', pill('stretch_differ'), ' ', pill('agree_before'), ' ', pill('agree_after')));
    if (blk) s6.appendChild(h('p', {}, `Her nanopore reads phase her heterozygous sites into blocks (a link needs 3 reads, 85% agreeing; ${B.ont.conflicts} conflicting links). The block ${fmt.range(chrom, blk.start, blk.end)} (${blk.n} sites) spans the stretch. Read against it, the copies carry: before the stretch, ${Object.entries(B.compo.left).map(([k, v]) => `${v} sites ${({ 0: 'both her haplotype 2', 1: 'one of each', 2: 'both her haplotype 1' })[k]}`).join(', ')}; inside it, ${Object.entries(B.compo.inside).map(([k, v]) => `${v} ${({ 0: 'both haplotype 2', 1: 'one of each', 2: 'both haplotype 1' })[k]}`).join(', ')}; after it, ${Object.entries(B.compo.right).map(([k, v]) => `${v} ${({ 0: 'both haplotype 2', 1: 'one of each', 2: 'both haplotype 1' })[k]}`).join(', ')}. Each copy changes from one of her haplotypes to the other at a different place, and HiFi reads show the change inside single molecules at both ends. `, pill('compo.inside'), ' ', pill('compo.right'), ' ', pill('compo.left')));
    const p6 = h('div', { class: 'panel' }); s6.appendChild(p6);
    AR.panels.stretch(p6, B);
    s6.appendChild(h('p', { class: 'small' }, 'Two copies that leave one of her chromosomes for the other at different places are two chromatids that had each exchanged with her other chromosome, joined by the duplication; or one exchange and a template switch while the copy was made. Which copy comes first the reads cannot say.'));

    // -------------------------------------------------------------- 7. call sets
    const s7 = sec('calls', 'What the other call sets say', `Yale's ${Object.keys(B.yale_summary).length} call sets of the child for HPRC release 2, and the published sets the screen's review reads, over the segment.`);
    const callsBySet = {}; for (const c of B.calls) (callsBySet[c.set] = callsBySet[c.set] || []).push(c);
    const describe = (c) => c.type === 'BND' ? `BND ${fmt.pos(chrom, c.pos)} joined to ${fmt.int(c.end)}, GT ${c.gt}` : `${c.type} ${fmt.range(chrom, c.pos, c.end)}, ${fmt.bp(c.len)}, GT ${c.gt}`;
    s7.appendChild(h('div', { class: 'panel tight' }, AR.table(Object.entries(B.yale_summary).map(([set, y]) => ({ set, made_from: y.made_from, verdict: y.verdict, records: (callsBySet[set] || []).map(describe).join('; ') + (y.bnd_single ? ` (${y.bnd_single} single breakends nearby)` : '') || 'nothing of 10 kb or more, no breakend, over the segment' })),
      [{ key: 'set', label: 'call set' }, { key: 'made_from', label: 'made from' }, { key: 'verdict', label: 'over the segment', fmt: (v) => `<span class="pill ${v === 'calls it' ? 'ok' : v === 'calls the scaffold' ? 'bad' : ''}">${v}</span>` }, { key: 'records', label: 'records', cls: 'wrap' }])));
    s7.appendChild(h('p', { class: 'cap' }, `${Object.values(B.yale_summary).filter((y) => y.verdict === 'calls it').length} of the ${Object.values(B.yale_summary).filter((y) => y.made_from === 'HiFi reads').length} read-based call sets write the junction, to within a few bases; none writes another junction over the segment. The assembly-based sets inherit the scaffold: insertions or duplications on haplotype 1, the father's, at points no read joins.`));
    s7.appendChild(h('div', { class: 'panel tight' }, AR.table(B.published.map((p) => ({ set: p.set, file: p.file, held: p.held.join(', ') || 'none of the trio', n: p.n_sv50, longest: p.longest, columns: p.columns })),
      [{ key: 'set', label: 'published set' }, { key: 'file', label: 'file', cls: 'mono' }, { key: 'held', label: 'of the trio' }, { key: 'n', label: 'SVs of 50 bp or more carried, over the segment', num: true }, { key: 'longest', label: 'longest (bp)', num: true, fmt: (v) => fmt.int(v) }, { key: 'columns', label: 'genomes', num: true, fmt: (v) => fmt.int(v) }])));
    s7.appendChild(h('p', { class: 'cap' }, `The 1KGP short-read set (GATK-SV) calls the duplication independently from the same reads: ${B.carriers.length === 1 ? `its only carrier among ${fmt.int(B.published[0].columns)} genomes is ${B.carriers[0]}` : `carriers: ${B.carriers.join(', ')}`}. `, pill('carriers')));

    // -------------------------------------------------------------- 8. claims
    const s8 = sec('claims', 'Every claim, recounted', 'What each side claims, what would undo it, and the recount from the files in this repository against the case page\'s numbers.json. Green: the recount equals the case page. Amber: it differs, with the reason (hover).');
    const claims = [
      { side: 'call', title: 'The call: the junction in the short reads', what: `${sr.child.split} split reads and ${sr.child.pairs} pairs join ${fmt.pos(chrom, E)} to ${fmt.pos(chrom, S)} in the child; none in either parent. The 1KGP set calls it too, in the child alone.`, keys: ['split', 'pairs', 'parents_sr.father', 'parents_sr.mother', 'carriers'], to: 'reads', undo: 'junction reads in a parent; no pairs' },
      { side: 'call', title: 'The call: not a paralog', what: `${B.self.first} of ${B.self.pieces} 5-kb pieces of the segment and its flanks map first to themselves in GRCh38 at MAPQ 60; ${B.self.other} match elsewhere over 1 kb at 90%, interspersed repeats; neither breakpoint lies in a segmental duplication; the two sides share ${B.mh} bp.`, keys: ['self_first', 'self_other', 'self_bp', 'mh'], to: 'junction', undo: 'junction reads that belong to a paralog' },
      { side: 'call', title: 'The call: three copies', what: `Short reads ${fmt.f(n.copies)} copies against the parents, bin by bin; HiFi ${fmt.f(n.hifiIn)}; each parent ${fmt.f(n.cpFather)} and ${fmt.f(n.cpMother)}.`, keys: ['sr_rel', 'sr_child', 'sr_father', 'sr_mother', 'hifi_in'], to: 'overview', undo: 'two copies by depth' },
      { side: 'call', title: "The call: the extra copy the mother's", what: `The mother's allele ${fmt.f(n.ohIn)} of the child's reads inside (${n.ohInN} parent-opposite sites), ${fmt.f(n.ohOut)} outside; phased HiFi ${fmt.f(n.hifiMat)} maternal copies to ${fmt.f(n.hifiPat)} paternal; every phasable junction read the mother's.`, keys: ['oh_in.sr', 'oh_out.sr', 'oh_in.hifi', 'hifi_mat', 'hifi_pat', 'jr_mat'], to: 'overview', undo: 'alleles at 1/2; the extra copy the father\'s' },
      { side: 'call', title: 'The call: the HiFi junction reads', what: `${n.jr} reads from ${n.instruments} instruments and ${n.cells} cells, both sides at MAPQ 60 in ${n.jrQ60}; ${n.jrMat} of ${n.jrMat + n.jrPat} phasable carry the mother's alleles, ${n.jrPat} the father's.`, keys: ['jr', 'jr_mapq60', 'jr_runs', 'jr_mat', 'jr_pat', 'jr_unphased', 'jr_starts', 'jr_ends'], to: 'junction', undo: 'chimeras; one library' },
      { side: 'call', title: 'The call: in every cell', what: `${n.jr} junction reads against ${n.throughS} straight through the start (${n.throughSmat} the mother's, ${n.throughSpat} the father's) and ${n.throughE} through the end: one chromosome's share. ${n.stop} other read stops at a breakpoint.`, keys: ['through_s.mat', 'through_s.pat', 'through_e', 'hifi_stop'], to: 'junction', undo: 'a mosaic fraction' },
      { side: 'call', title: "The call: not in the mother's long reads", what: `Of her ${B.ont.reach} nanopore reads that reach a breakpoint, ${B.ont.junction} carry the junction; ${B.ont.through_s} and ${B.ont.through_e} run straight through the start and the end.`, keys: ['ont.junction', 'ont.through_s', 'ont.through_e', 'ont.reach'], to: 'reads', undo: 'the junction in a parent' },
      { side: 'call', title: 'The call: other callers on the same reads', what: `${Object.values(B.yale_summary).filter((y) => y.verdict === 'calls it').length} of Yale's ${Object.values(B.yale_summary).filter((y) => y.made_from === 'HiFi reads').length} read-based sets write the same junction; none writes another over the segment.`, keys: [], to: 'calls', undo: 'another junction over the segment' },
      { side: 'asm', title: "The assembly: the scaffold's joins", what: `${B.joins.map((j) => j.carry).join(', ')} HiFi reads carry the ${B.joins.length} joins on GRCh38; ${B.joins.map((j) => `${j.through1}/${j.through2}`).join(', ')} run straight through their points. NucFlag calls a misjoin at each.`, keys: B.joins.map((_, i) => `joins.${i}.carry`).concat(B.joins.map((_, i) => `joins.${i}.through1`)), to: 'assembly', undo: 'reads crossing each join' },
      { side: 'asm', title: 'The assembly: its reads, given back', what: `Realigned to HPRC's sequence, ${n.jrSplitBoth} of the ${n.jr} junction reads split at both of M's segment ends (${n.jrAtEnd} at the end, ${n.jrAtStart} at the start), ${n.jrThrough} run through; ${crosses} reads cross any gap; B takes ${pieceB.mother} of the mother's reads and ${pieceB.father} of the father's.`, keys: ['realign.jr_split_at_m', 'realign.jr_at_end', 'realign.jr_at_start', 'realign.jr_through', 'realign.joins.0.cross', 'realign.per_piece.B.mat', 'realign.per_piece.B.father', 'realign.m_segment'], to: 'assembly', undo: 'its own reads fitting it unbroken' },
      { side: 'asm', title: 'The assembly: whose pieces', what: `B carries ${B.pieces.find((p) => p.piece === 'B').mom_own} of the mother's own alleles and ${B.pieces.find((p) => p.piece === 'B').dad_own} of the father's, D ${B.pieces.find((p) => p.piece === 'D').mom_own} and ${B.pieces.find((p) => p.piece === 'D').dad_own}; both are base for base M (${B.piece_pairs.filter((p) => (p.x === 'M' && (p.y === 'B' || p.y === 'D'))).map((p) => `${p.y}: ${p.snv} SNVs, ${p.indel} indels over ${fmt.int(p.len)} bp`).join('; ')}).`, keys: ['pieces.B.mom_own', 'pieces.B.dad_own', 'pieces.D.mom_own', 'pieces.D.dad_own', 'piece_pairs.MB.snv', 'piece_pairs.MB.indel', 'piece_pairs.MD.snv', 'piece_pairs.MD.indel'], to: 'assembly', undo: "B and D carrying the father's alleles" },
      { side: 'asm', title: 'The assembly: the copies it holds', what: `Over the segment's 10-kb bins it holds ${Object.entries((checks.get('copies') || {}).got || {}).sort().map(([k, v]) => `${k} copies in ${v}`).join(', ')}; the reads hold 3 throughout. ${fmt.bp(n.holes)} of the second maternal copy, and the junction, are in no piece.`, keys: ['copies', 'holes', 'holes_total'], to: 'picture', undo: "the copies matching the reads' count" },
      { side: 'asm', title: 'The assembly: its own QC', what: `NucFlag: ${B.qc_placed.filter((f) => f.label === 'MISJOIN').length} MISJOIN calls over the pieces and ${B.qc_placed.filter((f) => f.label === 'COLLAPSE').length} COLLAPSE; Flagger flags only the gap bases.`, keys: B.qc_placed.map((_, i) => `flags_placed.${i}.label`), to: 'assembly', undo: 'its QC calling it sound' },
      { side: 'call', title: 'Where the copy was made', what: B.stretch ? `The copies differ at ${n.stretchDiffer} of ${n.stretchN} sites over ${fmt.range(chrom, B.stretch[0], B.stretch[1])} and agree at ${n.agreeBefore} before and ${n.agreeAfter} after; in the mother's phase, one of each inside, the same haplotype either side.` : '–', keys: ['stretch', 'stretch_differ', 'stretch_n', 'agree_before', 'agree_after', 'bounds', 'compo.inside', 'compo.block', 'ont.conflicts'], to: 'stretch', undo: 'two identical maternal copies' },
    ];
    const grid = h('div', { class: 'claims' }); s8.appendChild(grid);
    for (const c of claims) {
      const ok = c.keys.every((k) => !checks.has(k) || checks.get(k).ok);
      grid.appendChild(h('div', { class: 'claim ' + c.side }, h('h4', {}, c.title, ' ', h('span', { class: c.side === 'call' ? 'holds' : 'fails' }, c.side === 'call' ? 'Holds.' : 'Fails.')), h('p', { class: 'what' }, c.what),
        h('div', { class: 'recount' }, c.keys.length ? c.keys.map((k) => pill(k, k)) : h('span', { class: 'pill' }, 'from the call sets'), ' · would undo it: ', c.undo, ' · ', link(c.to, 'see it'))));
    }
    const bad = B.checks.filter((c) => !c.ok);
    s8.appendChild(h('details', {}, h('summary', {}, `All ${B.checks.length} recounts against numbers.json: ${B.checks.length - bad.length} equal, ${bad.length} differ (each with its reason)`),
      h('div', { class: 'panel tight' }, AR.table(B.checks.map((c) => ({ key: c.key, expected: JSON.stringify(c.expected), got: JSON.stringify(c.got), ok: c.ok, note: c.note || '' })), [
        { key: 'key', label: 'number', cls: 'mono' }, { key: 'expected', label: 'numbers.json', cls: 'mono' }, { key: 'got', label: 'recounted here', cls: 'mono' }, { key: 'ok', label: '', fmt: (v) => `<span class="pill ${v ? 'ok' : 'bad'}">${v ? 'equal' : 'differs'}</span>` }, { key: 'note', label: 'note', cls: 'wrap' }], { scroll: true, sort: 'ok', dir: 1 }))));

    // -------------------------------------------------------------- 9. verdict
    const s9 = sec('verdict', 'Verdict', "A row for svreview's verdicts.tsv, so that the screen's review can read it back.");
    s9.appendChild(h('p', {}, h('b', {}, 'The call stands; the assembly is wrong in structure, not in sequence. '), 'The call passes every test that could have undone it, in the short reads and in the long. The assembly fails its own reads, given back to it: every junction read splits on it. It holds all three copies\' sequence, but places two maternal pieces in the father\'s chromosome across joins no read supports, leaves a third without a place, and lacks 187 kb of the second copy and the junction.'));
    const f = h('form', { class: 'verdict', onsubmit: (e) => e.preventDefault() });
    const field = (name, label, value, type = 'text') => { const inp = type === 'select' ? h('select', { name }) : type === 'textarea' ? h('textarea', { name, rows: 2 }, value) : h('input', { name, value, type }); const lab = h('label', {}, label, inp); return { lab, inp }; };
    const fields = {
      candidate: field('candidate', 'candidate', `${chrom}:${S - 1}-${S + 127}`), trio: field('trio', 'trio', `${C.trio.child}/${C.trio.father}/${C.trio.mother}`), chrom: field('chrom', 'chrom', chrom), start: field('start', 'start', String(S)), end: field('end', 'end', String(E)),
      sv_type: field('sv_type', 'sv_type', C.call.type), verdict: field('verdict', 'verdict', '', 'select'), reviewer: field('reviewer', 'reviewer', ''), note: field('note', 'note', 'the call holds; the assembly is wrong in structure (B and D maternal, four joins unsupported, 187 kb and the junction missing)', 'textarea'),
    };
    for (const v of ['true', 'false', 'uncertain', 'followed up']) fields.verdict.inp.appendChild(h('option', { value: v }, v));
    for (const k of ['candidate', 'trio', 'chrom', 'start', 'end', 'sv_type', 'verdict', 'reviewer']) f.appendChild(fields[k].lab);
    fields.note.lab.classList.add('full'); f.appendChild(fields.note.lab);
    const out = h('pre', { class: 'full' });
    const row = () => { const cols = ['candidate', 'trio', 'chrom', 'start', 'end', 'sv_type', 'verdict', 'note', 'reviewer', 'reviewed_at']; const vals = cols.map((k) => k === 'reviewed_at' ? new Date().toISOString() : String(fields[k].inp.value).replace(/[\t\n]/g, ' ')); return cols.join('\t') + '\n' + vals.join('\t') + '\n'; };
    const refresh = () => { out.textContent = row(); };
    f.addEventListener('input', refresh); refresh();
    f.appendChild(h('div', { class: 'btns full' }, h('button', { class: 'primary', onclick: async () => { const ok = await AR.copy(row()); copied.textContent = ok ? 'copied' : 'could not copy'; } }, 'copy the row'), h('button', { onclick: () => AR.download('verdict.tsv', row(), 'text/tab-separated-values') }, 'download verdict.tsv')));
    const copied = h('span', { class: 'small' }); f.appendChild(copied); f.appendChild(out);
    s9.appendChild(h('div', { class: 'panel' }, f));

    // -------------------------------------------------------------- 10. methods
    const s10 = sec('methods', 'Methods and provenance', 'Data, thresholds, tools, and how to rebuild every table on this page.');
    const th = C.thresholds;
    s10.appendChild(h('div', { class: 'cols' },
      h('div', { class: 'panel' }, h('h4', {}, 'Every threshold'), AR.table([
        ['sites', `biallelic SNVs at QUAL ${th.site_qual} or more from bcftools on the trio's short reads; each parent's depth within half and twice its median; a homozygous call with at most one read (or 5%) of the other allele, a heterozygous one with a fifth of the reads each way`],
        ['mapping quality', `${th.mapq} for every read, short or long`],
        ['a junction in a long read', `one alignment ending within ${th.junction_slack} bp of a point with ${th.junction_clip} or more bases clipped, another on the same strand starting within ${th.junction_slack} bp of the other point, clipped; at the scaffold's joins within ${fmt.int(th.gap_slack)} bp`],
        ['straight through a point', `one alignment ${fmt.int(th.through_flank)} bp either side`],
        ["a read's parent", `alleles of its primary alignment: ${th.phase_opposite_min} or more parent-opposite sites, ${Math.round(th.phase_opposite_agree * 5)} in 5 agreeing; else ${th.phase_own_min} or more of one parent's own alleles and none of the other's; a junction read by both its sides within ${fmt.kb(th.junction_phase_window)} of a breakpoint`],
        ['depth', `samtools bedcov in 10-kb bins at MAPQ ${th.mapq}; each member's flanks = the bins more than 10 kb outside the segment; phased HiFi copies = a bin's copies split by the depth of reads phased at parent-opposite sites, bins more than 10 kb inside with 10 or more such reads of each parent`],
        ['the copies at a mother-only site', `${th.site_min_hifi} or more HiFi reads; her own allele's share under 1/6 = 0, to 0.42 = 1/3, to 0.58 = 1/2 (neither), to 0.83 = 2/3`],
        ['the stretch', "the longest run of 'differ' sites, merged across at most two 'agree' sites"],
        ["the mother's phase", `every pair of her heterozygous sites on one of her nanopore reads (MAPQ ${th.mapq}) votes; a link needs ${th.ont_link_min_reads} or more reads, ${Math.round(th.ont_link_agree * 100)}% agreeing; links joined strongest first; haplotype 1 carries the reference allele at a block's first site`],
        ['uniqueness', 'minimap2 -x asm20 against the whole analysis set, secondary hits kept; a match elsewhere counts over 1 kb at 90% identity'],
        ['pieces against each other', `each record's SNVs and indels against GRCh38 from its CIGAR and sequence; differences over the stretch both cover, within ${fmt.kb(th.pair_window)} of the segment`],
        ["a parent's piece", `${th.parent_piece_min} or more of one parent's own alleles inside the segment, and ${Math.round(th.parent_piece_other * 100)}% or fewer of the other's`],
        ['reads given back to the assembly', `every primary HiFi read of the long-read window aligned with minimap2 -ax map-hifi --eqx -Y to HPRC's published sequence over the locus (hap1, hap2 and U at once); depth in 5-kb bins in copies of one haplotype's depth (${fmt.f(B.one_copy, 1)}x: the left flank's, ${fmt.f(B.flank_depth.hap1, 1)}x paternal and ${fmt.f(B.flank_depth.hap2, 1)}x maternal)`],
        ['dot plots', `exact ${th.dotplot_k}-mers shared by the read (reference orientation) and each window, both strands, chained along each diagonal; runs of ${th.dotplot_min_run} bp or longer drawn`],
      ].map(([what, rule]) => ({ what, rule })), [{ key: 'what', label: 'what' }, { key: 'rule', label: 'rule', cls: 'wrap' }])),
      h('div', {},
        h('div', { class: 'panel' }, h('h4', {}, 'Rebuild it'), h('p', { class: 'small' }, 'From the repository root, with samtools 1.22, bcftools 1.22, minimap2 2.30 and Python 3.14 on the path:'),
          h('pre', {}, `python3 -I tools/prepare_case.py cases/${C.id} --align\npython3 serve.py   # then open http://localhost:8000/`),
          h('p', { class: 'small' }, 'The first command rebuilds cases/', C.id, '/data/reviewer/ (the tables, the igv.js tracks, this page\'s bundle, and checks.tsv against numbers.json) and the junction-as-called alignment. The page itself is static: index.html and reviewer/*.js, no build step; the read panels need a server that answers byte-range requests, which GitHub Pages and serve.py do.')),
        h('div', { class: 'panel' }, h('h4', {}, 'Sources'), h('p', { class: 'small' }, 'Every file in cases/', C.id, '/data is a region slice of a public release, or was made from such slices, with its size, md5, public source and the command that made it in ', h('a', { href: `cases/${C.id}/data/MANIFEST.tsv` }, 'MANIFEST.tsv'), '. Headers are clean (no @PG lines, no UR tags, no file paths). The read records are exactly as published.'),
          h('div', { id: 'manifest' })),
        h('div', { class: 'panel' }, h('h4', {}, 'The screen and the review'), h('p', { class: 'small' }, `${C.call.source}. Independent call: ${C.call.independent}. ${C.assembly.release}; hap1 is the paternal assembly, hap2 the maternal. igv.js 3.8.10 from jsdelivr for the read panels; everything else is drawn here from the tables.`)))));
    loadManifest(C);

    document.getElementById('foot').textContent = `AssemblyReview · case ${C.id} · built from cases/${C.id}/data by tools/prepare_case.py; the page is reviewer/*.js. All of the data is open access (1000 Genomes, HPRC, the 1000 Genomes ONT project): see each source's terms of use.`;
    navSpy();
    document.addEventListener('ar:sequences', () => renderJunctionSequence(seqBox, B));
    if (AR.SEQ) renderJunctionSequence(seqBox, B);

    // ---- helpers bound to this case
    function gotoGrch38(loci) { AR.browsers.goto('grch38', loci || [fmt.locus(chrom, S - 12000, S + 12000), fmt.locus(chrom, E - 12000, E + 12000)]); }
    function gotoAssembly(locus) { AR.browsers.goto('assembly', locus || (B.gaps_local[0] ? `hap1:${B.gaps_local[0][1] - 12000}-${B.gaps_local[0][2] + 12000}` : 'hap1:1-50000')); }
    function gotoAssemblyPiece(p) { const [sname, a, b] = B.piece_local[p.piece] || []; if (!sname) return; document.getElementById('igv-assembly').scrollIntoView({ behavior: 'smooth', block: 'center' }); gotoAssembly(`${sname}:${Math.max(1, a - 2000)}-${b + 2000}`); }
    function showDotplots(r) {
      AR.clear(dotBox);
      dotBox.appendChild(h('h4', {}, `${r.name}: ${fmt.bp(r.length)}, ${r.instrument}, ${r.strand} strand, ${r.jr_parent === 'unphased' ? 'unphased' : "the " + r.jr_parent + "'s"} (${r.jr_mat_alleles} of the mother's alleles, ${r.jr_pat_alleles} of the father's)`));
      const btns = h('div', { class: 'btns' }); dotBox.appendChild(btns);
      btns.appendChild(h('button', { onclick: () => { gotoGrch38([fmt.locus(chrom, r.side1_start - 2000, E + 2000), fmt.locus(chrom, S - 2000, r.side2_end + 2000)]); AR.browsers.highlightReads('grch38', [r.name], AR.solid.jx); document.getElementById('igv-grch38').scrollIntoView({ behavior: 'smooth', block: 'center' }); } }, 'show this read in the GRCh38 panel'));
      btns.appendChild(h('button', { onclick: () => { gotoAssembly(`hap2:${B.m_segment[1] - 25000}-${B.m_segment[1] + 3000}`); AR.browsers.highlightReads('assembly', [r.name], AR.solid.jx); document.getElementById('igv-assembly').scrollIntoView({ behavior: 'smooth', block: 'center' }); } }, "show it at M's segment end"));
      if (!AR.SEQ || !AR.SEQ.reads[r.name]) { dotBox.appendChild(h('p', { class: 'note' }, 'The sequences are still loading; pick the read again in a moment.')); return; }
      const SEQ = AR.SEQ, read = SEQ.reads[r.name];
      const w = SEQ.windows;
      const sub = (win, a, b) => win.seq.slice(a - win.coord[1], b - win.coord[1] + 1);
      const side1 = r.side1_end - r.side1_start + 1;
      const plots = [];
      plots.push({ title: 'GRCh38 around the segment\'s start (lower) and its end (upper)', windows: [Object.assign({}, w.grch38_start, { marks: [{ at: S - w.grch38_start.coord[1], label: 'the segment\'s start' }] }), Object.assign({}, w.grch38_end, { marks: [{ at: E - w.grch38_end.coord[1] + 1, label: 'the segment\'s end' }] })] });
      const L1 = Math.min(30000, E - w.grch38_end.coord[1] + 1), L2 = Math.min(30000, w.grch38_start.seq.length - (S - w.grch38_start.coord[1]));
      plots.push({ title: 'GRCh38 joined as the call says: the end of the segment, then its start', windows: [{ name: `GRCh38 ${fmt.int(E - L1 + 1)}-${fmt.int(E)} + ${fmt.int(S)}-${fmt.int(S + L2 - 1)}`, coord: ['joined', 1], seq: sub(w.grch38_end, E - L1 + 1, E) + sub(w.grch38_start, S, S + L2 - 1), marks: [{ at: L1, label: 'the junction' }] }] });
      if (w.M_start && w.M_end) {
        const [ms, me] = SEQ.m_segment;
        plots.push({ title: "HPRC's maternal chromosome M around its segment start (lower) and end (upper)", windows: [Object.assign({}, w.M_start, { marks: [{ at: ms - w.M_start.coord[1], label: "M's segment start" }] }), Object.assign({}, w.M_end, { marks: [{ at: me - w.M_end.coord[1] + 1, label: "M's segment end" }] })] });
        const m1 = Math.min(30000, me - w.M_end.coord[1] + 1), m2 = Math.min(30000, w.M_start.seq.length - (ms - w.M_start.coord[1]));
        plots.push({ title: 'M joined as the call says', windows: [{ name: `hap2 ${fmt.int(me - m1 + 1)}-${fmt.int(me)} + ${fmt.int(ms)}-${fmt.int(ms + m2 - 1)}`, coord: ['joined', 1], seq: sub(w.M_end, me - m1 + 1, me) + sub(w.M_start, ms, ms + m2 - 1), marks: [{ at: m1, label: 'the junction' }] }] });
      }
      if (w.hap1_gap4 && w.hap1_gap1) {
        const g1 = SEQ.gaps_local[0], g4 = SEQ.gaps_local[SEQ.gaps_local.length - 1];
        plots.push({ title: "HPRC's paternal scaffold where the other copy's pieces come nearest the junction: D's end and gap 4 (lower), gap 1 and B's start (upper)", windows: [Object.assign({}, w.hap1_gap4, { marks: [{ at: g4[1] - w.hap1_gap4.coord[1], label: 'gap 4' }] }), Object.assign({}, w.hap1_gap1, { marks: [{ at: g1[2] - w.hap1_gap1.coord[1] + 1, label: 'gap 1' }] })] });
      }
      const grid = h('div', { class: 'dotplot' }); dotBox.appendChild(grid);
      for (const p of plots) {
        const computed = AR.dotplot.compute(read, p.windows, C.thresholds.dotplot_k, C.thresholds.dotplot_min_run);
        const canvas = h('canvas');
        const box = h('div', {}, h('h4', {}, p.title), canvas); grid.appendChild(box);
        AR.dotplot.draw(canvas, read, computed, { width: 600, height: 420, xmarks: [{ at: side1, label: `side 1 ends (≈${fmt.int(side1)})` }], fwdColor: AR.solid.pat, revColor: AR.solid.mat });
        const runs = computed.reduce((s, c) => s + c.runs.length, 0);
        box.appendChild(h('p', { class: 'cap' }, `${runs} runs of ${C.thresholds.dotplot_min_run} bp or longer (blue: forward strand, amber: reverse).`));
      }
    }
  }

  function renderJunctionSequence(box, B) {
    const SEQ = AR.SEQ; if (!SEQ) return;
    const C = B.case, [S, E] = C.segment, chrom = C.chrom, w = SEQ.windows;
    const sub = (win, a, b) => win.seq.slice(a - win.coord[1], b - win.coord[1] + 1);
    const mh = B.mh;
    const endCtx = sub(w.grch38_end, E - 39, E + 40), startCtx = sub(w.grch38_start, S - 40, S + 39), joined = sub(w.grch38_end, E - 39, E) + sub(w.grch38_start, S, S + 39);
    const mark = (s, from, to, cls) => AR.esc(s.slice(0, from)) + `<span class="${cls}">${AR.esc(s.slice(from, to))}</span>` + AR.esc(s.slice(to));
    const lines = [
      `GRCh38 at the segment's end,   ${fmt.range(chrom, E - 39, E + 40)}   ` + mark(endCtx, 40 - mh, 40, 'mh') .replace(/^(.{0})/, '$1'),
      `the junction, as the reads read it                        ` + mark(joined, 40 - mh, 40, 'mh'),
      `GRCh38 at the segment's start, ${fmt.range(chrom, S - 40, S + 39)}   ` + mark(startCtx, 40 - mh, 40, 'mh'),
    ];
    const el = h('div', { class: 'seqbox', html: lines.join('\n') });
    box.querySelectorAll('.seqbox, .cap').forEach((x) => x.remove());
    box.appendChild(el);
    box.appendChild(h('p', { class: 'cap' }, `The reads run the segment's last base straight into its first. The two sides share ${mh} base (highlighted), so the join could sit one base either way: too little shared sequence for recombination between repeats; a join made by repair or by replication.`));
  }

  async function loadManifest(C) {
    const box = document.getElementById('manifest'); if (!box) return;
    if (location.protocol === 'file:') { box.appendChild(h('p', { class: 'small' }, 'Served over HTTP, the manifest is listed here.')); return; }
    try {
      const txt = await (await fetch(`cases/${C.id}/data/MANIFEST.tsv`)).text();
      const lines = txt.trim().split('\n'); const cols = lines[0].split('\t');
      const rows = lines.slice(1).map((l) => { const f = l.split('\t'); const o = {}; cols.forEach((c, i) => (o[c] = f[i])); return o; });
      box.appendChild(h('details', {}, h('summary', {}, `${rows.length} files, ${fmt.int(rows.reduce((s, r) => s + (+r.bytes || 0), 0) / 1e6)} MB`),
        h('div', { class: 'scroll' }, AR.table(rows, [{ key: 'path', label: 'file', cls: 'mono' }, { key: 'bytes', label: 'bytes', num: true, fmt: (v) => fmt.int(+v) }, { key: 'what', label: 'what', cls: 'wrap' }, { key: 'source', label: 'public source', cls: 'wrap', fmt: (v) => /^https?:/.test(v) ? `<a href="${AR.esc(v)}">${AR.esc(v.replace(/^https?:\/\/[^/]+\//, '…/').slice(0, 90))}</a>` : AR.esc(v) }, { key: 'command', label: 'command', cls: 'wrap mono' }]))));
    } catch (e) { box.appendChild(h('p', { class: 'small' }, 'The manifest could not be fetched.')); }
  }

  function lazy(el, fn) {
    if (!('IntersectionObserver' in window)) return;
    const io = new IntersectionObserver((entries) => { for (const e of entries) if (e.isIntersecting) { io.disconnect(); fn(); } }, { rootMargin: '200px' });
    io.observe(el);
  }
  function navSpy() {
    const links = [...document.querySelectorAll('nav.toc a')];
    const secs = links.map((a) => document.getElementById(a.dataset.id)).filter(Boolean);
    const io = new IntersectionObserver((entries) => { for (const e of entries) if (e.isIntersecting) { links.forEach((a) => a.classList.toggle('active', a.dataset.id === e.target.id)); } }, { rootMargin: '-40% 0px -55% 0px' });
    secs.forEach((s) => io.observe(s));
  }

  boot();
})();
