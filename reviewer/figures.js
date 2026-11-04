/* AssemblyReview: the figures. Each is an SVG drawn from the case's tables, with lettered panels, axes and a
   caption written by the page. One palette throughout: paternal blue, maternal orange, junction and missing sequence
   magenta, maternal haplotype 1 teal and haplotype 2 purple, unphased grey. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});
  const F = (AR.fig = {});
  const { svg, h, fmt } = AR;
  const W = 1100;
  const C = (AR.pal = {
    pat: '#2b6cb0', mat: '#d98c1f', jx: '#c2185b', h1: '#138d75', h2: '#7e57c2', grey: '#9e9e9e', ink: '#222222', muted: '#6b6b6b',
    grid: '#e4e4e4', axis: '#444444', band: '#f6efe4', band2: '#efe9f7', white: '#ffffff', gap: '#111111',
  });
  const parentColor = (p) => (/^(mother|maternal|mat)$/.test(p) ? C.mat : /^(father|paternal|pat)$/.test(p) ? C.pat : C.grey);
  const FONT = 'font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif';

  // ------------------------------------------------------------------------------------------- helpers
  function figure(container, height, width = W) {
    const el = svg('svg', { class: 'fig', viewBox: `0 0 ${width} ${height}`, preserveAspectRatio: 'xMidYMin meet', style: FONT });
    container.appendChild(el);
    return el;
  }
  const text = (x, y, s, a = {}) => svg('text', Object.assign({ x, y, text: s, 'font-size': 11, fill: C.ink }, a));
  const line = (x1, y1, x2, y2, a = {}) => svg('line', Object.assign({ x1, y1, x2, y2, stroke: C.axis, 'stroke-width': 1 }, a));
  function panelLabel(g, x, y, letter, title) {
    g.appendChild(text(x, y, letter, { 'font-weight': 700, 'font-size': 13 }));
    if (title) g.appendChild(text(x + 18, y, title, { 'font-size': 11.5, fill: C.ink }));
  }
  function axisY(g, scale, x, ticks, opts = {}) {
    const [r0, r1] = scale.range;
    g.appendChild(line(x, Math.min(r0, r1), x, Math.max(r0, r1)));
    for (const t of ticks) {
      const y = scale(t);
      g.appendChild(line(x - 4, y, x, y));
      g.appendChild(text(x - 6, y + 3.5, opts.fmt ? opts.fmt(t) : String(t), { 'text-anchor': 'end', 'font-size': 10, fill: C.muted }));
      if (opts.grid) g.appendChild(line(opts.grid[0], y, opts.grid[1], y, { stroke: C.grid, 'stroke-dasharray': '2 3' }));
    }
    if (opts.label) g.appendChild(text(0, 0, opts.label, { transform: `translate(${x - (opts.labelOffset || 34)},${(r0 + r1) / 2}) rotate(-90)`, 'text-anchor': 'middle', 'font-size': 10.5, fill: C.muted }));
  }
  function axisX(g, scale, y, ticks, opts = {}) {
    const [d0, d1] = scale.range;
    g.appendChild(line(d0, y, d1, y));
    for (const t of ticks) {
      const x = scale(t);
      g.appendChild(line(x, y, x, y + 4));
      g.appendChild(text(x, y + 15, opts.fmt ? opts.fmt(t) : String(t), { 'text-anchor': 'middle', 'font-size': 10, fill: C.muted }));
    }
    if (opts.label) g.appendChild(text((d0 + d1) / 2, y + 30, opts.label, { 'text-anchor': 'middle', 'font-size': 10.5, fill: C.muted }));
  }
  const mb = (p) => (p / 1e6).toFixed(2);
  const hover = (node, html) => { AR.hover(node, html); node.classList.add('hit'); return node; };
  const arc = (g, xa, xb, yb, lift, color, label, tip) => {
    const p = svg('path', { d: `M${xa},${yb} Q${(xa + xb) / 2},${yb - lift} ${xb},${yb}`, fill: 'none', stroke: color, 'stroke-width': 1.3 });
    if (tip) hover(p, tip);
    g.appendChild(p);
    if (label) g.appendChild(text((xa + xb) / 2, yb - lift / 2 - 3, label, { 'text-anchor': 'middle', 'font-size': 10.5, fill: color }));
  };
  F.legend = function (container, items) {
    const el = h('div', { class: 'legend' });
    for (const [label, color, shape] of items) el.appendChild(h('span', {}, h('i', { class: shape || 'sq', style: { background: color, borderColor: color } }), label));
    container.appendChild(el);
    return el;
  };
  F.siteTip = function (B, s) {
    const chrom = B.case.chrom;
    const kind = { OH: 'parents homozygous for opposite alleles', MH: 'mother heterozygous, father homozygous', FH: 'father heterozygous, mother homozygous', BH: 'both parents heterozygous' }[s.kind];
    return `<b>${fmt.pos(chrom, s.pos)}</b> ${s.ref}/${s.alt} · ${kind}${s.in_segment ? ' · within the segment' : ' · outside the segment'}${s.in_stretch ? ', within the discordant interval' : ''}<br>` +
      `genotypes (ref/alt reads): proband ${s.child_gt} (${s.child_ref}/${s.child_alt}), father ${s.father_gt} (${s.father_ref}/${s.father_alt}), mother ${s.mother_gt} (${s.mother_ref}/${s.mother_alt})<br>proband HiFi ${s.hifi_ref}/${s.hifi_alt}` +
      (s.share_sr != null ? `<br>maternal allele ${s.mother_allele}: fraction ${fmt.f(s.share_sr)} (short reads), ${fmt.f(s.share_hifi)} (HiFi)` : '') +
      (s.copies != null ? `<br>copies carrying the maternal-specific allele: ${s.copies} of 3 → the maternal copies ${s.copies === 1 ? 'differ' : 'agree'}` : '') +
      (s.ref_on_hap ? `<br>maternal phase: ${s.ref} on haplotype ${s.ref_on_hap} (block ${s.block})` : '') +
      ['A', 'B', 'C', 'D', 'E', 'M', 'U'].filter((l) => s['base_' + l]).map((l) => ` ${l}:${s['base_' + l]}`).join('').replace(/^(.)/, '<br>assembly bases:$1');
  };

  // ------------------------------------------------------------------------------------------- Figure 1: the locus
  F.locus = function (container, B, opts = {}) {
    const Cs = B.case, [S, E] = Cs.segment, chrom = Cs.chrom, [W0, W1] = Cs.windows.WIN;
    const depth = AR.rows(B.depth_grch38), sites = AR.rows(B.sites), pieces = B.pieces, calls = B.calls, reads = AR.rows(B.reads);
    const Mletter = Cs.assembly.maternal_chromosome;
    const jr = reads.filter((r) => r.junction);
    const jrMat = AR.count(jr, (r) => r.jr_parent === 'mother'), jrPat = AR.count(jr, (r) => r.jr_parent === 'father');
    const thrS = AR.count(reads, (r) => r.through_S), thrE = AR.count(reads, (r) => r.through_E);
    const perPiece = {};
    for (const r of reads) { if (!r.asm_piece) continue; const c = (perPiece[r.asm_piece] = perPiece[r.asm_piece] || { reads: 0, mother: 0, father: 0, unphased: 0, q0: 0 }); c.reads++; c[r.parent]++; if (r.asm_mapq === 0) c.q0++; }
    const vsM = {}; for (const a of B.asm_vs_asm || []) vsM[a.piece] = a;
    const x0 = 70, x1 = 1080;
    const panels = [
      { id: 'sr', h: 118, letter: 'a', title: 'Read-depth copy number, short reads (trio; 10-kb bins, flank = 2 copies)' },
      { id: 'hifi', h: 118, letter: 'b', title: 'Read-depth copy number, HiFi (proband): all reads, and by parental origin of the reads' },
      { id: 'oh', h: 84, letter: 'c', title: 'Maternal allele fraction at SNVs where the parents are homozygous for opposite alleles (proband short reads)' },
      { id: 'mh', h: 78, letter: 'd', title: 'Maternal-specific allele fraction at SNVs heterozygous only in the mother (proband HiFi), within the segment' },
      { id: 'asmref', h: 176, letter: 'e', title: 'HPRC assembly aligned to GRCh38 (Winnowmap): contig position against reference position; each gap-free segment a diagonal' },
      { id: 'asmreads', h: 86, letter: 'f', title: 'Parental origin of the proband HiFi reads realigned to each contig segment (n maternal / n paternal)' },
      { id: 'asmasm', h: 74, letter: 'g', title: 'Identity of each hap1 segment and U to hap2, aligned directly (minimap2 asm20)' },
      { id: 'calls', h: 70, letter: 'h', title: 'Structural-variant records over the segment, 22 call sets' },
    ];
    let H = 14; const py = {}; for (const p of panels) { py[p.id] = H; H += p.h; } H += 44;
    const el = figure(container, H);
    let d0 = W0, d1 = W1;
    const body = svg('g'); el.appendChild(body);
    const brush = svg('rect', { class: 'brush', y: 0, height: H, x: 0, width: 0, visibility: 'hidden' }); el.appendChild(brush);
    const cross = svg('line', { class: 'cross', y1: 0, y2: H - 40, x1: 0, x2: 0, visibility: 'hidden' }); el.appendChild(cross);
    const X = () => AR.scale(d0, d1, x0, x1);
    const clipId = `clip-${Math.random().toString(36).slice(2)}`;
    el.appendChild(svg('defs', {}, svg('clipPath', { id: clipId }, svg('rect', { x: x0, y: 0, width: x1 - x0, height: H }))));
    const rowOf = { A: 0, C: 0, E: 0, B: 1, D: 1, M: 2, U: 3 };
    const pieceTip = (p) => { const c = perPiece[p.piece] || {}; const a = vsM[p.piece];
      return `<b>Segment ${p.piece}</b> (${p.parent} origin by informative alleles) ${p.contig}:${fmt.int(p.qs)}–${fmt.int(p.qe)}<br>GRCh38 ${fmt.range(chrom, p.rs, p.re)}<br>parent-specific alleles within the duplicated segment: maternal ${p.mom_own}, paternal ${p.dad_own}` +
        `<br>realigned HiFi reads: ${fmt.int(c.reads)} (maternal ${fmt.int(c.mother)}, paternal ${fmt.int(c.father)}, unphased ${fmt.int(c.unphased)}; MAPQ 0: ${fmt.int(c.q0)})` +
        (a ? `<br>identity to hap2 (direct alignment): ${a.identity === 1 ? '100.000%' : fmt.pct(a.identity, 3)} over ${fmt.bp(a.blen)}` : '') + '<br><span class="small">click: open in the assembly view</span>'; };

    function draw() {
      AR.clear(body);
      const x = X();
      const inView = (a, b) => b >= d0 && a <= d1;
      const g = svg('g', { 'clip-path': `url(#${clipId})` });
      // bands first: the duplicated segment and the discordant interval, across every panel
      const bandTop = py.sr, bandBot = H - 40;
      g.appendChild(svg('rect', { x: x(S), y: bandTop, width: Math.max(1, x(E) - x(S)), height: bandBot - bandTop, fill: C.band }));
      if (B.stretch) g.appendChild(svg('rect', { x: x(B.stretch[0]), y: bandTop, width: Math.max(1, x(B.stretch[1]) - x(B.stretch[0])), height: bandBot - bandTop, fill: C.band2 }));
      body.appendChild(g);
      for (const p of panels) {
        panelLabel(body, 8, py[p.id] + 12, p.letter, p.title);
        body.appendChild(line(x0, py[p.id] + p.h - 0.5, x1, py[p.id] + p.h - 0.5, { stroke: C.grid }));
      }
      const vis = depth.filter((d) => inView(d.bin0 + 1, d.bin1));
      // a. short reads
      let top = py.sr, base = top + panels[0].h - 10, Y = AR.scale(0, 4, base, top + 36);
      axisY(body, Y, x0, [0, 1, 2, 3, 4], { grid: [x0, x1], label: 'copies', labelOffset: 30 });
      const path = (key, Yf, stroke, dash, width = 1.3) => { let dd = ''; for (const d of vis) { if (d[key] == null) continue; dd += (dd ? 'L' : 'M') + `${x((d.bin0 + d.bin1) / 2).toFixed(1)},${Yf(Math.min(4, d[key])).toFixed(1)}`; } g.appendChild(svg('path', { d: dd, fill: 'none', stroke, 'stroke-width': width, 'stroke-dasharray': dash || null })); };
      let sd = ''; for (const d of vis) { if (d.asm_copies == null) continue; sd += `${sd ? 'L' : 'M'}${x(d.bin0 + 1).toFixed(1)},${Y(d.asm_copies).toFixed(1)} L${x(d.bin1).toFixed(1)},${Y(d.asm_copies).toFixed(1)}`; }
      g.appendChild(svg('path', { d: sd, fill: 'none', stroke: C.jx, 'stroke-width': 1.2, 'stroke-opacity': .8, 'stroke-dasharray': '1 0' }));
      path('cp_father', Y, C.pat, null, 1.1); path('cp_mother', Y, C.mat, null, 1.1); path('cp_child', Y, C.ink, null, 1.5);
      const srC = B.sr.child, srP = B.sr.father.pairs + B.sr.mother.pairs + B.sr.father.split + B.sr.mother.split;
      arc(g, x(E), x(S), top + 34, 20, C.ink, `${srC.pairs} discordant pairs and ${srC.split} split reads join the two breakpoints in the proband; ${srP} in the parents`,
        `<b>Junction in short reads</b>: ${srC.pairs} pairs in tandem-duplication orientation and ${srC.split} split alignments in ${Cs.trio.child}; none in either parent. Reads clipped exactly at the start / end: ${srC.clip_s} / ${srC.clip_e}.`);
      // b. HiFi
      top = py.hifi; base = top + panels[1].h - 10; const Yh = AR.scale(0, 4, base, top + 36);
      axisY(body, Yh, x0, [0, 1, 2, 3, 4], { grid: [x0, x1], label: 'copies', labelOffset: 30 });
      path('cp_hifi', Yh, C.ink, '4 3', 1.2);
      for (const d of vis) if (d.cp_hifi_mat != null) { g.appendChild(svg('circle', { cx: x((d.bin0 + d.bin1) / 2), cy: Yh(Math.min(4, d.cp_hifi_mat)), r: 2.4, fill: C.mat })); g.appendChild(svg('circle', { cx: x((d.bin0 + d.bin1) / 2), cy: Yh(Math.min(4, d.cp_hifi_pat)), r: 2.4, fill: C.pat })); }
      arc(g, x(E), x(S), top + 34, 20, C.jx, `${jr.length} junction-spanning HiFi reads: ${jrMat} maternal, ${jrPat} paternal, ${jr.length - jrMat - jrPat} unphased`,
        `<b>Junction-spanning HiFi reads</b>: ${jr.length} reads align up to ${fmt.pos(chrom, E)} and continue from ${fmt.pos(chrom, S)}, both alignments MAPQ 60. Reads continuous across the start: ${thrS}; across the end: ${thrE}.`);
      for (const [p, n, lab] of [[S, thrS, 'start'], [E, thrE, 'end']]) g.appendChild(hover(text(x(p) + (p === S ? 4 : -4), base - 3, `${n} reads continuous across the ${lab}`, { 'text-anchor': p === S ? 'start' : 'end', 'font-size': 9.5, fill: C.muted }), `${n} HiFi reads align continuously across ${fmt.pos(chrom, p)} with 1 kb on either side${p === S ? ` (${AR.count(reads, (r) => r.through_S && r.parent === 'mother')} maternal, ${AR.count(reads, (r) => r.through_S && r.parent === 'father')} paternal)` : ''}.`));
      for (const d of vis) {
        const r = svg('rect', { x: x(d.bin0 + 1), y: py.sr + 36, width: Math.max(1, x(d.bin1) - x(d.bin0 + 1)), height: panels[0].h + panels[1].h - 46, fill: 'transparent' });
        AR.hover(r, () => `<b>${fmt.range(chrom, d.bin0 + 1, d.bin1)}</b>${d.inside ? ' · within the segment' : d.flank ? ' · flank' : ''}<br>short-read copy number: proband ${fmt.f(d.cp_child)} (${fmt.f(d.cp_rel)} relative to the parents), father ${fmt.f(d.cp_father)}, mother ${fmt.f(d.cp_mother)}<br>` +
          `HiFi copy number ${fmt.f(d.cp_hifi)}${d.cp_hifi_mat != null ? `; maternal ${fmt.f(d.cp_hifi_mat)}, paternal ${fmt.f(d.cp_hifi_pat)} (${d.phased_reads_mat} + ${d.phased_reads_pat} phased reads)` : ''}<br>assembly: ${d.asm_copies} ${d.asm_copies === 1 ? 'copy' : 'copies'} of this bin`);
        g.appendChild(r);
      }
      // c. maternal allele fraction
      top = py.oh; base = top + panels[2].h - 8; const Yo = AR.scale(0, 1, base, top + 22);
      axisY(body, Yo, x0, [0, 0.5, 1], { fmt: (v) => v.toFixed(1), grid: [x0, x1], label: 'fraction', labelOffset: 34 });
      g.appendChild(line(x0, Yo(2 / 3), x1, Yo(2 / 3), { stroke: C.mat, 'stroke-dasharray': '2 3', 'stroke-opacity': .6 }));
      for (const s of sites) if (inView(s.pos, s.pos) && s.kind === 'OH' && s.share_sr != null) g.appendChild(hover(svg('circle', { cx: x(s.pos), cy: Yo(s.share_sr), r: 2.2, fill: s.in_segment ? C.mat : C.grey, 'fill-opacity': .85 }), () => F.siteTip(B, s)));
      // d. maternal-specific allele fraction
      top = py.mh; base = top + panels[3].h - 8; const Ym = AR.scale(0, 1, base, top + 22);
      axisY(body, Ym, x0, [0, 1 / 3, 2 / 3, 1], { fmt: (v) => (v === 1 / 3 ? '1/3' : v === 2 / 3 ? '2/3' : String(v)), grid: [x0, x1], label: 'fraction', labelOffset: 34 });
      for (const s of sites) if (inView(s.pos, s.pos) && s.kind === 'MH' && s.share_hifi != null && s.in_segment) g.appendChild(hover(svg('circle', { cx: x(s.pos), cy: Ym(s.share_hifi), r: 2.1, fill: s.copies === 1 ? C.jx : s.copies == null ? C.grey : C.mat, 'fill-opacity': .85 }), () => F.siteTip(B, s)));
      // e. assembly vs reference
      top = py.asmref; const aH = panels[4].h;
      const hap1 = pieces.filter((p) => p.hap === 'hap1').sort((a, b) => a.qs - b.qs), hap2 = pieces.filter((p) => p.hap !== 'hap1');
      const b1top = top + 24, b1bot = top + aH - 58, b2top = top + aH - 50, b2bot = top + aH - 8;
      const qlo = Math.min(...hap1.map((p) => p.qs)), qhi = Math.max(...hap1.map((p) => p.qe));
      const Y1 = AR.scale(qlo, qhi, b1bot, b1top);
      axisY(body, Y1, x0, AR.ticks(qlo, qhi, 4), { fmt: mb, label: 'hap1 (Mb)', labelOffset: 44 });
      const M = hap2.find((p) => p.piece === Mletter);
      const Y2 = AR.scale(M.qs, M.qe, b2bot, b2top);
      axisY(body, Y2, x0, [M.qs, M.qe], { fmt: mb, label: 'hap2', labelOffset: 44 });
      body.appendChild(line(x0, b2top - 4, x1, b2top - 4, { stroke: C.grid }));
      hap1.forEach((p, i) => {
        const ln = svg('line', { x1: x(p.rs), y1: Y1(p.qs), x2: x(p.re), y2: Y1(p.qe), stroke: parentColor(p.parent), 'stroke-width': 3, 'stroke-linecap': 'round' });
        hover(ln, () => pieceTip(p)); ln.addEventListener('click', (e) => { e.stopPropagation(); opts.onPiece && opts.onPiece(p); }); g.appendChild(ln);
        g.appendChild(text((x(p.rs) + x(p.re)) / 2, (Y1(p.qs) + Y1(p.qe)) / 2 - 6, p.piece, { 'text-anchor': 'middle', 'font-weight': 700, 'font-size': 10.5, fill: parentColor(p.parent) }));
        const q = hap1[i + 1];
        if (q) {
          const j = B.joins[i];
          const cn = svg('line', { x1: x(p.re), y1: Y1(p.qe), x2: x(q.rs), y2: Y1(q.qs), stroke: C.gap, 'stroke-width': 1, 'stroke-dasharray': '3 3' });
          hover(cn, j ? `<b>Gap ${i + 1}</b> (${j.left}→${j.right}): 100 N in the scaffold; on GRCh38 the next segment starts ${fmt.bp(j.back)} upstream (${fmt.pos(chrom, j.p1)} → ${fmt.pos(chrom, j.p2)}).<br>HiFi reads supporting the join: ${j.carry}; continuous across its two points: ${j.through1} / ${j.through2}; realigned reads spanning the gap: ${j.realign.cross}` : ''); g.appendChild(cn);
          g.appendChild(svg('circle', { cx: (x(p.re) + x(q.rs)) / 2, cy: (Y1(p.qe) + Y1(q.qs)) / 2, r: 3, fill: C.gap }));
        }
      });
      for (const p of hap2) {
        const isM = p.piece === Mletter;
        const Yp = isM ? Y2 : AR.scale(p.qs, p.qe, b2bot - 2, b2bot - 12);
        const ln = svg('line', { x1: x(p.rs), y1: Yp(p.qs), x2: x(p.re), y2: Yp(p.qe), stroke: parentColor(p.parent), 'stroke-width': isM ? 3 : 2.2, 'stroke-dasharray': isM ? null : '4 3', 'stroke-linecap': 'round' });
        hover(ln, () => pieceTip(p)); ln.addEventListener('click', (e) => { e.stopPropagation(); opts.onPiece && opts.onPiece(p); }); g.appendChild(ln);
        g.appendChild(text(isM ? x(Math.max(p.rs, d0)) + 6 : (x(p.rs) + x(p.re)) / 2, isM ? b2top + 10 : b2bot - 14, isM ? 'M (hap2): one copy, continuous' : 'U (unplaced)', { 'font-weight': 700, 'font-size': 10, fill: parentColor(p.parent), 'text-anchor': isM ? 'start' : 'middle' }));
      }
      // f. read origin per segment
      top = py.asmreads;
      for (const p of pieces) {
        const c = perPiece[p.piece]; if (!c || !inView(p.rs, p.re)) continue;
        const r = rowOf[p.piece] ?? 0, yy = top + 22 + r * 14;
        const xa = x(Math.max(p.rs, d0)), xb = x(Math.min(p.re, d1)), w = Math.max(2, xb - xa);
        const grp = svg('g'); let cx = xa;
        for (const [k, cc] of [['mother', C.mat], ['father', C.pat], ['unphased', C.grey]]) { const ww = w * (c[k] / c.reads); if (ww > 0) grp.appendChild(svg('rect', { x: cx, y: yy, width: ww, height: 10, fill: cc })); cx += ww; }
        if (w > 64) grp.appendChild(text(xa + 3, yy + 8.5, `${p.piece}: ${fmt.int(c.mother)} / ${fmt.int(c.father)}`, { fill: C.white, 'font-weight': 700, 'font-size': 8.5 }));
        hover(grp, () => `<b>${p.piece}</b>: ${fmt.int(c.reads)} realigned HiFi reads: maternal ${fmt.int(c.mother)}, paternal ${fmt.int(c.father)}, unphased ${fmt.int(c.unphased)}; MAPQ 0: ${fmt.int(c.q0)}${p.piece === Mletter ? `<br>junction-spanning reads: ${AR.count(jr, (r) => r.asm_split_end && r.asm_split_start)} of ${jr.length} split at both ends of M's copy; ${AR.count(jr, (r) => r.asm_through_M)} continuous` : ''}`);
        g.appendChild(grp);
      }
      // g. identity to hap2
      top = py.asmasm;
      for (const p of pieces) {
        const a = vsM[p.piece]; if (!a || !inView(p.rs, p.re)) continue;
        const r = rowOf[p.piece] ?? 0, yy = top + 22 + Math.min(r, 1) * 14 + (p.piece === 'U' ? 28 : 0);
        const xa = x(Math.max(p.rs, d0)), xb = x(Math.min(p.re, d1));
        const same = a.identity >= 0.99999;
        const rect = svg('rect', { x: xa, y: yy, width: Math.max(2, xb - xa), height: 10, fill: same ? C.mat : C.pat, 'fill-opacity': same ? 1 : .45 });
        hover(rect, `<b>${p.piece} aligned to hap2</b> (minimap2 asm20): ${same ? '100.000% identity' : fmt.pct(a.identity, 3) + ' identity'} over ${fmt.bp(a.blen)}<br>${a.q}:${fmt.int(a.qs)}–${fmt.int(a.qe)} ↔ hap2:${fmt.int(a.ts)}–${fmt.int(a.te)}`);
        g.appendChild(rect);
        if (xb - xa > 80) g.appendChild(text(xa + 3, yy + 8.5, same ? `${p.piece} = hap2, 100.000% over ${fmt.kb(a.blen)}` : `${p.piece}: ${fmt.pct(a.identity, 2)}`, { fill: same ? C.white : C.ink, 'font-weight': 700, 'font-size': 8.5 }));
      }
      // h. calls
      top = py.calls;
      const rowc = { 'the junction': 0, 'the scaffold': 1, other: 2 };
      const labc = { 'the junction': 'junction', 'the scaffold': 'scaffold artefact', other: 'other' };
      for (const [k, i] of Object.entries(rowc)) body.appendChild(text(x0 - 6, top + 30 + i * 14, labc[k], { 'text-anchor': 'end', 'font-size': 9.5, fill: C.muted }));
      for (const c of calls) {
        const a = c.pos, b = Math.max(c.end, c.pos + 1); if (!inView(a, b)) continue;
        const yy = top + 22 + (rowc[c.reflects] ?? 2) * 14;
        const color = c.reflects === 'the junction' ? C.jx : c.reflects === 'the scaffold' ? C.pat : C.grey;
        g.appendChild(hover(svg('rect', { x: x(a), y: yy, width: Math.max(2, x(b) - x(a)), height: 9, fill: color, 'fill-opacity': .7 }), `<b>${c.set}</b> (${c.made_from}) ${c.type} ${fmt.range(chrom, c.pos, c.end)} · ${fmt.bp(c.len)} · GT ${c.gt || '–'}<br>reflects ${c.reflects === 'the scaffold' ? 'the scaffold structure' : c.reflects}`));
      }
      // breakpoints through every panel; x axis
      for (const p of [S, E]) g.appendChild(line(x(p), py.sr + 20, x(p), H - 40, { stroke: C.jx, 'stroke-width': .8, 'stroke-dasharray': '3 3', 'pointer-events': 'none' }));
      axisX(body, x, H - 40, AR.ticks(d0, d1, 8), { fmt: (t) => (t / 1e6).toFixed(d1 - d0 < 300000 ? 3 : 2), label: `GRCh38 ${chrom} (Mb)` });
    }
    draw();
    let drag = null;
    const pt = (e) => { const r = el.getBoundingClientRect(); return ((e.clientX - r.left) / r.width) * W; };
    el.addEventListener('mousedown', (e) => { const px = pt(e); if (px < x0 || px > x1) return; drag = { a: px, moved: false }; brush.setAttribute('x', px); brush.setAttribute('width', 0); brush.setAttribute('visibility', 'visible'); e.preventDefault(); });
    el.addEventListener('mousemove', (e) => { const px = Math.min(x1, Math.max(x0, pt(e))); cross.setAttribute('x1', px); cross.setAttribute('x2', px); cross.setAttribute('visibility', 'visible'); if (!drag) return; drag.moved = drag.moved || Math.abs(px - drag.a) > 3; brush.setAttribute('x', Math.min(drag.a, px)); brush.setAttribute('width', Math.abs(px - drag.a)); });
    el.addEventListener('mouseleave', () => cross.setAttribute('visibility', 'hidden'));
    window.addEventListener('mouseup', (e) => { if (!drag) return; const px = Math.min(x1, Math.max(x0, pt(e))); const x = X(); brush.setAttribute('visibility', 'hidden');
      if (drag.moved) { const a = x.invert(Math.min(drag.a, px)), b = x.invert(Math.max(drag.a, px)); if (b - a > 500) { d0 = a; d1 = b; draw(); } } else if (opts.onLocus) opts.onLocus(Math.round(x.invert(px))); drag = null; });
    el.addEventListener('dblclick', () => { d0 = W0; d1 = W1; draw(); });
    return { el, zoom(a, b) { d0 = Math.max(W0, a); d1 = Math.min(W1, b); draw(); } };
  };

  // ------------------------------------------------------------------------------------------- Figure 2b: the junction-spanning reads
  F.junctionReads = function (container, B, opts = {}) {
    const Cs = B.case, [S, E] = Cs.segment, chrom = Cs.chrom;
    const reads = AR.rows(B.reads).filter((r) => r.junction);
    const siteByPos = new Map(AR.rows(B.sites).map((s) => [s.pos, s]));
    const order = { mother: 0, father: 1, unphased: 2 };
    reads.sort((a, b) => (order[a.jr_parent] - order[b.jr_parent]) || ((b.side1_end - b.side1_start) + (b.side2_end - b.side2_start)) - ((a.side1_end - a.side1_start) + (a.side2_end - a.side2_start)));
    const L = Math.max(...reads.map((r) => r.side1_end - r.side1_start + 1)), R = Math.max(...reads.map((r) => r.side2_end - r.side2_start + 1));
    const x0 = 70, x1 = 980, rowH = 9, top = 46;
    const X = AR.scale(-L - 500, R + 500, x0, x1);
    const H = top + reads.length * rowH + 60;
    const el = figure(container, H);
    panelLabel(el, 8, 14, 'b', `The ${reads.length} junction-spanning HiFi reads, aligned at the junction: bases to the left of 0 align up to ${fmt.pos(chrom, E)}, bases to the right from ${fmt.pos(chrom, S)}`);
    el.appendChild(svg('rect', { x: X(-L - 500), y: top - 4, width: X(0) - X(-L - 500), height: reads.length * rowH + 8, fill: C.band, 'fill-opacity': .6 }));
    el.appendChild(svg('rect', { x: X(0), y: top - 4, width: X(R + 500) - X(0), height: reads.length * rowH + 8, fill: C.band2, 'fill-opacity': .5 }));
    el.appendChild(text(X(0) - 6, top - 8, 'segment end', { 'text-anchor': 'end', 'font-size': 10, fill: C.muted }));
    el.appendChild(text(X(0) + 6, top - 8, 'segment start', { 'font-size': 10, fill: C.muted }));
    let y = top + 4, lastGroup = null;
    for (const r of reads) {
      if (r.jr_parent !== lastGroup) { lastGroup = r.jr_parent; el.appendChild(text(x1 + 8, y + 3, r.jr_parent === 'mother' ? `maternal (${AR.count(reads, (q) => q.jr_parent === 'mother')})` : r.jr_parent === 'father' ? 'paternal' : `unphased (${AR.count(reads, (q) => q.jr_parent === 'unphased')})`, { 'font-size': 10, fill: parentColor(r.jr_parent), 'font-weight': 600 })); }
      const col = parentColor(r.jr_parent);
      const g = svg('g');
      const l1 = r.side1_end - r.side1_start + 1, l2 = r.side2_end - r.side2_start + 1;
      g.appendChild(svg('rect', { x: X(-l1), y: y - 2.5, width: X(0) - X(-l1), height: 5, fill: col, 'fill-opacity': .35 }));
      g.appendChild(svg('rect', { x: X(0), y: y - 2.5, width: X(l2) - X(0), height: 5, fill: col, 'fill-opacity': .35 }));
      for (const [pos, base] of (B.read_alleles[r.name] || [])) {
        const s = siteByPos.get(pos); if (!s) continue;
        let xx = null;
        if (pos <= E && pos >= E - l1) xx = pos - E - 1; else if (pos >= S && pos <= S + l2) xx = pos - S + 1;
        if (xx == null) continue;
        const who = base === s.mother_allele ? 'mother' : base === s.father_allele ? 'father' : null;
        if (!who) continue;
        g.appendChild(svg('rect', { x: X(xx) - 1, y: y - 4, width: 2, height: 8, fill: parentColor(who) }));
      }
      hover(g, `<b>${r.name}</b> · ${r.instrument}, ${r.cell} · ${r.strand} strand · ${fmt.bp(r.length)}<br>to the end: ${fmt.range(chrom, r.side1_start, r.side1_end)} (MAPQ ${r.side1_mapq}); from the start: ${fmt.range(chrom, r.side2_start, r.side2_end)} (MAPQ ${r.side2_mapq})<br>informative alleles within 30 kb of a breakpoint: maternal ${r.jr_mat_alleles}, paternal ${r.jr_pat_alleles} → ${r.jr_parent === 'unphased' ? 'unphased' : r.jr_parent === 'mother' ? 'maternal' : 'paternal'}<br>on the assembly: ${r.asm_split_end && r.asm_split_start ? "split at both ends of M's copy" : r.asm_through_M ? 'continuous' : 'split at one end'}<br><span class="small">click: dot plots of this read</span>`);
      g.addEventListener('click', () => opts.onRead && opts.onRead(r));
      el.appendChild(g);
      y += rowH;
    }
    el.appendChild(line(X(0), top - 4, X(0), y, { stroke: C.jx, 'stroke-width': 1.2 }));
    axisX(el, X, y + 6, AR.ticks(-L, R, 8), { fmt: (t) => `${t / 1000}`, label: 'distance from the junction (kb)' });
    return el;
  };

  // ------------------------------------------------------------------------------------------- Figure 3a: structure schematic
  F.structure = function (container, B, opts = {}) {
    const Cs = B.case, [S, E] = Cs.segment, L = E - S + 1, chrom = Cs.chrom;
    const pieces = B.pieces, reads = AR.rows(B.reads);
    const perPiece = {}; for (const r of reads) { if (!r.asm_piece) continue; const c = (perPiece[r.asm_piece] = perPiece[r.asm_piece] || { reads: 0, mother: 0, father: 0, unphased: 0, q0: 0 }); c.reads++; c[r.parent]++; if (r.asm_mapq === 0) c.q0++; }
    const FL = 100000, x0 = 250, x1 = 1080, k = (x1 - x0) / (2 * FL + 2 * L);
    const X = (p, copy) => x0 + (FL + (copy - 1) * L + (p - S)) * k;
    const H = 392, bh = 18;
    const Y = { M: 40, mom: 104, U: 160, dad: 220, scaf: 300 };
    const el = figure(container, H);
    panelLabel(el, 8, 14, 'a', 'Structure of the locus: the proband haplotypes inferred from the reads, and the HPRC assembly of them, to scale');
    const lab = (y, main, sub) => { el.appendChild(text(x0 - 14, y + 9, main, { 'text-anchor': 'end', 'font-weight': 600, 'font-size': 11 })); if (sub) el.appendChild(text(x0 - 14, y + 22, sub, { 'text-anchor': 'end', 'font-size': 10, fill: C.muted })); };
    const bar = (xa, xb, y, fill, txt, opacity = 1) => { const g = svg('g'); g.appendChild(svg('rect', { x: Math.max(x0, xa), y, width: Math.max(1, Math.min(x1, xb) - Math.max(x0, xa)), height: bh, fill, 'fill-opacity': opacity })); if (txt && xb - xa > 40) g.appendChild(text((Math.max(x0, xa) + Math.min(x1, xb)) / 2, y + 12.5, txt, { 'text-anchor': 'middle', fill: C.white, 'font-weight': 600, 'font-size': 10 })); return g; };
    const ribbon = (x0a, x1a, ya, x0b, x1b, yb, color, op = .18) => { const cy = (ya + yb) / 2; return svg('path', { d: `M${x0a},${ya} C${x0a},${cy} ${x0b},${cy} ${x0b},${yb} L${x1b},${yb} C${x1b},${cy} ${x1a},${cy} ${x1a},${ya} Z`, fill: color, 'fill-opacity': op }); };
    lab(Y.M, 'HPRC hap2 (maternal)', Cs.assembly.sequences.hap2.contig);
    lab(Y.mom, 'Proband, maternal haplotype', 'inferred from the reads');
    lab(Y.U, 'HPRC unplaced contig U', Cs.assembly.sequences.U.contig);
    lab(Y.dad, 'Proband, paternal haplotype', 'inferred from the reads');
    lab(Y.scaf, 'HPRC hap1 (paternal) scaffold', Cs.assembly.sequences.hap1.contig);
    const M = pieces.find((p) => p.piece === Cs.assembly.maternal_chromosome), U = pieces.find((p) => p.piece === 'U');
    const scaf = pieces.filter((p) => p.hap === 'hap1').sort((a, b) => a.qs - b.qs);
    const pieceTip = (p) => { const c = perPiece[p.piece] || {}; return `<b>Segment ${p.piece}</b> · ${p.parent} origin by alleles<br>${p.contig}:${fmt.int(p.qs)}–${fmt.int(p.qe)} = GRCh38 ${fmt.range(chrom, p.rs, p.re)} (${fmt.kb(p.re - p.rs + 1)})<br>parent-specific alleles in the segment: maternal ${p.mom_own}, paternal ${p.dad_own}<br>realigned HiFi reads: ${fmt.int(c.reads)} (maternal ${fmt.int(c.mother)}, paternal ${fmt.int(c.father)}, unphased ${fmt.int(c.unphased)}; MAPQ 0: ${fmt.int(c.q0)})<br><span class="small">click: open in the assembly view</span>`; };
    const rib = svg('g'); el.appendChild(rib);
    rib.appendChild(ribbon(X(S, 1), X(E, 1), Y.M + bh, X(S, 1), X(E, 1), Y.mom, C.mat));
    rib.appendChild(ribbon(X(E, 1), X(E + FL, 1), Y.M + bh, X(E, 2), X(E + FL, 2), Y.mom, C.grey, .15));
    if (U) rib.appendChild(ribbon(X(U.rs, 2), X(U.re, 2), Y.mom + bh, X(U.rs, 2), X(U.re, 2), Y.U, C.mat));
    const A = scaf[0];
    const localOfS = A.local_start + (S - A.rs) * ((A.local_end - A.local_start) / (A.re - A.rs));
    const XS = (local) => X(S, 1) + (local - localOfS) * k;
    for (const p of scaf) {
      const copy = p.parent === 'maternal' ? 2 : 1, yT = p.parent === 'maternal' ? Y.mom + bh : Y.dad + bh;
      const lo = Math.max(p.rs, copy === 2 ? S : S - FL), hi = Math.min(p.re, copy === 2 ? E : E + FL);
      rib.appendChild(ribbon(Math.max(x0, XS(p.local_start)), Math.min(x1, XS(p.local_end)), Y.scaf, X(lo, copy), X(hi, copy), yT, p.parent === 'maternal' ? C.mat : C.pat));
    }
    el.appendChild(bar(X(S - FL, 1), X(S, 1), Y.M, C.grey, '', .5)); el.appendChild(hover(bar(X(S, 1), X(E, 1), Y.M, C.mat, 'M: one copy of the segment'), () => pieceTip(M))); el.appendChild(bar(X(E, 1), X(E + FL, 1), Y.M, C.grey, '', .5));
    el.appendChild(bar(X(S - FL, 1), X(S, 1), Y.mom, C.grey, '', .5));
    el.appendChild(hover(bar(X(S, 1), X(E, 1), Y.mom, C.mat, 'copy 1'), `<b>Proband maternal haplotype</b>: two tandem copies of ${fmt.bp(L)}, joined at the junction. Which copy is proximal cannot be determined from the reads.`));
    el.appendChild(hover(bar(X(S, 2), X(E, 2), Y.mom, C.mat, 'copy 2'), `<b>The second maternal copy</b>: represented in the assembly only by segments B, U and D; ${fmt.bp(B.holes.reduce((s, [a, b]) => s + b - a + 1, 0))} of it, the junction included, in no contig.`));
    el.appendChild(bar(X(E, 2), X(E + FL, 2), Y.mom, C.grey, '', .5));
    for (const [a, b] of B.holes) el.appendChild(hover(svg('rect', { x: X(a, 2), y: Y.mom - 3, width: Math.max(2, X(b, 2) - X(a, 2)), height: bh + 6, fill: C.jx, 'fill-opacity': .8 }), `<b>Absent from the assembly</b>: ${fmt.range(chrom, a, b)} of the second maternal copy (${fmt.bp(b - a + 1)})`));
    const jx = X(E, 1);
    el.appendChild(hover(svg('rect', { x: jx - 1.5, y: Y.mom - 8, width: 3, height: bh + 16, fill: C.jx }), `<b>Junction</b>: ${fmt.pos(chrom, E)} joined to ${fmt.pos(chrom, S)}; ${B.numbers.jr} junction-spanning HiFi reads.`));
    el.appendChild(text(jx, Y.mom - 12, 'junction', { 'text-anchor': 'middle', fill: C.jx, 'font-size': 10 }));
    if (U) { el.appendChild(hover(bar(X(U.rs, 2), X(U.re, 2), Y.U, C.mat, 'U'), () => pieceTip(U))); el.appendChild(text(X(U.re, 2) + 8, Y.U + 12.5, `${fmt.bp(U.qe - U.qs + 1)}, no chromosomal placement`, { 'font-size': 10, fill: C.muted })); }
    el.appendChild(bar(X(S - FL, 1), X(S, 1), Y.dad, C.grey, '', .5)); el.appendChild(hover(bar(X(S, 1), X(E, 1), Y.dad, C.pat, 'one copy'), `<b>Proband paternal haplotype</b>: one copy of the segment (${fmt.bp(L)}), represented by segments A, C and E of the hap1 scaffold.`)); el.appendChild(bar(X(E, 1), X(E + FL, 1), Y.dad, C.grey, '', .5));
    const gaps = B.gaps_local.filter((g) => g[0] === 'hap1');
    for (const p of scaf) { const g = hover(bar(XS(p.local_start), XS(p.local_end), Y.scaf, p.parent === 'maternal' ? C.mat : C.pat, p.piece), () => pieceTip(p)); g.addEventListener('click', () => opts.onPiece && opts.onPiece(p)); el.appendChild(g); }
    gaps.forEach((gp, i) => { const x = XS(gp[1]); const j = B.joins[i]; el.appendChild(hover(svg('rect', { x: x - 1.5, y: Y.scaf - 5, width: 3, height: bh + 10, fill: C.gap }), `<b>Gap ${i + 1}</b>: ${gp[2] - gp[1] + 1} N between ${j ? j.left : '?'} and ${j ? j.right : '?'}${j ? `; on GRCh38 the next segment starts ${fmt.bp(j.back)} upstream. Reads supporting the join: ${j.carry}; continuous across its points: ${j.through1} / ${j.through2}` : ''}`)); el.appendChild(text(x, Y.scaf + bh + 14, `gap ${i + 1}`, { 'text-anchor': 'middle', 'font-size': 9.5, fill: C.muted })); });
    for (const [n, x, y] of [[1, X((S + E) / 2, 1), Y.M - 10], [2, (XS(scaf[1].local_start) + X(scaf[1].rs, 2)) / 2, (Y.scaf + Y.mom + bh) / 2 + 6], [3, XS(gaps[0][1]) + 14, Y.scaf - 12], [4, X((B.holes[1][0] + B.holes[1][1]) / 2, 2), Y.mom - 12], [5, (XS(scaf[0].local_end) + XS(scaf[2].local_start)) / 2, Y.dad + bh + 16]]) {
      el.appendChild(svg('circle', { cx: x, cy: y, r: 7.5, fill: C.white, stroke: C.ink, 'stroke-width': 1 })); el.appendChild(text(x, y + 3.5, String(n), { 'text-anchor': 'middle', 'font-size': 10, 'font-weight': 700 }));
    }
    el.appendChild(text(x0, H - 10, 'Flanks (100 kb) in grey. The hap1 scaffold is drawn in its own order; segments coloured by the parental origin of their alleles; bands join each segment to the copy it represents.', { 'font-size': 10, fill: C.muted }));
    return el;
  };

  // ------------------------------------------------------------------------------------------- Figure 3b: hap1 against hap2
  F.asmDotplot = function (container, B) {
    const Cs = B.case, al = B.asm_vs_asm || [], len = B.local_len;
    if (!al.length) return;
    const x0 = 96, x1 = 1060, top = 26, h1 = 330, gap = 30, hU = 64, H = top + h1 + gap + hU + 60;
    const el = figure(container, H);
    panelLabel(el, 8, 14, 'b', 'Each gap-free segment of hap1, and U, aligned directly to hap2 (minimap2 asm20); dot-plot representation with identity');
    const X = AR.scale(1, len.hap2, x0, x1), Y1 = AR.scale(1, len.hap1, top + h1, top), YU = AR.scale(1, len.U || 1, top + h1 + gap + hU, top + h1 + gap);
    el.appendChild(svg('rect', { x: x0, y: top, width: x1 - x0, height: h1, fill: '#fafafa', stroke: C.grid }));
    if (len.U) el.appendChild(svg('rect', { x: x0, y: top + h1 + gap, width: x1 - x0, height: hU, fill: '#fafafa', stroke: C.grid }));
    for (const [letter, [sname, a, b]] of Object.entries(B.piece_local)) {
      const p = B.pieces.find((q) => q.piece === letter);
      if (sname === 'hap1') { el.appendChild(hover(svg('rect', { x: x0 - 14, y: Y1(b), width: 8, height: Y1(a) - Y1(b), fill: parentColor(p.parent) }), `<b>${letter}</b> hap1:${fmt.int(a)}–${fmt.int(b)} (${p.parent} origin)`)); el.appendChild(text(x0 - 18, (Y1(a) + Y1(b)) / 2 + 4, letter, { 'text-anchor': 'end', 'font-weight': 700, 'font-size': 11 })); }
    }
    B.gaps_local.forEach((gp, i) => { if (gp[0] !== 'hap1') return; el.appendChild(line(x0, Y1(gp[1]), x1, Y1(gp[1]), { stroke: C.gap, 'stroke-dasharray': '3 3', 'stroke-opacity': .6 })); el.appendChild(text(x1 + 4, Y1(gp[1]) + 3, `gap ${i + 1}`, { 'font-size': 9.5, fill: C.muted })); });
    if (B.m_segment) { const [ms, me] = B.m_segment; for (const p of [ms, me]) el.appendChild(line(X(p), top, X(p), top + h1 + gap + hU, { stroke: C.jx, 'stroke-dasharray': '4 3' })); el.appendChild(hover(svg('rect', { x: X(ms), y: top + h1 + gap + hU + 4, width: X(me) - X(ms), height: 5, fill: C.jx }), `<b>M's single copy of the segment</b>: hap2:${fmt.int(ms)}–${fmt.int(me)}`)); }
    axisX(el, X, top + h1 + gap + hU + 12, AR.ticks(1, len.hap2, 8), { fmt: mb, label: 'hap2 position (Mb)' });
    axisY(el, Y1, x0 - 22, AR.ticks(1, len.hap1, 6), { fmt: mb, label: 'hap1 position (Mb)', labelOffset: 48 });
    if (len.U) axisY(el, YU, x0 - 22, [1, len.U], { fmt: mb, label: 'U', labelOffset: 48 });
    for (const a of al) {
      const Y = a.q === 'U' ? YU : Y1, same = a.identity >= 0.99999;
      const ln = svg('line', { x1: X(a.ts), y1: Y(a.qs), x2: X(a.te), y2: Y(a.qe), stroke: same ? C.mat : C.pat, 'stroke-width': same ? 3.5 : 2.5, 'stroke-linecap': 'round' });
      hover(ln, `<b>${a.piece} to hap2</b>: ${same ? '100.000%' : fmt.pct(a.identity, 3)} identity over ${fmt.bp(a.blen)} (${fmt.int(a.matches)} matches)<br>${a.q}:${fmt.int(a.qs)}–${fmt.int(a.qe)} ↔ hap2:${fmt.int(a.ts)}–${fmt.int(a.te)}, ${a.strand} strand, MAPQ ${a.mapq}`);
      el.appendChild(ln);
      el.appendChild(text((X(a.ts) + X(a.te)) / 2 + 8, (Y(a.qs) + Y(a.qe)) / 2 - 7, `${a.piece}: ${same ? '100.000%' : fmt.pct(a.identity, 2)}`, { 'font-weight': 700, 'font-size': 10.5, fill: same ? C.mat : C.pat }));
    }
    return el;
  };

  // ------------------------------------------------------------------------------------------- Figure 3c: reads realigned to the assembly
  F.assemblyDepth = function (container, B, opts = {}) {
    const Cs = B.case, rows = AR.rows(B.depth_assembly), len = B.local_len, seqs = ['hap1', 'hap2', 'U'].filter((s) => len[s]);
    const x0 = 160, x1 = 1070, k = (x1 - x0) / Math.max(...seqs.map((s) => len[s]));
    const rowH = 124, H = 20 + seqs.length * rowH + 30;
    const el = figure(container, H);
    panelLabel(el, 8, 14, 'c', `Proband HiFi reads realigned to the assembly: depth by parental origin of the reads (5-kb bins, in multiples of the haploid depth, ${fmt.f(B.one_copy, 1)}×)`);
    const labels = { hap1: ['hap1 (paternal)', Cs.assembly.sequences.hap1.contig], hap2: ['hap2 (maternal)', Cs.assembly.sequences.hap2.contig], U: ['U (unplaced)', Cs.assembly.sequences.U.contig] };
    const M = B.pieces.find((p) => p.piece === Cs.assembly.maternal_chromosome);
    seqs.forEach((seq, si) => {
      const top = 26 + si * rowH, base = top + 72, Y = AR.scale(0, 2.4, base, top + 6);
      const X = (p) => x0 + (p - 1) * k;
      el.appendChild(text(x0 - 30, top + 30, labels[seq][0], { 'text-anchor': 'end', 'font-weight': 600, 'font-size': 10.5 }));
      el.appendChild(text(x0 - 30, top + 43, labels[seq][1], { 'text-anchor': 'end', 'font-size': 9.5, fill: C.muted }));
      axisY(el, Y, x0 - 2, [0, 1, 2], { grid: [x0, X(len[seq])] });
      if (seq === 'hap2' && M) for (const [a, b] of B.holes) { if (b - a < 10000) continue; const la = M.local_start + (a - M.rs) * ((M.local_end - M.local_start) / (M.re - M.rs)), lb = M.local_start + (b - M.rs) * ((M.local_end - M.local_start) / (M.re - M.rs));
        el.appendChild(hover(svg('rect', { x: X(la), y: top + 2, width: Math.max(2, X(lb) - X(la)), height: base - top - 2, fill: C.jx, 'fill-opacity': .1 }), `<b>No contig represents the second copy here</b> (GRCh38 ${fmt.range(Cs.chrom, a, b)}): its reads align to the best available sequence, M or C.`)); }
      for (const q of B.qc_rows.filter((r) => r[0] === 'nucflag' && Cs.assembly.sequences[seq].contig === r[1])) { const off = Cs.assembly.sequences[seq].local_offset, a = q[2] + 1 - off, b = q[3] - off; if (b < 1 || a > len[seq]) continue;
        el.appendChild(hover(svg('rect', { x: X(Math.max(1, a)), y: base + 13, width: Math.max(2, X(Math.min(len[seq], b)) - X(Math.max(1, a))), height: 3, fill: q[4] === 'COLLAPSE' ? '#b71c1c' : '#f59e0b' }), `<b>NucFlag ${q[4]}</b> ${q[1]}:${fmt.int(q[2] + 1)}–${fmt.int(q[3])}`)); }
      for (const r of rows) {
        if (r.seq !== seq) continue;
        const xa = X(r.bin0 + 1), w = Math.max(0.8, X(r.bin1) - xa);
        const g = svg('g'); let y = base;
        for (const [v, c] of [[Math.min(2.4, r.mat), C.mat], [Math.min(2.4, r.pat), C.pat], [Math.min(2.4, r.unphased), C.grey]]) { if (v <= 0) continue; const hh = base - Y(v); g.appendChild(svg('rect', { x: xa, y: y - hh, width: w, height: hh, fill: c })); y -= hh; }
        if (r.q0 > 0.02) g.appendChild(svg('rect', { x: xa, y: Y(Math.min(2.4, r.q0)) - 1, width: w, height: 2, fill: C.ink }));
        hover(g, `<b>${seq}:${fmt.int(r.bin0 + 1)}–${fmt.int(r.bin1)}</b> · segment ${r.piece || '–'}<br>depth (× haploid): maternal ${fmt.f(r.mat)}, paternal ${fmt.f(r.pat)}, unphased ${fmt.f(r.unphased)} (total ${fmt.f(r.mat + r.pat + r.unphased)}); MAPQ 0: ${fmt.f(r.q0)}`);
        g.addEventListener('click', () => opts.onLocus && opts.onLocus(seq, r.bin0 + 1, r.bin1));
        el.appendChild(g);
      }
      for (const [letter, [sname, a, b]] of Object.entries(B.piece_local)) {
        if (sname !== seq) continue; const p = B.pieces.find((q) => q.piece === letter);
        const r = hover(svg('rect', { x: X(a), y: base + 2, width: Math.max(2, X(b) - X(a)), height: 9, fill: parentColor(p.parent) }), `<b>Segment ${letter}</b> ${seq}:${fmt.int(a)}–${fmt.int(b)} (${p.parent} origin)`);
        r.addEventListener('click', () => opts.onLocus && opts.onLocus(seq, a, b)); el.appendChild(r);
        if (X(b) - X(a) > 12) el.appendChild(text((X(a) + X(b)) / 2, base + 9.5, letter, { 'text-anchor': 'middle', fill: C.white, 'font-weight': 700, 'font-size': 8.5 }));
      }
      B.gaps_local.forEach((gp, i) => { if (gp[0] !== seq) return; const j = B.joins[i];
        el.appendChild(line(X(gp[1]), top, X(gp[1]), base + 11, { stroke: C.gap, 'stroke-dasharray': '3 2' }));
        el.appendChild(hover(text(X(gp[1]), top - 2, `gap ${i + 1}: ${j ? j.realign.cross : 0} reads span`, { 'text-anchor': 'middle', 'font-size': 9.5 }), j ? `<b>Gap ${i + 1}</b> (${j.left}→${j.right}), ${seq}:${fmt.int(gp[1])}–${fmt.int(gp[2])}<br>realigned reads spanning it: ${j.realign.cross}<br>depth within 2 kb of its edges: ${fmt.f(j.realign.edge_left)} / ${fmt.f(j.realign.edge_right)} × haploid` : '')); });
      if (seq === 'hap2' && B.m_segment) { const [ms, me] = B.m_segment; const jr = AR.rows(B.reads).filter((r) => r.junction);
        for (const p of [ms, me]) el.appendChild(line(X(p), top + 6, X(p), base, { stroke: C.jx, 'stroke-width': 1.3 }));
        el.appendChild(hover(svg('path', { d: `M${X(ms)},${top + 6} C${X(ms)},${top - 14} ${X(me)},${top - 14} ${X(me)},${top + 6}`, fill: 'none', stroke: C.jx, 'stroke-width': 1.1, 'stroke-dasharray': '5 3' }), `<b>${jr.length} junction-spanning reads realigned to M</b>: ${AR.count(jr, (r) => r.asm_split_end && r.asm_split_start)} split at both ends of M's copy, ${AR.count(jr, (r) => r.asm_through_M)} continuous.`));
        el.appendChild(text((X(ms) + X(me)) / 2, top - 4, `${jr.length} junction-spanning reads: split at both ends of the single copy`, { 'text-anchor': 'middle', 'font-size': 10, fill: C.jx })); }
      axisX(el, AR.scale(1, len[seq], x0, X(len[seq])), base + 14, AR.ticks(1, len[seq], seq === 'U' ? 3 : 8), { fmt: mb });
    });
    return el;
  };

  // ------------------------------------------------------------------------------------------- Figure 4: haplotype composition of the maternal copies
  F.composition = function (container, B) {
    const Cs = B.case, chrom = Cs.chrom, [S, E] = Cs.segment, sites = AR.rows(B.sites), reads = AR.rows(B.reads);
    const readByName = new Map(reads.map((r) => [r.name, r]));
    const st = B.stretch, bounds = B.bounds;
    const blk = B.blocks.find((b) => st && b.start <= st[0] && b.end >= st[1]);
    const x0 = 70, x1 = 1080;
    // a. along the segment: copies carrying the maternal-specific allele
    const el = figure(container, 230);
    panelLabel(el, 8, 14, 'a', 'Copies carrying the maternal-specific allele, at SNVs heterozygous only in the mother, along the segment (proband HiFi; 3 copies in all)');
    const X = AR.scale(S - 20000, E + 20000, x0, x1), Ya = AR.scale(-0.5, 2.5, 190, 40);
    if (st) el.appendChild(svg('rect', { x: X(st[0]), y: 30, width: X(st[1]) - X(st[0]), height: 170, fill: C.band2 }));
    axisY(el, Ya, x0, [0, 1, 2], { grid: [x0, x1], label: 'copies', labelOffset: 30 });
    for (const s of sites) { if (s.kind !== 'MH' || !s.in_segment || s.copies == null) continue;
      el.appendChild(hover(svg('circle', { cx: X(s.pos), cy: Ya(s.copies) + (Math.random() - 0.5) * 10, r: 2.4, fill: s.copies === 1 ? C.jx : C.mat, 'fill-opacity': .8 }), () => F.siteTip(B, s))); }
    if (st) { el.appendChild(text(X((st[0] + st[1]) / 2), 28, `discordant interval ${fmt.range(chrom, st[0], st[1])}: the copies differ at ${(B.checks.find((c) => c.key === 'stretch_differ') || {}).got} of ${(B.checks.find((c) => c.key === 'stretch_n') || {}).got} sites`, { 'text-anchor': 'middle', 'font-size': 10, fill: C.h2 }));
      el.appendChild(text(X(st[0]) - 6, 46, `${(B.checks.find((c) => c.key === 'agree_before') || {}).got} sites concordant`, { 'text-anchor': 'end', 'font-size': 10, fill: C.muted })); el.appendChild(text(X(st[1]) + 6, 46, `${(B.checks.find((c) => c.key === 'agree_after') || {}).got} sites concordant`, { 'font-size': 10, fill: C.muted })); }
    axisX(el, X, 196, AR.ticks(S, E, 8), { fmt: mb, label: `GRCh38 ${chrom} (Mb)` });
    if (!blk) return el;
    // b. against the mother's phased haplotypes
    const pad = 6000, xb0 = 150, Xb = AR.scale(blk.start - pad, blk.end + pad, xb0, x1);
    const el2 = figure(container, 128);
    panelLabel(el2, 8, 14, 'b', `Maternal haplotype carried by each copy, at the mother's phased sites within her ONT phase block ${fmt.range(chrom, blk.start, blk.end)} (${blk.n} sites)`);
    const rowsY = { 2: 40, 1: 62, 0: 84 }, rowLab = { 2: 'both copies: haplotype 1', 1: 'one copy each', 0: 'both copies: haplotype 2' }, rowCol = { 2: C.h1, 1: C.jx, 0: C.h2 };
    if (st) el2.appendChild(svg('rect', { x: Xb(st[0]), y: 30, width: Xb(st[1]) - Xb(st[0]), height: 64, fill: C.band2 }));
    for (const k of [2, 1, 0]) { el2.appendChild(text(xb0 - 8, rowsY[k] + 3.5, rowLab[k], { 'text-anchor': 'end', 'font-size': 10 })); el2.appendChild(line(xb0, rowsY[k], x1, rowsY[k], { stroke: C.grid, 'stroke-dasharray': '2 3' })); }
    for (const b of bounds) if (b) el2.appendChild(line(Xb(b), 30, Xb(b), 94, { stroke: C.muted, 'stroke-dasharray': '3 3' }));
    for (const s of sites) { if (s.hap1_copies == null) continue;
      el2.appendChild(hover(svg('rect', { x: Xb(s.pos) - 1.2, y: rowsY[s.hap1_copies] - 7, width: 2.4, height: 14, fill: rowCol[s.hap1_copies] }), () => `${F.siteTip(B, s)}<br>${rowLab[s.hap1_copies]}`)); }
    axisX(el2, Xb, 100, AR.ticks(blk.start, blk.end, 6), { fmt: mb });
    // c. reads across each transition
    const ends = [['left', bounds[0], bounds[1], 'proximal transition'], ['right', bounds[2], bounds[3], 'distal transition']];
    const grid = h('div', { class: 'cols' }); container.appendChild(grid);
    const keyLabel = (k) => ({ '2>2': 'haplotype 2 throughout', '2>1': 'haplotype 2 → 1', '1>1': 'haplotype 1 throughout', '1>2': 'haplotype 1 → 2', dad: 'paternal (carry paternal-specific alleles)' })[k] || k;
    ends.forEach(([side, lo, hi, title], idx) => {
      if (!lo || !hi) return;
      const list = (B.switch_reads[side] || []).slice();
      const order = ['2>2', '2>1', '1>1', '1>2', 'dad'];
      list.sort((a, b) => order.indexOf(a[1]) - order.indexOf(b[1]) || (readByName.get(a[0]) || {}).start - (readByName.get(b[0]) || {}).start);
      const groups = {}; for (const [n, k] of list) (groups[k] = groups[k] || []).push(n);
      const rowPx = 7, groupGap = 20, Hh = 50 + list.length * rowPx + Object.keys(groups).length * groupGap + 30;
      const box = h('div', {}); const sv = figure(box, Hh, 540);
      const xs0 = 14, xs1 = 528, lo2 = lo - 24000, hi2 = hi + 24000, XX = AR.scale(lo2, hi2, xs0, xs1);
      panelLabel(sv, 6, 14, idx === 0 ? 'c' : 'd', `${title}: ${list.length} HiFi reads spanning ${fmt.range(chrom, lo, hi)}`);
      sv.appendChild(svg('rect', { x: XX(lo), y: 24, width: Math.max(2, XX(hi) - XX(lo)), height: Hh - 56, fill: '#eeeeee' }));
      let y = 34;
      for (const k of order) { const names = groups[k]; if (!names) continue;
        sv.appendChild(text(xs0, y + 4, `${names.length} ${keyLabel(k)}`, { 'font-size': 10, fill: C.muted, 'font-weight': 600 })); y += 12;
        for (const n of names) { const r = readByName.get(n); if (!r) continue;
          const g = svg('g');
          g.appendChild(line(XX(Math.max(lo2, r.start)), y, XX(Math.min(hi2, r.end)), y, { stroke: k === 'dad' ? C.pat : '#bbbbbb', 'stroke-width': k === 'dad' ? 1.1 : 0.8 }));
          for (const [pos, base] of (B.read_alleles[n] || [])) { const s = sites.find((q) => q.pos === pos); if (!s || !s.ref_on_hap || s.block !== blk.id) continue; const hap = base === s.ref ? s.ref_on_hap : 3 - s.ref_on_hap; g.appendChild(svg('circle', { cx: XX(pos), cy: y, r: 2, fill: hap === 1 ? C.h1 : C.h2 })); }
          hover(g, `<b>${n}</b> ${fmt.range(chrom, r.start, r.end)} (${fmt.bp(r.length)}) · ${r.parent === 'unphased' ? 'unphased' : r.parent + ' origin'}<br>informative alleles: maternal ${r.oh_mat + r.own_mat}, paternal ${r.oh_pat + r.own_pat} · ${keyLabel(k)}`);
          sv.appendChild(g); y += rowPx; }
        y += groupGap - 12; }
      axisX(sv, XX, Hh - 24, AR.ticks(lo2, hi2, 5), { fmt: (t) => (t / 1e6).toFixed(3) });
      grid.appendChild(box);
    });
    return el;
  };
})();
