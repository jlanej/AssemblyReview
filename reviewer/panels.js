/* AssemblyReview: the custom panels, drawn as SVG from the case's tables. Each takes a container element and the bundle. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});
  const P = (AR.panels = {});
  const W = 1200; // viewBox width of every panel
  const { svg, h, fmt } = AR;
  const col = AR.colors;
  const parentColor = (p) => (p === 'mother' || p === 'maternal' || p === 'mat' ? col.mat : p === 'father' || p === 'paternal' || p === 'pat' ? col.pat : col.grey);

  function frame(container, height, cls, width) {
    const el = svg('svg', { class: 'panelsvg ' + (cls || ''), viewBox: `0 0 ${width || W} ${height}`, preserveAspectRatio: 'xMidYMin meet' });
    container.appendChild(el);
    return el;
  }
  function text(x, y, s, attrs = {}) { return svg('text', Object.assign({ x, y, text: s }, attrs)); }
  function ribbon(x0a, x1a, ya, x0b, x1b, yb, color, opacity = 0.22) {
    const cy = (ya + yb) / 2;
    const d = `M${x0a},${ya} C${x0a},${cy} ${x0b},${cy} ${x0b},${yb} L${x1b},${yb} C${x1b},${cy} ${x1a},${cy} ${x1a},${ya} Z`;
    return svg('path', { d, fill: color, 'fill-opacity': opacity, stroke: color, 'stroke-opacity': 0.5, 'stroke-width': 0.6 });
  }

  // ---------------------------------------------------------------------------------- 1. the picture
  P.picture = function (container, B, opts = {}) {
    const C = B.case, [S, E] = C.segment, L = E - S + 1, chrom = C.chrom;
    const pieces = B.pieces, reads = AR.rows(B.reads);
    const perPiece = {};
    for (const r of reads) { if (!r.asm_piece) continue; const c = (perPiece[r.asm_piece] = perPiece[r.asm_piece] || { reads: 0, mother: 0, father: 0, unphased: 0, q0: 0 }); c.reads++; c[r.parent]++; if (r.asm_mapq === 0) c.q0++; }
    const FL = 100000, x0 = 250, x1 = 1180, k = (x1 - x0) / (2 * FL + 2 * L);
    const X = (p, copy) => x0 + (FL + (copy - 1) * L + (p - S)) * k;
    const H = 450, bh = 22;
    const Y = { M: 42, mom: 118, U: 186, dad: 256, scaf: 336 };
    const el = frame(container, H, 'picture');
    const hover = (node, html) => { AR.hover(node, html); node.classList.add('hit'); return node; };
    const label = (y, main, sub) => { el.appendChild(text(x0 - 14, y + 13, main, { 'text-anchor': 'end', class: 'lab' })); if (sub) el.appendChild(text(x0 - 14, y + 27, sub, { 'text-anchor': 'end', class: 'mut' })); };
    const bar = (xa, xb, y, fill, txt, opacity = 1, rx = 4) => {
      const g = svg('g');
      g.appendChild(svg('rect', { x: xa, y, width: Math.max(1, xb - xa), height: bh, fill, 'fill-opacity': opacity, rx }));
      if (txt && xb - xa > 30) g.appendChild(text((xa + xb) / 2, y + 15, txt, { 'text-anchor': 'middle', fill: '#fff', 'font-weight': 600 }));
      return g;
    };
    const section = (y, s) => el.appendChild(text(20, y, s, { class: 'mut', 'font-weight': 600, 'letter-spacing': '.08em', 'font-size': 10 }));
    section(22, "HPRC'S ASSEMBLY, AS PUBLISHED");
    section(100, 'THE CHILD, AS ITS READS HOLD IT');
    section(318, "HPRC'S ASSEMBLY, AS PUBLISHED");
    label(Y.M, "HPRC's maternal chromosome", `hap2, ${C.assembly.sequences.hap2.contig}`);
    label(Y.mom, "the mother's chromosome 2", 'in the child: the reads');
    label(Y.U, 'an unplaced contig', `hap2, ${C.assembly.sequences.U.contig}`);
    label(Y.dad, "the father's chromosome 2", 'in the child: the reads');
    label(Y.scaf, "HPRC's paternal chromosome", `hap1, ${C.assembly.sequences.hap1.contig}`);
    const M = pieces.find((p) => p.piece === C.assembly.maternal_chromosome), U = pieces.find((p) => p.piece === 'U');
    const scaf = pieces.filter((p) => p.hap === 'hap1').sort((a, b) => a.qs - b.qs);
    const pieceTip = (p) => {
      const c = perPiece[p.piece] || {};
      return `<b>Piece ${p.piece}</b> · ${p.parent} by the trio's alleles<br>${p.contig}:${fmt.int(p.qs)}-${fmt.int(p.qe)}<br>GRCh38 ${fmt.range(chrom, p.rs, p.re)} (${fmt.kb(p.re - p.rs + 1)})<br>` +
        `mother's / father's own alleles inside the segment: ${p.mom_own} / ${p.dad_own}<br>realigned HiFi reads: ${fmt.int(c.reads)} (mother's ${fmt.int(c.mother)}, father's ${fmt.int(c.father)}, unphased ${fmt.int(c.unphased)}; ${fmt.int(c.q0)} at MAPQ 0)<br><span class="small">click to open it in the assembly browser</span>`;
    };
    // ribbons first (under the bars)
    const ribbons = svg('g');
    const mClipLo = S - FL, mClipHi = E + FL;
    // M -> copy 1 (+ flanks); its right flank to the mother's right flank
    ribbons.appendChild(ribbon(X(S, 1), X(E, 1), Y.M + bh, X(S, 1), X(E, 1), Y.mom, col.mat, 0.18));
    ribbons.appendChild(ribbon(X(E, 1), X(E + FL, 1), Y.M + bh, X(E, 2), X(E + FL, 2), Y.mom, col.grey, 0.15));
    // U -> copy 2
    if (U) ribbons.appendChild(ribbon(X(U.rs, 2), X(U.re, 2), Y.mom + bh, X(U.rs, 2), X(U.re, 2), Y.U, col.mat, 0.18));
    // scaffold pieces -> the copy each holds
    const scafK = k; // same scale
    const A = scaf[0];
    const localOfS = A.local_start + (S - A.rs) * ((A.local_end - A.local_start) / (A.re - A.rs));
    const XS = (local) => Math.max(x0, Math.min(x1, X(S, 1) + (local - localOfS) * scafK));
    for (const p of scaf) {
      const xa = XS(p.local_start), xb = XS(p.local_end);
      const copy = p.parent === 'maternal' ? 2 : 1;
      const yTarget = p.parent === 'maternal' ? Y.mom + bh : Y.dad + bh;
      const lo = Math.max(p.rs, copy === 2 ? S : mClipLo), hi = Math.min(p.re, copy === 2 ? E : mClipHi);
      ribbons.appendChild(ribbon(xa, xb, Y.scaf, X(lo, copy), X(hi, copy), yTarget, p.parent === 'maternal' ? col.mat : col.pat, 0.18));
    }
    el.appendChild(ribbons);
    // bars: M
    el.appendChild(bar(X(S - FL, 1), X(S, 1), Y.M, col.grey, '', 0.5));
    el.appendChild(hover(bar(X(S, 1), X(E, 1), Y.M, col.mat, 'M: the segment once'), () => pieceTip(M)));
    el.appendChild(bar(X(E, 1), X(E + FL, 1), Y.M, col.grey, '', 0.5));
    // the mother's chromosome in the child
    el.appendChild(bar(X(S - FL, 1), X(S, 1), Y.mom, col.grey, '', 0.5));
    el.appendChild(hover(bar(X(S, 1), X(E, 1), Y.mom, col.mat, 'one copy'), `<b>The mother's chromosome in the child</b>: the segment twice in tandem, ${fmt.bp(L)} each, joined at the junction. Which copy is first the reads cannot say.`));
    el.appendChild(hover(bar(X(S, 2), X(E, 2), Y.mom, col.mat, 'the other copy'), `<b>The other maternal copy</b>: HPRC's assembly holds it only in pieces (B, U, D), with ${fmt.bp(B.holes.reduce((s, [a, b]) => s + b - a + 1, 0))} in no piece.`));
    el.appendChild(bar(X(E, 2), X(E + FL, 2), Y.mom, col.grey, '', 0.5));
    for (const [a, b] of B.holes) {
      el.appendChild(hover(svg('rect', { x: X(a, 2), y: Y.mom - 3, width: Math.max(2, X(b, 2) - X(a, 2)), height: bh + 6, fill: col.jx, 'fill-opacity': 0.75, rx: 2 }),
        `<b>Missing from the assembly</b>: ${fmt.range(chrom, a, b)} of the second maternal copy (${fmt.bp(b - a + 1)}) is in no piece.`));
    }
    const jx = X(E, 1);
    el.appendChild(hover(svg('rect', { x: jx - 2, y: Y.mom - 8, width: 4, height: bh + 16, fill: col.jx }), `<b>The junction</b>: ${fmt.pos(chrom, E)} joined to ${fmt.pos(chrom, S)}; ${B.numbers.jr} HiFi reads and the child's short reads hold it.`));
    el.appendChild(text(jx, Y.mom - 12, 'junction', { 'text-anchor': 'middle', fill: col.jx, 'font-weight': 600 }));
    // U
    if (U) el.appendChild(hover(bar(X(U.rs, 2), X(U.re, 2), Y.U, col.mat, 'U'), () => pieceTip(U)));
    if (U) el.appendChild(text(X(U.re, 2) + 8, Y.U + 15, `${U.contig}, ${fmt.bp(U.qe - U.qs + 1)}: no chromosome position`, { class: 'mut' }));
    // the father's chromosome
    el.appendChild(bar(X(S - FL, 1), X(S, 1), Y.dad, col.grey, '', 0.5));
    el.appendChild(hover(bar(X(S, 1), X(E, 1), Y.dad, col.pat, 'his copy'), `<b>The father's chromosome in the child</b>: the segment once (${fmt.bp(L)}). Pieces A, C and E of HPRC's paternal scaffold hold it.`));
    el.appendChild(bar(X(E, 1), X(E + FL, 1), Y.dad, col.grey, '', 0.5));
    // the scaffold
    const gaps = B.gaps_local.filter((g) => g[0] === 'hap1');
    for (const p of scaf) {
      const xa = XS(p.local_start), xb = XS(p.local_end);
      el.appendChild(hover(bar(xa, xb, Y.scaf, p.parent === 'maternal' ? col.mat : col.pat, p.piece), () => pieceTip(p)));
    }
    gaps.forEach((g, i) => {
      const x = XS(g[1]);
      const j = B.joins[i];
      el.appendChild(hover(svg('rect', { x: x - 2, y: Y.scaf - 6, width: 4, height: bh + 12, fill: col.gap }),
        `<b>Gap ${i + 1}</b>: ${g[2] - g[1] + 1} N between ${j ? j.left : '?'} and ${j ? j.right : '?'}${j ? `<br>on GRCh38 it jumps back ${fmt.bp(j.back)}: ${fmt.pos(chrom, j.p1)} to ${fmt.pos(chrom, j.p2)}<br>HiFi reads carrying the join: ${j.carry}; straight through its points: ${j.through1} / ${j.through2}` : ''}`));
    });
    // callouts
    const callout = (x, y, n, tip) => { const g = svg('g', { class: 'hit' }); g.appendChild(svg('circle', { cx: x, cy: y, r: 9, fill: col.jx })); g.appendChild(text(x, y + 4, String(n), { 'text-anchor': 'middle', fill: '#fff', 'font-weight': 700, 'font-size': 11 })); el.appendChild(hover(g, tip)); };
    callout(X((S + E) / 2, 1), Y.M - 14, 1, '<b>1</b> The two maternal copies collapsed into one, M: they are identical outside one 114-kb stretch, and no read spans a 1.09-Mb copy.');
    const Bp = scaf.find((p) => p.piece === 'B');
    if (Bp) callout((XS(Bp.local_start) + X(Bp.rs, 2)) / 2, (Y.scaf + Y.mom + bh) / 2, 2, "<b>2</b> Pieces of the other copy, B and D, were put in the father's chromosome; U was left without a place.");
    if (gaps[0]) callout(XS(gaps[0][1]), Y.scaf + bh + 18, 3, '<b>3</b> Joined across 100-N gaps that no read crosses.');
    if (B.holes[1]) callout(X((B.holes[1][0] + B.holes[1][1]) / 2, 2), Y.mom - 14, 4, `<b>4</b> ${fmt.kb(B.holes.reduce((s, [a, b]) => s + b - a + 1, 0))} of the other copy, and the junction, are in no piece.`);
    const Cp = scaf.find((p) => p.piece === 'C'), Ap = scaf[0];
    if (Cp && Ap) callout((XS(Ap.local_end) + XS(Cp.local_start)) / 2, Y.dad + bh + 22, 5, "<b>5</b> The father's pieces overlap: A and C over 20.5 kb, C and E over 17.8 kb, so two stretches of his sequence are there twice.");
    el.appendChild(text(20, H - 34, '1 The two maternal copies collapsed into one, M.  2 Pieces of the other copy, B and D, put in the father\'s chromosome; U left without a place.  3 Joined across 100-N gaps no read crosses.', { class: 'mut' }));
    el.appendChild(text(20, H - 18, '4 187 kb of the other copy, and the junction, in no piece.  5 The father\'s pieces overlap.   To scale, the flanks (100 kb) faint; the paternal chromosome in its own order, a black mark at each gap; colours by each piece\'s alleles.', { class: 'mut' }));
    // click handlers
    if (opts.onPiece) for (const g of el.querySelectorAll('g.hit')) { /* callouts: nothing */ }
    el.querySelectorAll('g').forEach((g) => { });
    const clickable = [...scaf, M, U].filter(Boolean);
    let i = 0;
    for (const node of el.querySelectorAll('g.hit')) {
      const t = node.querySelector('text');
      const letter = t ? t.textContent.split(':')[0].trim() : '';
      const p = clickable.find((q) => q.piece === letter);
      if (p && opts.onPiece) node.addEventListener('click', () => opts.onPiece(p));
      i++;
    }
    return el;
  };

  // ---------------------------------------------------------------------------------- 2. the evidence board on GRCh38
  P.overview = function (container, B, opts = {}) {
    const C = B.case, [S, E] = C.segment, chrom = C.chrom, [W0, W1] = C.windows.WIN;
    const depth = AR.rows(B.depth_grch38), sites = AR.rows(B.sites), pieces = B.pieces, calls = B.calls, reads = AR.rows(B.reads);
    const Mletter = C.assembly.maternal_chromosome;
    const jr = reads.filter((r) => r.junction);
    const jrMat = AR.count(jr, (r) => r.jr_parent === 'mother'), jrPat = AR.count(jr, (r) => r.jr_parent === 'father');
    const thrS = AR.count(reads, (r) => r.through_S), thrE = AR.count(reads, (r) => r.through_E);
    const perPiece = {};
    for (const r of reads) { if (!r.asm_piece) continue; const c = (perPiece[r.asm_piece] = perPiece[r.asm_piece] || { reads: 0, mother: 0, father: 0, unphased: 0, q0: 0 }); c.reads++; c[r.parent]++; if (r.asm_mapq === 0) c.q0++; }
    const vsM = {}; for (const a of B.asm_vs_asm || []) vsM[a.piece] = a;
    const x0 = 205, x1 = 1180;
    const tracks = [
      { id: 'ruler', h: 40 },
      { id: 'sr', h: 108, label: 'SHORT READS', sub: 'copies: child (black); father, mother; the assembly (pink)' },
      { id: 'hifi', h: 108, label: 'LONG READS vs GRCh38', sub: 'HiFi copies: all (dashed), by parent (dots)' },
      { id: 'oh', h: 62, label: "the mother's allele", sub: 'parent-opposite sites, short reads' },
      { id: 'mh', h: 52, label: 'her own allele, HiFi', sub: 'mother-only sites: 0, 1/3, 2/3' },
      { id: 'asmref', h: 160, label: 'ASSEMBLY vs GRCh38', sub: "HPRC's alignment: contig coordinate against GRCh38" },
      { id: 'asm', h: 80, label: 'pieces, whose by alleles', sub: 'A, C, E his; B, D, M, U hers' },
      { id: 'asmreads', h: 62, label: 'ASSEMBLY vs LONG READS', sub: 'whose reads each piece takes, realigned to it' },
      { id: 'asmasm', h: 56, label: 'ASSEMBLY vs ASSEMBLY', sub: 'each piece aligned directly to M' },
      { id: 'calls', h: 60, label: 'CALL SETS', sub: 'records over the segment' },
    ];
    let H = 10; const ty = {}; for (const t of tracks) { ty[t.id] = H; H += t.h; } H += 10;
    const el = frame(container, H, 'overview');
    let d0 = W0, d1 = W1;
    const body = svg('g'); el.appendChild(body);
    const brush = svg('rect', { class: 'brush', y: 0, height: H, x: 0, width: 0, visibility: 'hidden' }); el.appendChild(brush);
    const cross = svg('line', { class: 'cross', y1: 0, y2: H, x1: 0, x2: 0, visibility: 'hidden' }); el.appendChild(cross);
    const X = () => AR.scale(d0, d1, x0, x1);
    const clipId = `clip-${Math.random().toString(36).slice(2)}`;
    el.appendChild(svg('defs', {}, svg('clipPath', { id: clipId }, svg('rect', { x: x0, y: 0, width: x1 - x0, height: H }))));
    const rowOf = { A: 0, C: 0, E: 0, B: 1, D: 1, M: 2, U: 3 };
    const siteTip = (s) => `<b>${fmt.pos(chrom, s.pos)}</b> ${s.ref}/${s.alt} · ${({ OH: 'parent-opposite', MH: 'mother-only', FH: 'father-only', BH: 'both heterozygous' })[s.kind]}${s.in_segment ? ' · inside' : ' · outside'}${s.in_stretch ? ', in the stretch' : ''}<br>` +
      `child ${s.child_gt} (${s.child_ref}/${s.child_alt}), father ${s.father_gt} (${s.father_ref}/${s.father_alt}), mother ${s.mother_gt} (${s.mother_ref}/${s.mother_alt})<br>HiFi ${s.hifi_ref}/${s.hifi_alt}` +
      (s.share_sr != null ? `<br>the mother's allele ${s.mother_allele}: ${fmt.f(s.share_sr)} of the child's short reads, ${fmt.f(s.share_hifi)} of its HiFi reads` : '') +
      (s.copies != null ? `<br>copies carrying her own allele: ${s.copies} of 3 → the maternal copies ${s.copies === 1 ? 'differ' : 'agree'}` : '') +
      (s.ref_on_hap ? `<br>her phase: ${s.ref} on her haplotype ${s.ref_on_hap} (block ${s.block})` : '') +
      ['A', 'B', 'C', 'D', 'E', 'M', 'U'].filter((l) => s['base_' + l]).map((l) => ` ${l}:${s['base_' + l]}`).join('').replace(/^/, '<br>pieces:');
    const pieceTip = (p) => { const c = perPiece[p.piece] || {}; const a = vsM[p.piece];
      return `<b>Piece ${p.piece}</b> (${p.parent} by alleles) ${p.contig}:${fmt.int(p.qs)}-${fmt.int(p.qe)}<br>GRCh38 ${fmt.range(chrom, p.rs, p.re)}<br>own alleles inside the segment: mother's ${p.mom_own}, father's ${p.dad_own}` +
        `<br>realigned reads: ${fmt.int(c.reads)} (mother's ${fmt.int(c.mother)}, father's ${fmt.int(c.father)}, unphased ${fmt.int(c.unphased)}; ${fmt.int(c.q0)} at MAPQ 0)` +
        (a ? `<br>against M, aligned directly: ${a.identity === 1 ? 'identical' : fmt.pct(a.identity, 3) + ' identity'} over ${fmt.bp(a.blen)}` : ''); };
    const arc = (g, xa, xb, yb, lift, color, label, tip) => {
      const d = `M${xa},${yb} Q${(xa + xb) / 2},${yb - lift} ${xb},${yb}`;
      const p = svg('path', { d, fill: 'none', stroke: color, 'stroke-width': 1.6, 'stroke-dasharray': '5 3' }); if (tip) AR.hover(p, tip); g.appendChild(p);
      if (label) g.appendChild(text((xa + xb) / 2, yb - lift / 2 - 2, label, { 'text-anchor': 'middle', fill: color, 'font-weight': 600, 'font-size': 10.5 }));
    };

    function draw() {
      AR.clear(body);
      const x = X();
      const inView = (a, b) => b >= d0 && a <= d1;
      for (const t of tracks) {
        if (!t.label) continue;
        const caps = t.label === t.label.toUpperCase();
        body.appendChild(text(x0 - 10, ty[t.id] + 16, t.label, { 'text-anchor': 'end', class: 'lab', 'letter-spacing': caps ? '.06em' : null, 'font-size': caps ? 10.5 : 11.5 }));
        if (t.sub) body.appendChild(text(x0 - 10, ty[t.id] + 30, t.sub, { 'text-anchor': 'end', class: 'mut', 'font-size': 10 }));
        body.appendChild(svg('line', { x1: x0, x2: x1, y1: ty[t.id] + t.h - 0.5, y2: ty[t.id] + t.h - 0.5, stroke: 'var(--line)' }));
      }
      const g = svg('g', { 'clip-path': `url(#${clipId})` }); body.appendChild(g);
      // ruler
      const ry = ty.ruler;
      for (const tk of AR.ticks(d0, d1, 8)) { const px = x(tk); g.appendChild(svg('line', { x1: px, x2: px, y1: ry + 24, y2: ry + 30, stroke: 'var(--muted)' })); g.appendChild(text(px, ry + 20, (tk / 1e6).toFixed(d1 - d0 < 200000 ? 3 : 2) + ' Mb', { 'text-anchor': 'middle', class: 'mut' })); }
      g.appendChild(svg('line', { x1: x0, x2: x1, y1: ry + 30, y2: ry + 30, stroke: 'var(--muted)' }));
      g.appendChild(AR.hover(svg('rect', { x: x(S), y: ry + 32, width: Math.max(2, x(E) - x(S)), height: 6, fill: col.jx, 'fill-opacity': .8 }), `<b>The segment</b> ${fmt.range(chrom, S, E)}, ${fmt.bp(E - S + 1)}`));
      if (B.stretch) g.appendChild(AR.hover(svg('rect', { x: x(B.stretch[0]), y: ry + 32, width: Math.max(2, x(B.stretch[1]) - x(B.stretch[0])), height: 6, fill: col.h2 }), `<b>The stretch</b> ${fmt.range(chrom, B.stretch[0], B.stretch[1])}: where the two maternal copies differ`));
      for (const [a, b] of B.holes) g.appendChild(svg('rect', { x: x(a), y: ry + 34, width: Math.max(1, x(b) - x(a)), height: 2, fill: '#fff' }));
      const vis = depth.filter((d) => inView(d.bin0 + 1, d.bin1));
      // short reads: copies, the assembly's copies, the pairs' arc
      const sy = ty.sr, sbase = sy + tracks[1].h - 12, Ys = AR.scale(0, 4.2, sbase, sy + 34);
      for (const v of [1, 2, 3, 4]) { g.appendChild(svg('line', { x1: x0, x2: x1, y1: Ys(v), y2: Ys(v), stroke: 'var(--line)', 'stroke-dasharray': '2 3' })); body.appendChild(text(x0 - 4, Ys(v) + 3, String(v), { 'text-anchor': 'end', class: 'mut', 'font-size': 9 })); }
      const line = (key, Y, stroke, dash, width = 1.4) => {
        let dd = ''; for (const d of vis) { if (d[key] == null) continue; const px = x((d.bin0 + d.bin1) / 2); dd += (dd ? 'L' : 'M') + `${px.toFixed(1)},${Y(Math.min(4.2, d[key])).toFixed(1)}`; }
        g.appendChild(svg('path', { d: dd, fill: 'none', stroke, 'stroke-width': width, 'stroke-dasharray': dash || null }));
      };
      let sd = '';
      for (const d of vis) { if (d.asm_copies == null) continue; sd += `${sd ? 'L' : 'M'}${x(d.bin0 + 1).toFixed(1)},${Ys(d.asm_copies).toFixed(1)} L${x(d.bin1).toFixed(1)},${Ys(d.asm_copies).toFixed(1)}`; }
      g.appendChild(svg('path', { d: sd, fill: 'none', stroke: col.jx, 'stroke-width': 1.6, 'stroke-opacity': .7 }));
      line('cp_father', Ys, col.pat, null, 1); line('cp_mother', Ys, col.mat, null, 1); line('cp_child', Ys, 'var(--ink)', null, 1.6);
      const srC = B.sr.child, srP = B.sr.father.pairs + B.sr.mother.pairs + B.sr.father.split + B.sr.mother.split;
      arc(g, x(E), x(S), sy + 30, 22, 'var(--ink)', `${srC.pairs} pairs and ${srC.split} split reads join the end to the start in the child; ${srP} in the parents`,
        `<b>The junction in the short reads</b>: ${srC.pairs} pairs in a tandem duplication's orientation and ${srC.split} split reads in ${C.trio.child}; none in either parent. Clipped exactly at the start / end: ${srC.clip_s} / ${srC.clip_e}.`);
      // HiFi: copies by parent, the junction arc, straight-through counts
      const hy = ty.hifi, hbase = hy + tracks[2].h - 12, Yh = AR.scale(0, 4.2, hbase, hy + 34);
      for (const v of [1, 2, 3, 4]) { g.appendChild(svg('line', { x1: x0, x2: x1, y1: Yh(v), y2: Yh(v), stroke: 'var(--line)', 'stroke-dasharray': '2 3' })); body.appendChild(text(x0 - 4, Yh(v) + 3, String(v), { 'text-anchor': 'end', class: 'mut', 'font-size': 9 })); }
      line('cp_hifi', Yh, col.h1, '4 3', 1.4);
      for (const d of vis) if (d.cp_hifi_mat != null) { g.appendChild(svg('circle', { cx: x((d.bin0 + d.bin1) / 2), cy: Yh(Math.min(4.2, d.cp_hifi_mat)), r: 2.6, fill: col.mat })); g.appendChild(svg('circle', { cx: x((d.bin0 + d.bin1) / 2), cy: Yh(Math.min(4.2, d.cp_hifi_pat)), r: 2.6, fill: col.pat })); }
      arc(g, x(E), x(S), hy + 30, 22, col.jx, `${jr.length} HiFi molecules read the junction: ${jrMat} the mother's, ${jrPat} the father's, ${jr.length - jrMat - jrPat} unphased`,
        `<b>The junction reads</b>: ${jr.length} reads align up to ${fmt.pos(chrom, E)} and go on from ${fmt.pos(chrom, S)}, both sides at MAPQ 60. Straight through the start: ${thrS}; through the end: ${thrE}: one junction against two chromosomes' worth of straight reads.`);
      for (const [p, n, lab] of [[S, thrS, 'start'], [E, thrE, 'end']]) g.appendChild(AR.hover(text(x(p), hbase - 2, `${n} straight through the ${lab}`, { 'text-anchor': p === S ? 'start' : 'end', class: 'mut', 'font-size': 9.5, dx: p === S ? 4 : -4 }), `${n} HiFi reads run straight through ${fmt.pos(chrom, p)}, 1 kb either side${p === S ? ` (${AR.count(reads, (r) => r.through_S && r.parent === 'mother')} the mother's, ${AR.count(reads, (r) => r.through_S && r.parent === 'father')} the father's)` : ''}.`));
      for (const d of vis) {
        const r = svg('rect', { x: x(d.bin0 + 1), y: sy + 34, width: Math.max(1, x(d.bin1) - x(d.bin0 + 1)), height: tracks[1].h + tracks[2].h - 34, fill: 'transparent' });
        AR.hover(r, () => `<b>${fmt.range(chrom, d.bin0 + 1, d.bin1)}</b>${d.inside ? ' · inside the segment' : d.flank ? ' · flank' : ''}<br>short reads, copies: child ${fmt.f(d.cp_child)} (${fmt.f(d.cp_rel)} against the parents), father ${fmt.f(d.cp_father)}, mother ${fmt.f(d.cp_mother)}<br>` +
          `HiFi copies: ${fmt.f(d.cp_hifi)}${d.cp_hifi_mat != null ? ` · maternal ${fmt.f(d.cp_hifi_mat)}, paternal ${fmt.f(d.cp_hifi_pat)} (${d.phased_reads_mat} + ${d.phased_reads_pat} phased reads)` : ''}<br>the assembly holds ${d.asm_copies} ${d.asm_copies === 1 ? 'copy' : 'copies'} here`);
        g.appendChild(r);
      }
      // the mother's allele at parent-opposite sites
      const oy = ty.oh, oh = tracks[3].h - 10, Yo = AR.scale(0, 1, oy + oh, oy + 6);
      for (const v of [0.5, 2 / 3]) g.appendChild(svg('line', { x1: x0, x2: x1, y1: Yo(v), y2: Yo(v), stroke: 'var(--line)', 'stroke-dasharray': '2 3' }));
      body.appendChild(text(x0 - 4, Yo(0.5) + 3, '1/2', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 })); body.appendChild(text(x0 - 4, Yo(2 / 3) + 3, '2/3', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 }));
      for (const s of sites) if (inView(s.pos, s.pos) && s.kind === 'OH' && s.share_sr != null) g.appendChild(AR.hover(svg('circle', { cx: x(s.pos), cy: Yo(s.share_sr), r: 2.4, fill: s.in_segment ? col.mat : col.grey, 'fill-opacity': .85 }), () => siteTip(s)));
      // her own allele at mother-only sites
      const my = ty.mh, mh = tracks[4].h - 10, Ym = AR.scale(0, 1, my + mh, my + 6);
      for (const v of [1 / 3, 2 / 3]) g.appendChild(svg('line', { x1: x0, x2: x1, y1: Ym(v), y2: Ym(v), stroke: 'var(--line)', 'stroke-dasharray': '2 3' }));
      body.appendChild(text(x0 - 4, Ym(1 / 3) + 3, '1/3', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 })); body.appendChild(text(x0 - 4, Ym(2 / 3) + 3, '2/3', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 }));
      for (const s of sites) if (inView(s.pos, s.pos) && s.kind === 'MH' && s.share_hifi != null && s.in_segment) g.appendChild(AR.hover(svg('circle', { cx: x(s.pos), cy: Ym(s.share_hifi), r: 2.2, fill: s.copies === 1 ? col.jx : s.copies == null ? col.grey : col.mat, 'fill-opacity': .85 }), () => siteTip(s)));
      // the assembly against GRCh38: contig coordinate against reference coordinate, each piece a diagonal
      const ay0 = ty.asmref, aH = tracks[5].h;
      const hap1 = pieces.filter((p) => p.hap === 'hap1').sort((a, b) => a.qs - b.qs), hap2 = pieces.filter((p) => p.hap !== 'hap1');
      const b1top = ay0 + 8, b1bot = ay0 + aH - 56, b2top = ay0 + aH - 48, b2bot = ay0 + aH - 6;
      const qlo = Math.min(...hap1.map((p) => p.qs)), qhi = Math.max(...hap1.map((p) => p.qe));
      const Y1 = AR.scale(qlo, qhi, b1bot, b1top);
      g.appendChild(svg('rect', { x: x0, y: b1top, width: x1 - x0, height: b1bot - b1top, fill: 'var(--shade)', 'fill-opacity': .5 }));
      g.appendChild(svg('rect', { x: x0, y: b2top, width: x1 - x0, height: b2bot - b2top, fill: 'var(--shade)', 'fill-opacity': .5 }));
      body.appendChild(text(x1 + 4, b1top + 10, 'hap1', { class: 'mut', 'font-size': 9.5 })); body.appendChild(text(x1 + 4, b2top + 10, 'hap2', { class: 'mut', 'font-size': 9.5 }));
      hap1.forEach((p, i) => {
        const ln = svg('line', { x1: x(p.rs), y1: Y1(p.qs), x2: x(p.re), y2: Y1(p.qe), stroke: parentColor(p.parent), 'stroke-width': 3, 'stroke-linecap': 'round' });
        AR.hover(ln, () => pieceTip(p)); ln.classList.add('hit'); ln.addEventListener('click', (e) => { e.stopPropagation(); if (opts.onPiece) opts.onPiece(p); }); g.appendChild(ln);
        g.appendChild(text((x(p.rs) + x(p.re)) / 2, (Y1(p.qs) + Y1(p.qe)) / 2 - 6, p.piece, { 'text-anchor': 'middle', 'font-weight': 700, 'font-size': 10, fill: parentColor(p.parent) }));
        const q = hap1[i + 1];
        if (q) {
          const j = B.joins[i];
          const cn = svg('line', { x1: x(p.re), y1: Y1(p.qe), x2: x(q.rs), y2: Y1(q.qs), stroke: 'var(--ink)', 'stroke-width': 1, 'stroke-dasharray': '3 3' });
          AR.hover(cn, j ? `<b>Gap ${i + 1}</b> (${j.left} to ${j.right}): 100 N on the contig; on GRCh38 the scaffold jumps back ${fmt.bp(j.back)}, from ${fmt.pos(chrom, j.p1)} to ${fmt.pos(chrom, j.p2)}.<br>HiFi reads carrying the join: ${j.carry}; straight through its points: ${j.through1} / ${j.through2}; realigned reads across the gap: ${j.realign.cross}` : ''); g.appendChild(cn);
          g.appendChild(svg('circle', { cx: (x(p.re) + x(q.rs)) / 2, cy: (Y1(p.qe) + Y1(q.qs)) / 2, r: 3, fill: 'var(--ink)' }));
        }
      });
      for (const p of hap2) {
        const isM = p.piece === Mletter;
        const Y2 = isM ? AR.scale(p.qs, p.qe, b2bot, b2top) : AR.scale(p.qs, p.qe, b2bot - 2, b2bot - 14);
        const ln = svg('line', { x1: x(p.rs), y1: Y2(p.qs), x2: x(p.re), y2: Y2(p.qe), stroke: parentColor(p.parent), 'stroke-width': isM ? 3 : 2.5, 'stroke-dasharray': isM ? null : '4 3', 'stroke-linecap': 'round' });
        AR.hover(ln, () => pieceTip(p)); ln.classList.add('hit'); ln.addEventListener('click', (e) => { e.stopPropagation(); if (opts.onPiece) opts.onPiece(p); }); g.appendChild(ln);
        g.appendChild(text(isM ? x(Math.max(p.rs, d0)) + 6 : (x(p.rs) + x(p.re)) / 2, isM ? b2top + 11 : b2bot - 16, isM ? `M: the segment once, straight through` : 'U, unplaced', { 'font-weight': 700, 'font-size': 10, fill: parentColor(p.parent), 'text-anchor': isM ? 'start' : 'middle' }));
      }
      // pieces, whose by alleles
      const py = ty.asm;
      for (const p of pieces) {
        const r = rowOf[p.piece] ?? 0, yy = py + 6 + r * 17;
        if (!inView(p.rs, p.re)) continue;
        const rect = svg('rect', { x: x(Math.max(p.rs, d0 - 1e6)), y: yy, width: Math.max(1.5, x(Math.min(p.re, d1 + 1e6)) - x(Math.max(p.rs, d0 - 1e6))), height: 13, rx: 2, fill: parentColor(p.parent) });
        AR.hover(rect, () => pieceTip(p)); rect.classList.add('hit'); rect.addEventListener('click', (e) => { e.stopPropagation(); if (opts.onPiece) opts.onPiece(p); });
        g.appendChild(rect);
        if (x(p.re) - x(p.rs) > 14) g.appendChild(text((x(Math.max(p.rs, d0)) + x(Math.min(p.re, d1))) / 2, yy + 10, p.piece, { 'text-anchor': 'middle', fill: '#fff', 'font-weight': 700, 'font-size': 10 }));
      }
      for (const [a, b] of B.holes) if (inView(a, b)) g.appendChild(AR.hover(svg('rect', { x: x(a), y: py + 6 + 17, width: Math.max(1.5, x(b) - x(a)), height: 13, fill: col.jx, 'fill-opacity': .5 }), `<b>Missing</b>: no piece of the second maternal copy holds ${fmt.range(chrom, a, b)} (${fmt.bp(b - a + 1)})`));
      for (const f of B.qc_placed) if (inView(f.rs, f.re)) g.appendChild(AR.hover(svg('rect', { x: x(f.rs) - 1, y: py + 6 + (rowOf[f.piece] ?? 0) * 17 + 14, width: Math.max(2, x(f.re) - x(f.rs) + 2), height: 3, fill: f.label === 'COLLAPSE' ? col.bad : '#f59e0b' }), `<b>NucFlag ${f.label}</b> on ${f.piece}: ${f.contig}:${fmt.int(f.qs)}-${fmt.int(f.qe)} = ${fmt.range(chrom, f.rs, f.re)}`));
      // the assembly against the long reads: whose reads each piece takes
      const ry2 = ty.asmreads;
      for (const p of pieces) {
        const c = perPiece[p.piece]; if (!c || !inView(p.rs, p.re)) continue;
        const r = rowOf[p.piece] ?? 0, yy = ry2 + 6 + r * 12;
        const xa = x(Math.max(p.rs, d0)), xb = x(Math.min(p.re, d1)), w = Math.max(2, xb - xa);
        const parts = [['mother', col.mat], ['father', col.pat], ['unphased', col.grey]];
        let cx = xa;
        const grp = svg('g');
        for (const [k, cc] of parts) { const ww = w * (c[k] / c.reads); if (ww > 0) grp.appendChild(svg('rect', { x: cx, y: yy, width: ww, height: 9, fill: cc })); cx += ww; }
        if (w > 60) grp.appendChild(text(xa + 3, yy + 8, `${p.piece}: ${fmt.int(c.mother)} / ${fmt.int(c.father)}`, { fill: '#fff', 'font-weight': 700, 'font-size': 8.5 }));
        AR.hover(grp, () => `<b>${p.piece}</b>: ${fmt.int(c.reads)} realigned HiFi reads: the mother's ${fmt.int(c.mother)}, the father's ${fmt.int(c.father)}, unphased ${fmt.int(c.unphased)}; ${fmt.int(c.q0)} at MAPQ 0${p.piece === Mletter ? `<br>the ${jr.length} junction reads: ${AR.count(jr, (r) => r.asm_split_end && r.asm_split_start)} split at both of M's segment ends, ${AR.count(jr, (r) => r.asm_through_M)} run through` : ''}`);
        g.appendChild(grp);
      }
      // the assembly against itself: each piece against M
      const zy = ty.asmasm;
      for (const p of pieces) {
        const a = vsM[p.piece]; if (!a || !inView(p.rs, p.re)) continue;
        const r = rowOf[p.piece] ?? 0, yy = zy + 6 + Math.min(r, 2) * 15 + (p.piece === 'U' ? 15 : 0);
        const xa = x(Math.max(p.rs, d0)), xb = x(Math.min(p.re, d1));
        const same = a.identity >= 0.99999;
        const rect = svg('rect', { x: xa, y: yy, width: Math.max(2, xb - xa), height: 11, rx: 2, fill: same ? col.mat : col.pat, 'fill-opacity': same ? 1 : .55 });
        AR.hover(rect, `<b>${p.piece} against M</b>, aligned directly (minimap2 asm20): ${same ? 'identical, base for base' : fmt.pct(a.identity, 3) + ' identity'} over ${fmt.bp(a.blen)}<br>${p.piece} ${a.q}:${fmt.int(a.qs)}-${fmt.int(a.qe)} ↔ M ${a.t}:${fmt.int(a.ts)}-${fmt.int(a.te)}`);
        g.appendChild(rect);
        if (xb - xa > 70) g.appendChild(text(xa + 3, yy + 9, same ? `${p.piece} = M, ${fmt.kb(a.blen)} identical` : `${p.piece}: ${fmt.pct(a.identity, 2)}`, { fill: '#fff', 'font-weight': 700, 'font-size': 8.5 }));
      }
      // calls
      const cy = ty.calls;
      const rowc = { 'the junction': 0, 'the scaffold': 1, other: 2 };
      for (const c of calls) {
        const a = c.pos, b = Math.max(c.end, c.pos + 1);
        if (!inView(a, b)) continue;
        const yy = cy + 8 + (rowc[c.reflects] ?? 2) * 16;
        const color = c.reflects === 'the junction' ? col.jx : c.reflects === 'the scaffold' ? col.pat : col.grey;
        const r = svg('rect', { x: x(a), y: yy, width: Math.max(2, x(b) - x(a)), height: 10, rx: 1, fill: color, 'fill-opacity': .7 });
        AR.hover(r, `<b>${c.set}</b> (${c.made_from}) ${c.type} ${fmt.range(chrom, c.pos, c.end)} · ${fmt.bp(c.len)} · GT ${c.gt || '–'}<br>reflects ${c.reflects}`);
        g.appendChild(r);
      }
      body.appendChild(text(x0 - 4, cy + 16, 'junction', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 })); body.appendChild(text(x0 - 4, cy + 32, 'scaffold', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 }));
      body.appendChild(text(x0 - 4, cy + 48, 'other', { 'text-anchor': 'end', class: 'mut', 'font-size': 9 }));
      for (const p of [S, E]) g.appendChild(svg('line', { x1: x(p), x2: x(p), y1: ty.sr + 34, y2: H - 10, stroke: col.jx, 'stroke-width': .8, 'stroke-dasharray': '3 3', 'pointer-events': 'none' }));
    }
    draw();
    let drag = null;
    const pt = (e) => { const r = el.getBoundingClientRect(); return ((e.clientX - r.left) / r.width) * W; };
    el.addEventListener('mousedown', (e) => { const px = pt(e); if (px < x0 || px > x1) return; drag = { a: px, moved: false }; brush.setAttribute('x', px); brush.setAttribute('width', 0); brush.setAttribute('visibility', 'visible'); e.preventDefault(); });
    el.addEventListener('mousemove', (e) => {
      const px = Math.min(x1, Math.max(x0, pt(e)));
      cross.setAttribute('x1', px); cross.setAttribute('x2', px); cross.setAttribute('visibility', 'visible');
      if (!drag) return; drag.moved = drag.moved || Math.abs(px - drag.a) > 3;
      brush.setAttribute('x', Math.min(drag.a, px)); brush.setAttribute('width', Math.abs(px - drag.a));
    });
    el.addEventListener('mouseleave', () => cross.setAttribute('visibility', 'hidden'));
    window.addEventListener('mouseup', (e) => {
      if (!drag) return; const px = Math.min(x1, Math.max(x0, pt(e))); const x = X();
      brush.setAttribute('visibility', 'hidden');
      if (drag.moved) { const a = x.invert(Math.min(drag.a, px)), b = x.invert(Math.max(drag.a, px)); if (b - a > 500) { d0 = a; d1 = b; draw(); if (opts.onRange) opts.onRange(d0, d1); } }
      else if (opts.onLocus) opts.onLocus(Math.round(x.invert(px)));
      drag = null;
    });
    el.addEventListener('dblclick', () => { d0 = W0; d1 = W1; draw(); if (opts.onRange) opts.onRange(d0, d1); });
    return { el, zoom(a, b) { d0 = Math.max(W0, a); d1 = Math.min(W1, b); draw(); if (opts.onRange) opts.onRange(d0, d1); }, domain: () => [d0, d1] };
  };

  // ---------------------------------------------------------------------------------- 2b. the assembly against itself, in its own coordinates
  P.asmDotplot = function (container, B) {
    const C = B.case, al = B.asm_vs_asm || [], len = B.local_len, Mletter = C.assembly.maternal_chromosome;
    if (!al.length) return;
    const x0 = 132, x1 = 1150, top = 30, h1 = 330, gap = 26, hU = 60, H = top + h1 + gap + hU + 44;
    const el = frame(container, H, 'asmdot');
    const X = AR.scale(1, len.hap2, x0, x1);
    const Y1 = AR.scale(1, len.hap1, top + h1, top);
    const YU = AR.scale(1, len.U || 1, top + h1 + gap + hU, top + h1 + gap);
    el.appendChild(svg('rect', { x: x0, y: top, width: x1 - x0, height: h1, fill: 'var(--shade)', 'fill-opacity': .5 }));
    if (len.U) el.appendChild(svg('rect', { x: x0, y: top + h1 + gap, width: x1 - x0, height: hU, fill: 'var(--shade)', 'fill-opacity': .5 }));
    el.appendChild(text(x0, 16, "hap1 (y) and U against hap2 (x), each gap-free piece aligned on its own with minimap2 asm20: amber = identical to M, blue = the father's sequence", { class: 'lab' }));
    // pieces on the y axis, gaps across
    for (const [letter, [sname, a, b]] of Object.entries(B.piece_local)) {
      const p = B.pieces.find((q) => q.piece === letter);
      if (sname === 'hap1') {
        el.appendChild(AR.hover(svg('rect', { x: x0 - 16, y: Y1(b), width: 10, height: Y1(a) - Y1(b), fill: parentColor(p.parent), rx: 2 }), `<b>${letter}</b> hap1:${fmt.int(a)}-${fmt.int(b)} (${p.parent} by alleles)`));
        el.appendChild(text(x0 - 20, (Y1(a) + Y1(b)) / 2 + 4, letter, { 'text-anchor': 'end', 'font-weight': 700 }));
      }
    }
    B.gaps_local.forEach((gp, i) => { if (gp[0] !== 'hap1') return; el.appendChild(svg('line', { x1: x0, x2: x1, y1: Y1(gp[1]), y2: Y1(gp[1]), stroke: 'var(--ink)', 'stroke-dasharray': '3 3', 'stroke-opacity': .6 })); el.appendChild(text(x1 + 4, Y1(gp[1]) + 3, `gap ${i + 1}`, { class: 'mut', 'font-size': 9.5 })); });
    // M's segment on the x axis
    if (B.m_segment) {
      const [ms, me] = B.m_segment;
      for (const p of [ms, me]) el.appendChild(svg('line', { x1: X(p), x2: X(p), y1: top, y2: top + h1 + gap + hU, stroke: col.jx, 'stroke-dasharray': '4 3' }));
      el.appendChild(AR.hover(svg('rect', { x: X(ms), y: H - 30, width: X(me) - X(ms), height: 6, fill: col.jx }), `<b>M's one copy of the segment</b>: hap2:${fmt.int(ms)}-${fmt.int(me)}`));
      el.appendChild(text((X(ms) + X(me)) / 2, H - 34, "M: the segment, once", { 'text-anchor': 'middle', fill: col.jx, 'font-weight': 600, 'font-size': 10 }));
    }
    for (const t of AR.ticks(1, len.hap2, 8)) { el.appendChild(text(X(t), H - 10, fmt.int(t), { 'text-anchor': 'middle', class: 'mut' })); }
    for (const t of AR.ticks(1, len.hap1, 6)) { el.appendChild(text(x0 - 40, Y1(t) + 3, fmt.int(t), { 'text-anchor': 'end', class: 'mut', 'font-size': 9.5 })); }
    el.appendChild(text(x1 + 4, top + 12, 'hap1', { class: 'mut' })); if (len.U) el.appendChild(text(x1 + 4, top + h1 + gap + 12, 'U', { class: 'mut' }));
    for (const a of al) {
      const Y = a.q === 'U' ? YU : Y1;
      const same = a.identity >= 0.99999;
      const ln = svg('line', { x1: X(a.ts), y1: Y(a.qs), x2: X(a.te), y2: Y(a.qe), stroke: same ? col.mat : col.pat, 'stroke-width': same ? 4 : 3, 'stroke-linecap': 'round' });
      AR.hover(ln, `<b>${a.piece} against M</b>: ${same ? 'identical, base for base' : fmt.pct(a.identity, 3) + ' identity'} over ${fmt.bp(a.blen)} (${fmt.int(a.matches)} matches)<br>${a.q}:${fmt.int(a.qs)}-${fmt.int(a.qe)} ↔ hap2:${fmt.int(a.ts)}-${fmt.int(a.te)}, ${a.strand} strand, MAPQ ${a.mapq}`);
      el.appendChild(ln);
      el.appendChild(text((X(a.ts) + X(a.te)) / 2 + 10, (Y(a.qs) + Y(a.qe)) / 2 - 8, `${a.piece}${same ? ' = M' : ' ' + fmt.pct(a.identity, 2)}`, { 'font-weight': 700, 'font-size': 10.5, fill: same ? col.mat : col.pat }));
    }
    return el;
  };

  // ---------------------------------------------------------------------------------- 3. depth by parent on the assembly
  P.assemblyDepth = function (container, B, opts = {}) {
    const C = B.case, rows = AR.rows(B.depth_assembly), len = B.local_len, seqs = ['hap1', 'hap2', 'U'].filter((s) => len[s]);
    const x0 = 150, x1 = 1180, k = (x1 - x0) / Math.max(...seqs.map((s) => len[s]));
    const rowH = 120, H = seqs.length * rowH + 20;
    const el = frame(container, H, 'asmdepth');
    const labels = { hap1: ["HPRC's paternal", 'chromosome (hap1)'], hap2: ["HPRC's maternal", 'chromosome (hap2)'], U: ['unplaced contig U', '(hap2)'] };
    const M = B.pieces.find((p) => p.piece === C.assembly.maternal_chromosome);
    seqs.forEach((seq, si) => {
      const top = 10 + si * rowH, base = top + 78, Y = AR.scale(0, 2.4, base, top + 8);
      const X = (p) => x0 + (p - 1) * k;
      el.appendChild(text(x0 - 10, top + 30, labels[seq][0], { 'text-anchor': 'end', class: 'lab' }));
      el.appendChild(text(x0 - 10, top + 44, labels[seq][1], { 'text-anchor': 'end', class: 'mut' }));
      el.appendChild(text(x0 - 10, top + 58, C.assembly.sequences[seq].contig, { 'text-anchor': 'end', class: 'mut', 'font-size': 10 }));
      for (const v of [1, 2]) { el.appendChild(svg('line', { x1: x0, x2: X(len[seq]), y1: Y(v), y2: Y(v), stroke: 'var(--line)', 'stroke-dasharray': '2 3' })); el.appendChild(text(X(len[seq]) + 6, Y(v) + 3, String(v), { class: 'mut', 'font-size': 10 })); }
      // where M holds the missing copy's reads, and C's collapse
      if (seq === 'hap2' && M) {
        for (const [a, b] of B.holes) { if (b - a < 10000) continue; const la = M.local_start + (a - M.rs) * ((M.local_end - M.local_start) / (M.re - M.rs)), lb = M.local_start + (b - M.rs) * ((M.local_end - M.local_start) / (M.re - M.rs));
          el.appendChild(AR.hover(svg('rect', { x: X(la), y: top + 4, width: Math.max(2, X(lb) - X(la)), height: base - top - 4, fill: col.jx, 'fill-opacity': .12 }), `<b>No piece holds the second copy here</b> (GRCh38 ${fmt.range(C.chrom, a, b)}): its reads must land on whatever matches best, M or C.`)); }
      }
      for (const q of B.qc_rows.filter((r) => r[0] === 'nucflag' && C.assembly.sequences[seq].contig === r[1])) {
        const off = C.assembly.sequences[seq].local_offset, a = q[2] + 1 - off, b = q[3] - off; if (b < 1 || a > len[seq]) continue;
        el.appendChild(AR.hover(svg('rect', { x: X(Math.max(1, a)), y: base + 14, width: Math.max(2, X(Math.min(len[seq], b)) - X(Math.max(1, a))), height: 4, fill: q[4] === 'COLLAPSE' ? col.bad : '#f59e0b' }), `<b>NucFlag ${q[4]}</b> ${q[1]}:${fmt.int(q[2] + 1)}-${fmt.int(q[3])}`));
      }
      // bars
      for (const r of rows) {
        if (r.seq !== seq) continue;
        const xa = X(r.bin0 + 1), w = Math.max(0.8, X(r.bin1) - xa);
        const mat = Math.min(2.4, r.mat), pat = Math.min(2.4, r.pat), un = Math.min(2.4, r.unphased), tot = Math.min(2.4, r.mat + r.pat + r.unphased);
        const g = svg('g');
        let y = base;
        for (const [v, c] of [[mat, col.mat], [pat, col.pat], [un, col.grey]]) { if (v <= 0) continue; const hh = base - Y(v); g.appendChild(svg('rect', { x: xa, y: y - hh, width: w, height: hh, fill: c })); y -= hh; }
        if (r.q0 > 0.02) g.appendChild(svg('rect', { x: xa, y: Y(Math.min(2.4, r.q0)) - 1, width: w, height: 2, fill: 'var(--ink)' }));
        AR.hover(g, `<b>${seq}:${fmt.int(r.bin0 + 1)}-${fmt.int(r.bin1)}</b> · piece ${r.piece || '–'}<br>copies of one haplotype's depth: mother's ${fmt.f(r.mat)}, father's ${fmt.f(r.pat)}, unphased ${fmt.f(r.unphased)} (${fmt.f(r.mat + r.pat + r.unphased)} in all)<br>at MAPQ 0: ${fmt.f(r.q0)}`);
        g.classList.add('hit'); g.addEventListener('click', () => opts.onLocus && opts.onLocus(seq, r.bin0 + 1, r.bin1));
        el.appendChild(g);
      }
      // pieces
      for (const [letter, [sname, a, b]] of Object.entries(B.piece_local)) {
        if (sname !== seq) continue; const p = B.pieces.find((q) => q.piece === letter);
        const r = svg('rect', { x: X(a), y: base + 2, width: Math.max(2, X(b) - X(a)), height: 10, rx: 2, fill: parentColor(p.parent) });
        AR.hover(r, `<b>Piece ${letter}</b> ${seq}:${fmt.int(a)}-${fmt.int(b)} (${p.parent} by alleles)`); r.classList.add('hit');
        r.addEventListener('click', () => opts.onLocus && opts.onLocus(seq, a, b));
        el.appendChild(r);
        if (X(b) - X(a) > 12) el.appendChild(text((X(a) + X(b)) / 2, base + 10, letter, { 'text-anchor': 'middle', fill: '#fff', 'font-weight': 700, 'font-size': 9 }));
      }
      // gaps
      B.gaps_local.forEach((gp, i) => {
        if (gp[0] !== seq) return; const j = B.joins[i];
        el.appendChild(svg('line', { x1: X(gp[1]), x2: X(gp[1]), y1: top + 2, y2: base + 12, stroke: 'var(--ink)', 'stroke-dasharray': '3 2' }));
        el.appendChild(AR.hover(text(X(gp[1]), top + 2, `gap ${i + 1}: ${j ? j.realign.cross : 0} cross`, { 'text-anchor': 'middle', 'font-weight': 600, 'font-size': 10 }), j ? `<b>Gap ${i + 1}</b> (${j.left} to ${j.right}), ${seq}:${fmt.int(gp[1])}-${fmt.int(gp[2])}<br>reads realigned across it: ${j.realign.cross}<br>depth within 2 kb of its edges: ${fmt.f(j.realign.edge_left)} / ${fmt.f(j.realign.edge_right)} of one haplotype's` : ''));
      });
      // M's segment ends and the junction reads
      if (seq === 'hap2' && B.m_segment) {
        const [ms, me] = B.m_segment;
        for (const p of [ms, me]) el.appendChild(svg('line', { x1: X(p), x2: X(p), y1: top + 20, y2: base, stroke: col.jx, 'stroke-width': 1.5 }));
        const d = `M${X(ms)},${top + 20} C${X(ms)},${top - 6} ${X(me)},${top - 6} ${X(me)},${top + 20}`;
        const jr = AR.rows(B.reads).filter((r) => r.junction);
        el.appendChild(AR.hover(svg('path', { d, fill: 'none', stroke: col.jx, 'stroke-width': 1.2, 'stroke-dasharray': '5 3' }), `<b>${jr.length} junction reads, realigned to M</b>: ${AR.count(jr, (r) => r.asm_split_end && r.asm_split_start)} split at both of M's segment ends, ${AR.count(jr, (r) => r.asm_through_M)} run through.`));
        el.appendChild(text((X(ms) + X(me)) / 2, top + 3, `${jr.length} junction reads: split, end to start`, { 'text-anchor': 'middle', 'font-weight': 600, 'font-size': 10, fill: col.jx }));
      }
    });
    return el;
  };

  // ---------------------------------------------------------------------------------- 4. the stretch
  P.stretch = function (container, B) {
    const C = B.case, chrom = C.chrom, sites = AR.rows(B.sites), reads = AR.rows(B.reads);
    const readByName = new Map(reads.map((r) => [r.name, r]));
    const st = B.stretch, bounds = B.bounds;
    const blk = B.blocks.find((b) => st && b.start <= st[0] && b.end >= st[1]);
    if (!blk) { container.appendChild(h('p', { class: 'note' }, 'No phase block spans the stretch.')); return; }
    const x0 = 170, x1 = 1180, pad = 6000;
    const X = AR.scale(blk.start - pad, blk.end + pad, x0, x1);
    const topH = 110;
    const el = frame(container, topH, 'stretchtop');
    el.appendChild(text(x0, 14, `Each of the mother's phased sites in her block ${fmt.range(chrom, blk.start, blk.end)} (${blk.n} sites, from her nanopore reads): which of her haplotypes the child's two maternal copies carry`, { class: 'lab' }));
    const rowsY = { 2: 36, 1: 58, 0: 80 }, rowLab = { 2: 'both her haplotype 1', 1: 'one of each', 0: 'both her haplotype 2' }, rowCol = { 2: col.h1, 1: col.jx, 0: col.h2 };
    for (const k of [2, 1, 0]) { el.appendChild(text(x0 - 10, rowsY[k] + 4, rowLab[k], { 'text-anchor': 'end' })); el.appendChild(svg('line', { x1: x0, x2: x1, y1: rowsY[k], y2: rowsY[k], stroke: 'var(--line)', 'stroke-dasharray': '2 3' })); }
    if (st) el.appendChild(svg('rect', { x: X(st[0]), y: 24, width: X(st[1]) - X(st[0]), height: 66, fill: col.jx, 'fill-opacity': .07 }));
    for (const b of bounds) if (b) el.appendChild(svg('line', { x1: X(b), x2: X(b), y1: 24, y2: 90, stroke: 'var(--muted)', 'stroke-dasharray': '3 3' }));
    for (const s of sites) {
      if (s.hap1_copies == null) continue;
      const r = svg('rect', { x: X(s.pos) - 1.2, y: rowsY[s.hap1_copies] - 7, width: 2.4, height: 14, fill: rowCol[s.hap1_copies] });
      AR.hover(r, `<b>${fmt.pos(chrom, s.pos)}</b> ${s.ref}/${s.alt}, her own allele ${s.mother_allele} on her haplotype ${s.own_on_hap}<br>the child's HiFi reads: ${fmt.f(s.share_hifi)} of ${s.hifi_dp} carry it → ${s.copies} of its 3 copies<br>${rowLab[s.hap1_copies]}`);
      el.appendChild(r);
    }
    for (const t of AR.ticks(blk.start - pad, blk.end + pad, 6)) el.appendChild(text(X(t), topH - 6, (t / 1e6).toFixed(2), { 'text-anchor': 'middle', class: 'mut' }));
    // the reads at each end
    const ends = [['left', bounds[0], bounds[1], "the stretch's left end"], ['right', bounds[2], bounds[3], "the stretch's right end"]];
    const grid = h('div', { class: 'cols' });
    container.appendChild(grid);
    const keyLabel = (k) => ({ '2>2': 'stay on her haplotype 2', '2>1': 'go from her haplotype 2 to 1', '1>1': 'stay on her haplotype 1', '1>2': 'go from her haplotype 1 to 2', dad: "the father's (his own alleles)" })[k] || k;
    for (const [side, lo, hi, title] of ends) {
      if (!lo || !hi) continue;
      const list = (B.switch_reads[side] || []).slice();
      const order = ['2>2', '2>1', '1>1', '1>2', 'dad'];
      list.sort((a, b) => order.indexOf(a[1]) - order.indexOf(b[1]) || (readByName.get(a[0]) || {}).start - (readByName.get(b[0]) || {}).start);
      const groups = {}; for (const [n, k] of list) (groups[k] = groups[k] || []).push(n);
      const nReads = list.length, rowPx = 7, groupGap = 22;
      const Hh = 46 + nReads * rowPx + Object.keys(groups).length * groupGap + 24;
      const box = h('div', { class: 'panel tight' });
      const sv = frame(box, Hh, 'stretchend', 600);
      const xs0 = 16, xs1 = 584;
      const lo2 = lo - 24000, hi2 = hi + 24000;
      const XX = AR.scale(lo2, hi2, xs0, xs1);
      sv.appendChild(text(xs0, 16, `${title}: ${nReads} HiFi reads across ${fmt.range(chrom, lo, hi)}`, { class: 'lab' }));
      sv.appendChild(svg('rect', { x: XX(lo), y: 24, width: Math.max(2, XX(hi) - XX(lo)), height: Hh - 48, fill: 'var(--shade)' }));
      let y = 36;
      for (const k of order) {
        const names = groups[k]; if (!names) continue;
        sv.appendChild(text(xs0, y + 6, `${names.length} ${keyLabel(k)}`, { class: 'mut', 'font-weight': 600 })); y += 14;
        for (const n of names) {
          const r = readByName.get(n); if (!r) continue;
          const a = Math.max(lo2, r.start), b = Math.min(hi2, r.end);
          const g = svg('g');
          g.appendChild(svg('line', { x1: XX(a), x2: XX(b), y1: y, y2: y, stroke: k === 'dad' ? col.pat : 'var(--muted)', 'stroke-width': k === 'dad' ? 1.2 : 0.8 }));
          for (const [pos, base] of (B.read_alleles[n] || [])) {
            const s = sites.find((q) => q.pos === pos); if (!s || !s.ref_on_hap || s.block !== blk.id) continue;
            const hap = base === s.ref ? s.ref_on_hap : 3 - s.ref_on_hap;
            g.appendChild(svg('circle', { cx: XX(pos), cy: y, r: 2.2, fill: hap === 1 ? col.h1 : col.h2 }));
          }
          AR.hover(g, `<b>${n}</b> ${fmt.range(chrom, r.start, r.end)} (${fmt.bp(r.length)}) · ${r.parent}${r.phased_by ? ` by ${r.phased_by === 'opposite' ? 'parent-opposite sites' : 'own alleles'}` : ''}<br>alleles: mother's ${r.oh_mat + r.own_mat}, father's ${r.oh_pat + r.own_pat} · ${keyLabel(k)}`);
          sv.appendChild(g); y += rowPx;
        }
        y += groupGap - 14;
      }
      for (const t of AR.ticks(lo2, hi2, 5)) sv.appendChild(text(XX(t), Hh - 6, (t / 1e6).toFixed(3), { 'text-anchor': 'middle', class: 'mut' }));
      grid.appendChild(box);
    }
    container.appendChild(h('p', { class: 'cap' }, 'Dots: green her haplotype 1, violet her haplotype 2, at her phased sites in the block. Blue lines: reads carrying the father\'s own alleles. Grey band: the end\'s interval, between the last site where the copies agree and the first where they differ.'));
    return el;
  };
})();
