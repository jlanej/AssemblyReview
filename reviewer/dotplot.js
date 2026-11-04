/* AssemblyReview: a k-mer dot plot drawn on a canvas, no aligner. x: the read (reference orientation); y: one or more
   windows of reference or assembly sequence, stacked. Every run of exact k-mer matches of minRun bp or longer is drawn,
   forward matches in one colour and reverse-strand ones in another. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});

  function kmerIndex(seq, k, maxOcc) {
    const map = new Map();
    for (let i = 0; i + k <= seq.length; i++) {
      const w = seq.substr(i, k);
      if (w.indexOf('N') >= 0) continue;
      let a = map.get(w); if (!a) { a = []; map.set(w, a); } a.push(i);
    }
    for (const [w, a] of map) if (a.length > maxOcc) map.delete(w);
    return map;
  }

  // runs of consecutive matches along each diagonal: [{x0, y0, len, strand}]
  function runs(readIdx, win, k, minRun, strand) {
    const out = [];
    const open = new Map(); // diagonal -> {x0, y0, lastX, len}
    const L = win.length;
    for (let j = 0; j + k <= L; j++) {
      const w = win.substr(j, k);
      const hits = readIdx.get(w); if (!hits) continue;
      for (const i of hits) {
        const d = j - i;
        const o = open.get(d);
        if (o && o.lastJ === j - 1) { o.lastJ = j; o.len++; }
        else { if (o && o.len + k - 1 >= minRun) out.push(fin(o, k, strand, L)); open.set(d, { x0: i, y0: j, lastJ: j, len: 1 }); }
      }
    }
    for (const o of open.values()) if (o.len + k - 1 >= minRun) out.push(fin(o, k, strand, L));
    return out;
  }
  function fin(o, k, strand, L) {
    const len = o.len + k - 1;
    if (strand === '+') return { x0: o.x0, x1: o.x0 + len, y0: o.y0, y1: o.y0 + len, strand };
    // the window was reverse-complemented: position j in it is L-1-j in the window itself
    return { x0: o.x0, x1: o.x0 + len, y0: L - o.y0, y1: L - o.y0 - len, strand };
  }

  AR.dotplot = {
    compute(read, windows, k = 16, minRun = 40) {
      const idx = kmerIndex(read.toUpperCase(), k, 8);
      return windows.map((w) => {
        const seq = w.seq.toUpperCase();
        const fwd = runs(idx, seq, k, minRun, '+');
        const rev = runs(idx, AR.revcomp(seq), k, minRun, '-');
        return { name: w.name, coord: w.coord, len: seq.length, runs: fwd.concat(rev), marks: w.marks || [] };
      });
    },
    // draw on a canvas; windows stacked bottom to top; marks: [{at, label, axis:'y'}] in window coordinates, read marks via opts.xmarks
    draw(canvas, read, computed, opts = {}) {
      const dpr = window.devicePixelRatio || 1;
      const W = opts.width || 560, H = opts.height || 440, mL = 64, mB = 34, mT = 18, mR = 14;
      canvas.width = W * dpr; canvas.height = H * dpr; canvas.style.width = `${W}px`; canvas.style.maxWidth = '100%';
      const ctx = canvas.getContext('2d'); ctx.scale(dpr, dpr);
      const css = getComputedStyle(document.documentElement);
      const ink = css.getPropertyValue('--ink').trim() || '#222', muted = css.getPropertyValue('--muted').trim() || '#666',
        line = css.getPropertyValue('--rule').trim() || '#d9d9d9', panel = '#ffffff';
      ctx.fillStyle = panel; ctx.fillRect(0, 0, W, H);
      const total = computed.reduce((s, w) => s + w.len, 0);
      const gapPx = 10, plotH = H - mT - mB - gapPx * (computed.length - 1), plotW = W - mL - mR;
      const sx = plotW / read.length;
      ctx.font = '11px -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif';
      // x axis
      ctx.strokeStyle = line; ctx.fillStyle = muted; ctx.textAlign = 'center';
      for (const t of AR.ticks(0, read.length, 6)) { const x = mL + t * sx; ctx.beginPath(); ctx.moveTo(x, H - mB); ctx.lineTo(x, H - mB + 4); ctx.stroke(); ctx.fillText(AR.fmt.int(t), x, H - mB + 15); }
      ctx.fillText(opts.xlabel || 'the read (bp, reference orientation)', mL + plotW / 2, H - 6);
      let yTop = mT;
      for (let wi = computed.length - 1; wi >= 0; wi--) {
        const w = computed[wi];
        const h = plotH * (w.len / total);
        const sy = h / w.len;
        const yOf = (j) => yTop + h - j * sy;
        ctx.strokeStyle = line; ctx.strokeRect(mL, yTop, plotW, h);
        // y ticks
        ctx.fillStyle = muted; ctx.textAlign = 'right';
        for (const t of AR.ticks(0, w.len, 3)) { const y = yOf(t); if (y < yTop + 6 || y > yTop + h - 2) continue; ctx.fillText(AR.fmt.int(t + (w.coord ? w.coord[1] : 0)), mL - 4, y + 4); }
        // label
        ctx.save(); ctx.fillStyle = ink; ctx.textAlign = 'left'; ctx.font = 'bold 11px -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif';
        ctx.fillText(w.name, mL + 4, yTop + 12); ctx.restore();
        // marks in the window
        for (const m of w.marks || []) {
          const y = yOf(m.at); ctx.save(); ctx.strokeStyle = opts.markColor || '#d1336e'; ctx.setLineDash([4, 3]);
          ctx.beginPath(); ctx.moveTo(mL, y); ctx.lineTo(mL + plotW, y); ctx.stroke(); ctx.restore();
          if (m.label) { ctx.fillStyle = opts.markColor || '#d1336e'; ctx.textAlign = 'right'; ctx.fillText(m.label, mL + plotW - 4, y - 3); }
        }
        // runs
        for (const r of w.runs) {
          ctx.strokeStyle = r.strand === '+' ? (opts.fwdColor || '#2f6fdb') : (opts.revColor || '#d9a84e'); ctx.lineWidth = 2;
          ctx.beginPath(); ctx.moveTo(mL + r.x0 * sx, yOf(r.y0)); ctx.lineTo(mL + r.x1 * sx, yOf(r.y1)); ctx.stroke();
        }
        ctx.lineWidth = 1;
        yTop += h + gapPx;
      }
      // read marks
      for (const m of opts.xmarks || []) {
        const x = mL + m.at * sx; ctx.save(); ctx.strokeStyle = opts.markColor || '#d1336e'; ctx.setLineDash([4, 3]);
        ctx.beginPath(); ctx.moveTo(x, mT); ctx.lineTo(x, H - mB); ctx.stroke(); ctx.restore();
        if (m.label) { ctx.fillStyle = opts.markColor || '#d1336e'; ctx.textAlign = 'left'; ctx.fillText(m.label, x + 4, H - mB - 6); }
      }
    },
  };
})();
