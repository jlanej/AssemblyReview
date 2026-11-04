/* AssemblyReview: small helpers shared by the page's modules. Classic scripts, one namespace: window.AR. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});
  const SVGNS = 'http://www.w3.org/2000/svg';

  AR.solid = { pat: '#2b6cb0', mat: '#d98c1f', jx: '#c2185b', h1: '#138d75', h2: '#7e57c2', grey: '#9e9e9e', gap: '#111111' };

  AR.fmt = {
    int: (n) => (n == null || isNaN(n) ? '–' : Math.round(n).toLocaleString('en-US')),
    f: (n, d = 2) => (n == null || isNaN(n) ? '–' : Number(n).toFixed(d)),
    bp: (n) => (n == null ? '–' : `${Math.round(n).toLocaleString('en-US')} bp`),
    kb: (n) => (n == null ? '–' : `${(n / 1000).toFixed(n >= 100000 ? 0 : 1)} kb`),
    mb: (n) => (n == null ? '–' : `${(n / 1e6).toFixed(2)} Mb`),
    pos: (chrom, p) => `${chrom}:${Math.round(p).toLocaleString('en-US')}`,
    range: (chrom, a, b) => `${chrom}:${Math.round(a).toLocaleString('en-US')}-${Math.round(b).toLocaleString('en-US')}`,
    locus: (chrom, a, b) => `${chrom}:${Math.round(a)}-${Math.round(b)}`,
    pct: (x, d = 0) => (x == null ? '–' : `${(100 * x).toFixed(d)}%`),
  };

  // elements
  AR.h = function (tag, attrs, ...children) {
    const el = document.createElement(tag);
    if (attrs) for (const [k, v] of Object.entries(attrs)) {
      if (v == null || v === false) continue;
      if (k === 'class') el.className = v;
      else if (k === 'html') el.innerHTML = v;
      else if (k === 'text') el.textContent = v;
      else if (k.startsWith('on') && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else el.setAttribute(k, v === true ? '' : v);
    }
    for (const c of children.flat()) {
      if (c == null || c === false) continue;
      el.appendChild(typeof c === 'string' || typeof c === 'number' ? document.createTextNode(String(c)) : c);
    }
    return el;
  };
  AR.svg = function (tag, attrs, ...children) {
    const el = document.createElementNS(SVGNS, tag);
    if (attrs) for (const [k, v] of Object.entries(attrs)) {
      if (v == null || v === false) continue;
      if (k === 'text') el.textContent = v;
      else if (k.startsWith('on') && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else el.setAttribute(k, v);
    }
    for (const c of children.flat()) if (c != null && c !== false) el.appendChild(c);
    return el;
  };
  AR.clear = (el) => { while (el.firstChild) el.removeChild(el.firstChild); return el; };

  // a bundle table {cols, rows} -> array of objects
  AR.rows = function (t) {
    if (!t) return [];
    if (Array.isArray(t)) return t;
    const cols = t.cols;
    return t.rows.map((r) => { const o = {}; for (let i = 0; i < cols.length; i++) o[cols[i]] = r[i]; return o; });
  };

  // linear scale
  AR.scale = function (d0, d1, r0, r1) {
    const k = (r1 - r0) / (d1 - d0);
    const f = (x) => r0 + (x - d0) * k;
    f.invert = (y) => d0 + (y - r0) / k;
    f.domain = [d0, d1]; f.range = [r0, r1];
    return f;
  };
  AR.ticks = function (d0, d1, n = 6) {
    const span = d1 - d0; if (span <= 0) return [];
    const raw = span / n; const p = Math.pow(10, Math.floor(Math.log10(raw)));
    const steps = [1, 2, 2.5, 5, 10]; let step = p;
    for (const s of steps) { if (raw / (s * p) <= 1) { step = s * p; break; } }
    const out = []; for (let v = Math.ceil(d0 / step) * step; v <= d1; v += step) out.push(+v.toFixed(6));
    return out;
  };

  // tooltip
  const tipEl = () => document.getElementById('tip');
  AR.tip = {
    show(html, x, y) {
      const t = tipEl(); if (!t) return;
      t.innerHTML = html; t.hidden = false;
      const w = t.offsetWidth, hgt = t.offsetHeight;
      let left = x + 14, top = y + 14;
      if (left + w > window.innerWidth - 8) left = x - w - 14;
      if (top + hgt > window.innerHeight - 8) top = y - hgt - 14;
      t.style.left = `${Math.max(4, left)}px`; t.style.top = `${Math.max(4, top)}px`;
    },
    hide() { const t = tipEl(); if (t) t.hidden = true; },
  };
  AR.hover = function (el, html) {
    el.addEventListener('mousemove', (e) => AR.tip.show(typeof html === 'function' ? html() : html, e.clientX, e.clientY));
    el.addEventListener('mouseleave', () => AR.tip.hide());
    return el;
  };

  // scripts, downloads, misc
  AR.loadScript = function (url) {
    return new Promise((resolve, reject) => {
      const s = document.createElement('script'); s.src = url; s.async = true;
      s.onload = () => resolve(); s.onerror = () => reject(new Error(`could not load ${url}`));
      document.head.appendChild(s);
    });
  };
  AR.download = function (name, text, type = 'text/plain') {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  AR.copy = async function (text) { try { await navigator.clipboard.writeText(text); return true; } catch (e) { return false; } };
  AR.median = (xs) => { const a = xs.filter((x) => x != null).sort((p, q) => p - q); const n = a.length; return n ? (n % 2 ? a[(n - 1) / 2] : (a[n / 2 - 1] + a[n / 2]) / 2) : null; };
  AR.mean = (xs) => { const a = xs.filter((x) => x != null); return a.length ? a.reduce((p, q) => p + q, 0) / a.length : null; };
  AR.count = (xs, f) => xs.reduce((n, x) => n + (f(x) ? 1 : 0), 0);
  AR.esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  AR.revcomp = (s) => { const m = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' }; let o = ''; for (let i = s.length - 1; i >= 0; i--) o += m[s[i]] || 'N'; return o; };

  // a sortable table from rows and column specs [{key, label, num, fmt, cls}]
  AR.table = function (rows, cols, opts = {}) {
    const wrap = AR.h('div', { class: 'tbl' + (opts.scroll ? ' scroll' : '') });
    const table = AR.h('table', { class: opts.cls || null });
    const thead = AR.h('thead'); const tbody = AR.h('tbody');
    let sortKey = opts.sort || null, sortDir = opts.dir || 1;
    const render = () => {
      AR.clear(tbody);
      let rs = rows.slice();
      if (sortKey) {
        const c = cols.find((x) => x.key === sortKey);
        rs.sort((a, b) => { const x = a[sortKey], y = b[sortKey]; if (x == null) return 1; if (y == null) return -1;
          return (c && c.num ? x - y : String(x).localeCompare(String(y))) * sortDir; });
      }
      for (const r of rs) {
        const tr = AR.h('tr', { class: opts.onRow ? 'click' : null });
        for (const c of cols) {
          const v = r[c.key];
          const td = AR.h('td', { class: [c.num ? 'num' : '', c.cls || ''].join(' ').trim() || null });
          const out = c.fmt ? c.fmt(v, r) : (v == null ? '' : v);
          if (out instanceof Node) td.appendChild(out); else td.innerHTML = out;
          tr.appendChild(td);
        }
        if (opts.onRow) tr.addEventListener('click', () => { for (const x of tbody.children) x.classList.remove('sel'); tr.classList.add('sel'); opts.onRow(r, tr); });
        if (opts.rowAttrs) { const a = opts.rowAttrs(r); if (a) for (const [k, v] of Object.entries(a)) tr.setAttribute(k, v); }
        tbody.appendChild(tr);
      }
      if (opts.foot) for (const tr of opts.foot) tbody.appendChild(tr);
    };
    const trh = AR.h('tr');
    for (const c of cols) {
      const th = AR.h('th', { class: [c.num ? 'num' : '', c.cls || ''].join(' ').trim() || null, html: c.label + (sortKey === c.key ? (sortDir > 0 ? ' ▲' : ' ▼') : '') });
      th.addEventListener('click', () => { if (sortKey === c.key) sortDir = -sortDir; else { sortKey = c.key; sortDir = c.num ? -1 : 1; }
        for (const t of trh.children) t.innerHTML = t.innerHTML.replace(/ [▲▼]$/, ''); th.innerHTML = c.label + (sortDir > 0 ? ' ▲' : ' ▼'); render(); });
      trh.appendChild(th);
    }
    thead.appendChild(trh); table.appendChild(thead); table.appendChild(tbody); wrap.appendChild(table);
    render();
    return wrap;
  };
})();
