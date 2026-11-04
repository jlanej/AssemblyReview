/* AssemblyReview: the read panels, igv.js on three references: GRCh38 (igv's hosted hg38), HPRC's sequence over the
   locus (hprc_local.fa) and the junction as called. igv.js is loaded on demand from jsdelivr, one exact release. Reads are
   coloured by parent through a colour function keyed on the read name, so the BAMs stay exactly as published. */
(function () {
  'use strict';
  const AR = (window.AR = window.AR || {});
  const IGV_URL = 'https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js';
  const Bw = (AR.browsers = { igv: null, b: {}, tracks: {}, pending: {} });

  Bw.canServeRanges = () => location.protocol !== 'file:';

  Bw.load = async function () {
    if (window.igv) return window.igv;
    if (!Bw._loading) Bw._loading = AR.loadScript(IGV_URL).then(() => window.igv);
    return Bw._loading;
  };

  function note(container, html) { AR.clear(container).appendChild(AR.h('div', { class: 'igv-note', html })); }

  function parentMap(B) {
    if (Bw._pm) return Bw._pm;
    const m = new Map();
    for (const r of AR.rows(B.reads)) m.set(r.name, r.junction && r.jr_parent && r.jr_parent !== 'unphased' ? r.jr_parent : r.parent);
    Bw._pm = m; return m;
  }
  Bw.parentColor = function (B) {
    const m = parentMap(B);
    return (aln) => { const p = m.get(aln.readName); return p === 'mother' ? AR.solid.mat : p === 'father' ? AR.solid.pat : AR.solid.grey; };
  };

  Bw.base = (B) => `cases/${B.case.id}/data/`;
  // igv.js's BED decoder does not keep the itemRgb column for these tracks, so each feature is coloured by its name
  const bedColor = (f) => {
    if (f.color) return f.color;
    const n = (f.name || '').toLowerCase();
    if (n.includes('misjoin')) return '#f59e0b';
    if (n.includes('collapse')) return '#c0392b';
    if (n.includes('nnn') || n.startsWith('gap')) return AR.solid.gap;
    if (n.includes('maternal') || n.includes('mother')) return AR.solid.mat;
    if (n.includes('paternal') || n.includes('father')) return AR.solid.pat;
    if (n.includes('junction') || n.includes('missing') || n.includes('segment') || n.includes('parent-opposite')) return AR.solid.jx;
    if (n.includes('scaffold')) return AR.solid.pat;
    return AR.solid.grey;
  };
  const GENES_URL = 'https://genome-browser.s3.us-east-1.amazonaws.com/goldenPath/hg38/database/ncbiRefSeq.txt.gz';
  Bw.addGenes = async function () {
    const browser = Bw.b.grch38; if (!browser) return;
    if (Bw._genes) return;
    Bw._genes = true;
    try { await browser.loadTrack({ type: 'annotation', format: 'refgene', url: GENES_URL, name: 'RefSeq genes (UCSC, from its S3 mirror)', displayMode: 'EXPANDED', height: 80, order: 0, visibilityWindow: -1 }); }
    catch (e) { Bw._genes = false; console.warn('genes track', e); }
  };

  async function create(which, container, config, B) {
    if (!Bw.canServeRanges()) {
      note(container, 'This view needs the page served over HTTP (igv.js reads the BAM and FASTA files by byte range): open the GitHub Pages copy, or run <code>python3 serve.py</code> in the repository and open <code>http://localhost:8000/</code>. The figures need no server.');
      return null;
    }
    note(container, 'Loading igv.js…');
    try {
      const igv = await Bw.load();
      AR.clear(container);
      const div = AR.h('div', { class: 'igv-holder' }); container.appendChild(div);
      const browser = await igv.createBrowser(div, config);
      Bw.b[which] = browser;
      return browser;
    } catch (e) {
      note(container, `The read panel could not load: ${AR.esc(e.message || e)}. igv.js comes from jsdelivr (${AR.esc(IGV_URL)}); the BAM files are read by byte range from this site.`);
      return null;
    }
  }

  Bw.grch38 = async function (container, B, loci) {
    const C = B.case, base = Bw.base(B), [S, E] = C.segment, chrom = C.chrom, color = Bw.parentColor(B);
    const tr = (p, q) => `cases/${C.id}/data/reviewer/tracks/${p}`;
    const reference = C.files.ref_padded
      ? { id: 'GRCh38', name: `GRCh38 ${chrom} (N outside ${AR.fmt.range(chrom, C.windows.WIN[0], C.windows.WIN[1])})`, fastaURL: base + C.files.ref_padded, indexURL: base + C.files.ref_padded + '.fai', compressedIndexURL: base + C.files.ref_padded + '.gzi' }
      : null;
    const config = {
      locus: loci || [AR.fmt.locus(chrom, S - 12000, S + 12000), AR.fmt.locus(chrom, E - 12000, E + 12000)],
      showCenterGuide: false, showCursorTrackingGuide: true,
      tracks: [
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_segment.bed'), name: 'the segment', displayMode: 'COLLAPSED', height: 22, order: 1 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_junction_reads.bed'), name: `junction-spanning HiFi reads (${B.numbers.jr}), both alignments`, displayMode: 'SQUISHED', height: 60, order: 2 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_pieces.bed'), name: 'contig segments on GRCh38 (blue paternal, orange maternal, by alleles)', displayMode: 'EXPANDED', height: 70, order: 3 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_missing.bed'), name: 'second maternal copy: absent from the assembly', displayMode: 'COLLAPSED', height: 22, order: 4 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_nucflag.bed'), name: 'NucFlag calls placed on GRCh38', displayMode: 'COLLAPSED', height: 22, order: 5 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_sites.bed'), name: 'informative SNVs (magenta parents-opposite, orange mother-only, blue father-only, grey both heterozygous)', displayMode: 'COLLAPSED', height: 22, order: 6 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('grch38_calls.bed'), name: 'SV records over the segment', displayMode: 'EXPANDED', height: 60, order: 7 },
        { type: 'alignment', format: 'bam', url: base + C.files.hifi, indexURL: base + C.files.hifi + '.bai', name: `proband ${C.trio.child} HiFi (orange maternal, blue paternal, grey unphased)`, height: 320, colorBy: 'none', color, showSoftClips: true, displayMode: 'EXPANDED', visibilityWindow: 150000, order: 8 },
        { type: 'alignment', format: 'bam', url: base + C.files.sr.child, indexURL: base + C.files.sr.child + '.bai', name: `proband ${C.trio.child} short reads (pair orientation)`, height: 220, colorBy: 'pairOrientation', showSoftClips: true, displayMode: 'SQUISHED', visibilityWindow: 80000, order: 9 },
        { type: 'alignment', format: 'bam', url: base + C.files.sr.father, indexURL: base + C.files.sr.father + '.bai', name: `father ${C.trio.father} short reads`, height: 160, colorBy: 'pairOrientation', showSoftClips: true, displayMode: 'SQUISHED', visibilityWindow: 80000, order: 10 },
        { type: 'alignment', format: 'bam', url: base + C.files.sr.mother, indexURL: base + C.files.sr.mother + '.bai', name: `mother ${C.trio.mother} short reads`, height: 160, colorBy: 'pairOrientation', showSoftClips: true, displayMode: 'SQUISHED', visibilityWindow: 80000, order: 11 },
        { type: 'alignment', format: 'bam', url: base + C.files.ont_mother, indexURL: base + C.files.ont_mother + '.bai', name: `mother ${C.trio.mother} ONT reads`, height: 180, colorBy: 'none', color: AR.solid.grey, showSoftClips: true, displayMode: 'SQUISHED', visibilityWindow: 150000, order: 12 },
      ],
    };
    if (reference) config.reference = reference; else config.genome = 'hg38';
    const browser = await create('grch38', container, config, B);
    if (browser) Bw.highlightJunctionReads('grch38', B, true);
    return browser;
  };

  Bw.assembly = async function (container, B, locus) {
    const C = B.case, base = Bw.base(B), color = Bw.parentColor(B);
    const tr = (p) => `cases/${C.id}/data/reviewer/tracks/${p}`;
    const g1 = B.gaps_local[0];
    const config = {
      reference: { id: 'hprc_local', name: "HPRC's sequence over the locus", fastaURL: base + C.files.asm_local, indexURL: base + C.files.asm_local + '.fai' },
      locus: locus || (g1 ? `${g1[0]}:${g1[1] - 12000}-${g1[2] + 12000}` : 'hap1:580000-600000'), showCenterGuide: false, showCursorTrackingGuide: true,
      tracks: [
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('asm_pieces.bed'), name: 'contig segments (blue paternal, orange maternal, by alleles)', displayMode: 'COLLAPSED', height: 24, order: 1 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('asm_gaps.bed'), name: 'gaps (100 N)', displayMode: 'COLLAPSED', height: 20, order: 2 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('asm_m_segment.bed'), name: 'hap2: the single copy of the segment', displayMode: 'COLLAPSED', height: 20, order: 3 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('asm_nucflag.bed'), name: 'NucFlag (orange misjoin, red collapse)', displayMode: 'COLLAPSED', height: 20, order: 4 },
        { type: 'annotation', format: 'bed', color: bedColor, url: tr('asm_flagger.bed'), name: 'Flagger (HiFi and ONT)', displayMode: 'COLLAPSED', height: 20, order: 5 },
        { type: 'alignment', format: 'bam', url: base + C.files.realign, indexURL: base + C.files.realign + '.bai', name: 'proband HiFi reads realigned to the assembly (MAPQ 0 faint, not hidden)', height: 420, colorBy: 'none', color, showSoftClips: true, displayMode: 'EXPANDED', filter: { mq: 0 }, visibilityWindow: 150000, order: 6 },
      ],
    };
    const browser = await create('assembly', container, config, B);
    if (browser) Bw.highlightJunctionReads('assembly', B, true);
    return browser;
  };

  Bw.junction = async function (container, B) {
    const C = B.case, base = Bw.base(B), color = Bw.parentColor(B);
    if (!C.files.junction_ref) { note(container, 'This case has no junction-as-called reference.'); return null; }
    const config = {
      reference: { id: 'junction', name: 'the junction as called', fastaURL: base + C.files.junction_ref, indexURL: base + C.files.junction_ref + '.fai' },
      locus: 'junction_chr2_48027042-48057041_46970725-47000724:22000-38000', showCenterGuide: true, showCursorTrackingGuide: true,
      tracks: [
        { type: 'alignment', format: 'bam', url: base + C.files.hifi_on_junction, indexURL: base + C.files.hifi_on_junction + '.bai', name: 'proband HiFi reads on the junction as called (join point at 30,000)', height: 360, colorBy: 'none', color, showSoftClips: true, displayMode: 'EXPANDED', visibilityWindow: 150000, order: 1 },
      ],
    };
    const browser = await create('junction', container, config, B);
    if (browser) Bw.highlightJunctionReads('junction', B, true);
    return browser;
  };

  Bw.highlightJunctionReads = function (which, B, on) {
    const browser = Bw.b[which]; if (!browser) return;
    const names = on ? AR.rows(B.reads).filter((r) => r.junction).map((r) => r.name) : [];
    for (const t of browser.findTracks('type', 'alignment')) { try { t.setHighlightedReads(names, AR.solid.jx); } catch (e) { /* older igv */ } }
    browser.repaintViews && browser.repaintViews();
  };
  Bw.highlightReads = function (which, names, color) {
    const browser = Bw.b[which]; if (!browser) return;
    for (const t of browser.findTracks('type', 'alignment')) { try { t.setHighlightedReads(names, color || AR.solid.jx); } catch (e) { /* ignore */ } }
    browser.repaintViews && browser.repaintViews();
  };

  // go to a locus (string or array); creates the browser first if needed through the registered creator
  Bw.goto = async function (which, loci) {
    let browser = Bw.b[which];
    if (!browser && Bw.pending[which]) browser = await Bw.pending[which]();
    if (!browser) return;
    try { await browser.search(Array.isArray(loci) ? loci.join(' ') : loci); } catch (e) { console.warn('igv search failed', e); }
  };
})();
