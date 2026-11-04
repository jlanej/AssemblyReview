#!/usr/bin/env python3
"""Build a case's reviewer tables from the files in its data/ folder, and check them against numbers.json.

    python3 -I tools/prepare_case.py cases/hg01943 [--cache DIR] [--align] [--no-check]

Reads  cases/<case>/case.json (coordinates, file paths, thresholds) and cases/<case>/data/ (the slices of the public
releases). Writes cases/<case>/data/reviewer/ : the tables the page loads (TSV), the igv.js tracks (BED), the data
bundle (bundle.js, sequences.js) and checks.tsv, every number of numbers.json against its recount here.

Needs samtools, bcftools and minimap2 on PATH; Python's standard library only. --align also makes the junction as
called (realign/junction_as_called.fa) and the HiFi reads aligned to it (realign/hifi_on_junction.bam).
No path of this machine is written into any output: every path in a header is relative to the case's data folder.
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics
import subprocess
import sys

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
REF_OPS = set('MDN=X')
QRY_OPS = set('MIS=X')


# ----------------------------------------------------------------------------------------------------- helpers

def run(cmd, stderr=None):
    """Run a command and return its stdout as text."""
    return subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=stderr if stderr is not None else subprocess.PIPE,
                          text=True).stdout


def stream(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    try:
        for line in p.stdout:
            yield line.rstrip('\n')
    finally:
        p.stdout.close()
        p.wait()


def cigar_ops(cigar):
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def ref_span(ops):
    return sum(n for n, op in ops if op in REF_OPS)


def query_len(ops):
    return sum(n for n, op in ops if op in QRY_OPS or op == 'H')


def clips(ops):
    left = right = 0
    i = 0
    while i < len(ops) and ops[i][1] in 'SH':
        left += ops[i][0]
        i += 1
    j = len(ops) - 1
    while j >= 0 and ops[j][1] in 'SH':
        right += ops[j][0]
        j -= 1
    return left, right


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name:
        seqs[name] = ''.join(parts)
    return seqs


def fmt(x, nd=4):
    if isinstance(x, float):
        return f'{x:.{nd}f}'
    return str(x)


class Aln:
    __slots__ = ('name', 'flag', 'chrom', 'pos', 'end', 'mapq', 'ops', 'lclip', 'rclip', 'strand', 'qlen', 'seq',
                 'tags', 'primary')

    def __init__(self, fields, keep_seq=False):
        self.name = fields[0]
        self.flag = int(fields[1])
        self.chrom = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.ops = cigar_ops(fields[5])
        self.end = self.pos + ref_span(self.ops) - 1
        self.lclip, self.rclip = clips(self.ops)
        self.strand = '-' if self.flag & 16 else '+'
        self.qlen = query_len(self.ops)
        self.seq = fields[9] if keep_seq else None
        self.tags = fields[11:]
        self.primary = not (self.flag & 0x900)

    def tag(self, key):
        for t in self.tags:
            if t.startswith(key + ':'):
                return t.split(':', 2)[2]
        return None

    def through(self, p, flank):
        return self.pos <= p - flank and self.end >= p + flank


def group_alignments(lines, keep_seq_for=None, by_mate=False):
    reads = collections.OrderedDict()
    for line in lines:
        f = line.split('\t')
        keep = keep_seq_for is not None and (keep_seq_for is True or f[0] in keep_seq_for)
        a = Aln(f, keep_seq=keep)
        key = a.name
        if by_mate and a.flag & 1:
            key = f'{a.name}/{1 if a.flag & 0x40 else 2}'
        reads.setdefault(key, []).append(a)
    return reads


# ------------------------------------------------------------------------------------------------------- sites

def load_sites(cfg, data):
    """The trio's SNVs over WIN, classified: OH parent-opposite, MH mother-only, FH father-only, BH both heterozygous."""
    th = cfg['thresholds']
    samples = run(['bcftools', 'query', '-l', os.path.join(data, cfg['files']['sites_bcf'])]).split()
    trio = cfg['trio']
    col = {role: samples.index(trio[role]) for role in ('child', 'father', 'mother')}
    rows = run(['bcftools', 'query', '-i', f'QUAL>={th["site_qual"]} && TYPE="snp"',
                '-f', '%POS\t%REF\t%ALT\t%QUAL[\t%GT:%AD]\n', os.path.join(data, cfg['files']['sites_bcf'])])
    cand = []
    for line in rows.splitlines():
        f = line.split('\t')
        pos, ref, alt = int(f[0]), f[1], f[2]
        if len(ref) != 1 or len(alt) != 1:
            continue
        g = {}
        ok = True
        for role, i in col.items():
            gt, ad = f[4 + i].split(':')
            ads = [int(x) if x != '.' else 0 for x in ad.split(',')]
            if len(ads) != 2:
                ok = False
                break
            g[role] = (gt, ads[0], ads[1])
        if ok:
            cand.append((pos, ref, alt, float(f[3]), g))
    n_snvs = len(cand)
    med = {}
    for role in ('father', 'mother'):
        med[role] = statistics.median(g[role][1] + g[role][2] for _, _, _, _, g in cand)

    def call(ad_ref, ad_alt):
        dp = ad_ref + ad_alt
        if dp == 0:
            return None
        lo, hi = min(ad_ref, ad_alt), max(ad_ref, ad_alt)
        if lo <= max(1, 0.05 * dp):
            return 'ref' if ad_ref >= ad_alt else 'alt'
        if lo >= dp / 5:
            return 'het'
        return None

    sites = []
    for pos, ref, alt, qual, g in cand:
        calls = {}
        bad = False
        for role in ('father', 'mother'):
            dp = g[role][1] + g[role][2]
            if not (med[role] / 2 <= dp <= med[role] * 2):
                bad = True
                break
            c = call(g[role][1], g[role][2])
            if c is None:
                bad = True
                break
            calls[role] = c
        if bad:
            continue
        fa, mo = calls['father'], calls['mother']
        if fa == 'het' and mo == 'het':
            kind = 'BH'
        elif fa == 'het':
            kind = 'FH'
        elif mo == 'het':
            kind = 'MH'
        elif fa != mo:
            kind = 'OH'
        else:
            continue
        # whose is each allele: at OH each parent's; at MH her own is the one the father lacks; at FH his own
        mother_allele = father_allele = None
        if kind == 'OH':
            mother_allele = ref if mo == 'ref' else alt
            father_allele = ref if fa == 'ref' else alt
        elif kind == 'MH':
            mother_allele = alt if fa == 'ref' else ref
        elif kind == 'FH':
            father_allele = alt if mo == 'ref' else ref
        c = g['child']
        sites.append({'pos': pos, 'ref': ref, 'alt': alt, 'qual': qual, 'kind': kind,
                      'child_gt': c[0], 'child_ad': (c[1], c[2]),
                      'father_gt': g['father'][0], 'father_ad': (g['father'][1], g['father'][2]),
                      'mother_gt': g['mother'][0], 'mother_ad': (g['mother'][1], g['mother'][2]),
                      'mother_allele': mother_allele, 'father_allele': father_allele})
    return sites, n_snvs, med


# ----------------------------------------------------------------------------------------------- pileups

def parse_pileup_bases(s):
    out = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '^':
            i += 2
            continue
        if c == '$':
            i += 1
            continue
        if c in '+-':
            j = i + 1
            while j < n and s[j].isdigit():
                j += 1
            i = j + int(s[i + 1:j])
            continue
        out.append(c.upper())
        i += 1
    return out


def pileup(bam, bed, region, mapq, cache_path=None, min_bq=0, primary_only=False):
    """samtools mpileup at the BED's sites with read names: {pos: [(name, base), ...]}."""
    if cache_path and os.path.exists(cache_path):
        text = open(cache_path).read()
    else:
        ff = ['--ff', 'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY'] if primary_only else []
        text = run(['samtools', 'mpileup', '-l', bed, '-r', region, '-q', str(mapq), '-Q', str(min_bq),
                    '--output-QNAME'] + ff + [bam], stderr=subprocess.DEVNULL)
        if cache_path:
            with open(cache_path, 'w') as f:
                f.write(text)
    out = {}
    for line in text.splitlines():
        f = line.split('\t')
        if len(f) < 7 or int(f[3]) == 0:
            continue
        bases = parse_pileup_bases(f[4])
        names = f[6].split(',')
        if len(bases) != len(names):
            raise RuntimeError(f'pileup parse mismatch at {f[1]}: {len(bases)} bases, {len(names)} names')
        out[int(f[1])] = list(zip(names, bases))
    return out


# ----------------------------------------------------------------------------------------------- phasing

def phase_read(alleles, site_by_pos, th, restrict=None):
    """alleles: {pos: base}. Returns (parent, how, oh_m, oh_f, own_m, own_f)."""
    oh_m = oh_f = own_m = own_f = 0
    for pos, base in alleles.items():
        if restrict and not restrict(pos):
            continue
        s = site_by_pos[pos]
        if s['kind'] == 'OH':
            if base == s['mother_allele']:
                oh_m += 1
            elif base == s['father_allele']:
                oh_f += 1
        elif s['kind'] == 'MH' and base == s['mother_allele']:
            own_m += 1
        elif s['kind'] == 'FH' and base == s['father_allele']:
            own_f += 1
    n = oh_m + oh_f
    if n >= th['phase_opposite_min'] and max(oh_m, oh_f) >= th['phase_opposite_agree'] * n:
        return ('mother' if oh_m > oh_f else 'father'), 'opposite', oh_m, oh_f, own_m, own_f
    if own_m >= th['phase_own_min'] and own_f == 0:
        return 'mother', 'own', oh_m, oh_f, own_m, own_f
    if own_f >= th['phase_own_min'] and own_m == 0:
        return 'father', 'own', oh_m, oh_f, own_m, own_f
    return 'unphased', '', oh_m, oh_f, own_m, own_f


# ----------------------------------------------------------------------------------------------- the assembly

class PieceAln:
    """One record of HPRC's alignment of a haplotype to GRCh38, with contig <-> reference mapping."""

    def __init__(self, a, hap):
        self.hap = hap
        self.contig = a.name.split('#')[-1]
        self.full = a.name
        self.rs, self.re = a.pos, a.end
        if a.strand != '+':
            raise RuntimeError(f'reverse-strand assembly record not handled: {a.name} {a.pos}')
        self.qs = a.lclip + 1
        self.qe = a.lclip + sum(n for n, op in a.ops if op in QRY_OPS)
        self.mapq = a.mapq
        self.kind = 'supplementary' if a.flag & 0x800 else 'primary'
        self.ops = a.ops
        self.nm = int(a.tag('NM') or 0)
        self.aligned = sum(n for n, op in a.ops if op in 'M=X')
        self.seq = a.seq
        # aligned blocks (ref_start, ref_end, contig_start, seq_index) for M/=/X runs
        self.blocks = []
        r, q, si = a.pos, 1, 0
        for n, op in a.ops:
            if op in 'M=X':
                self.blocks.append((r, r + n - 1, q, si))
                r += n
                q += n
                si += n
            elif op in 'DN':
                r += n
            elif op == 'I':
                q += n
                si += n
            elif op == 'S':
                q += n
                si += n
            elif op == 'H':
                q += n
        self.block_starts = [b[0] for b in self.blocks]
        self.cblock_starts = [b[2] for b in self.blocks]

    def ref_to_contig(self, p):
        import bisect
        i = bisect.bisect_right(self.block_starts, p) - 1
        if i < 0:
            return self.blocks[0][2]
        rs, re_, cs, _ = self.blocks[i]
        if p <= re_:
            return cs + (p - rs)
        return cs + (re_ - rs)  # inside a deletion: the last aligned contig base before it

    def contig_to_ref(self, q):
        import bisect
        i = bisect.bisect_right(self.cblock_starts, q) - 1
        if i < 0:
            return self.blocks[0][0]
        rs, re_, cs, _ = self.blocks[i]
        ce = cs + (re_ - rs)
        if q <= ce:
            return rs + (q - cs)
        return re_

    def variants(self, ref_seq, ref_offset, lo, hi):
        """SNVs and indels against GRCh38 over [lo, hi] (1-based GRCh38), from the CIGAR and the sequence."""
        snv, ins, dele = {}, {}, {}
        r, si = self.rs, 0
        seq = self.seq
        for n, op in self.ops:
            if op in 'M=X':
                a, b = max(r, lo), min(r + n - 1, hi)
                if a <= b:
                    for p in range(a, b + 1):
                        rb = ref_seq[p - ref_offset]
                        qb = seq[si + (p - r)]
                        if rb != qb and rb.upper() != qb.upper():
                            snv[p] = (rb.upper(), qb.upper())
                r += n
                si += n
            elif op == 'D':
                if lo <= r <= hi:
                    dele[r] = n
                r += n
            elif op == 'N':
                r += n
            elif op == 'I':
                if lo <= r <= hi:
                    ins[r] = seq[si:si + n]
                si += n
            elif op == 'S':
                si += n
        return snv, ins, dele


def load_bed_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fl = line.rstrip('\n').split('\t')
            rows.append((fl[0].split('#')[-1], int(fl[1]), int(fl[2]), fl[3] if len(fl) > 3 else ''))
    return rows


# -------------------------------------------------------------------------------------------------- calls

def parse_vcf(path):
    header, samples, recs = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('##'):
                header.append(line.rstrip('\n'))
            elif line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]
            else:
                f_ = line.rstrip('\n').split('\t')
                info = {}
                for kv in f_[7].split(';'):
                    if '=' in kv:
                        k, v = kv.split('=', 1)
                        info[k] = v
                    else:
                        info[kv] = True
                gts, fmts = [], []
                if len(f_) > 9:
                    fmt_keys = f_[8].split(':')
                    for col in f_[9:]:
                        vals = col.split(':')
                        gts.append(vals[0] if fmt_keys and fmt_keys[0] == 'GT' else '.')
                        fmts.append(dict(zip(fmt_keys, vals)))
                recs.append({'chrom': f_[0], 'pos': int(f_[1]), 'id': f_[2], 'ref': f_[3], 'alt': f_[4],
                             'qual': f_[5], 'filter': f_[6], 'info': info, 'gts': gts, 'fmt': fmts, 'samples': samples})
    return samples, recs


def sv_of(rec):
    info = rec['info']
    svtype = info.get('SVTYPE')
    alt = rec['alt']
    if svtype is None:
        if alt.startswith('<') and alt.endswith('>'):
            svtype = alt[1:-1].split(':')[0]
        elif '[' in alt or ']' in alt:
            svtype = 'BND'
        else:
            svtype = 'INS' if len(alt) > len(rec['ref']) else ('DEL' if len(alt) < len(rec['ref']) else 'SNV')
    if svtype == 'BND' or '[' in alt or ']' in alt:
        m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', alt)
        mate = (m.group(1), int(m.group(2))) if m else None
        return 'BND', rec['pos'], rec['pos'], 0, mate
    svlen = info.get('SVLEN')
    try:
        svlen = abs(int(str(svlen).split(',')[0])) if svlen is not None else None
    except ValueError:
        svlen = None
    end = info.get('END')
    try:
        end = int(end) if end is not None else None
    except ValueError:
        end = None
    if svlen is None:
        if end is not None and end > rec['pos']:
            svlen = end - rec['pos']
        elif not alt.startswith('<'):
            svlen = max(abs(len(a) - len(rec['ref'])) for a in alt.split(','))
        else:
            svlen = 0
    if end is None:
        end = rec['pos'] + (svlen if svtype in ('DEL', 'DUP', 'INV', 'CNV') else 0)
    return svtype, rec['pos'], end, svlen, None


def carried_len(rec, svlen):
    """For a multiallelic record without SVLEN, the length of the allele the held samples carry (else the longest)."""
    if 'SVLEN' in rec['info'] or rec['alt'].startswith('<') or ',' not in rec['alt']:
        return svlen
    alts = rec['alt'].split(',')
    ks = {int(a) for g in rec['gts'] for a in g.replace('|', '/').split('/') if a.isdigit() and 0 < int(a) <= len(alts)}
    return max(abs(len(alts[k - 1]) - len(rec['ref'])) for k in ks) if ks else svlen


# the fields a record's own caller uses to state its evidence and precision, in the order shown
SUPPORT_INFO = ['SUPPORT', 'SUPPREAD', 'SR', 'PE', 'RE', 'VAF', 'SRQ', 'PRECISE', 'IMPRECISE', 'CIPOS', 'CIEND', 'HOMLEN', 'MAPQ', 'SRMAPQ',
                'CT', 'STRAND', 'STRANDS', 'MATEDIST', 'HAP', 'RNAMES', 'ALGORITHMS', 'EVIDENCE', 'BOTHSIDES_SUPPORT']
SUPPORT_FORMAT = ['GQ', 'DR', 'DV', 'RR', 'RV', 'AD', 'DP', 'RD_CN', 'RDCN', 'PE_GT', 'SR_GT', 'EV']


def support_of(rec):
    """The record's evidence as its caller states it: INFO fields, then each sample's FORMAT fields, as one string."""
    parts = []
    for k in SUPPORT_INFO:
        if k in rec['info']:
            v = rec['info'][k]
            parts.append(k if v is True else f'{k} {v}')
    for name, f in zip(rec['samples'], rec['fmt']):
        vals = [f'{k} {f[k]}' for k in SUPPORT_FORMAT if k in f and f[k].strip('.,')]
        if vals:
            parts.append((f'{name}: ' if len(rec['fmt']) > 1 else '') + ', '.join(vals))
    return '; '.join(parts)


# ================================================================================================= main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case', help='the case folder, e.g. cases/hg01943')
    ap.add_argument('--cache', help='folder for cached pileups (default: <case>/data/reviewer/.cache, not committed)')
    ap.add_argument('--align', action='store_true', help='also make the junction as called and align the HiFi reads to it')
    ap.add_argument('--no-check', action='store_true')
    args = ap.parse_args()

    case = args.case.rstrip('/')
    cfg = json.load(open(os.path.join(case, 'case.json')))
    data = os.path.join(case, 'data')
    out = os.path.join(data, 'reviewer')
    cache = args.cache or os.path.join(out, '.cache')
    os.makedirs(out, exist_ok=True)
    os.makedirs(os.path.join(out, 'tracks'), exist_ok=True)
    os.makedirs(cache, exist_ok=True)
    th = cfg['thresholds']
    chrom = cfg['chrom']
    S, E = cfg['segment']
    WIN, LWIN = cfg['windows']['WIN'], cfg['windows']['LWIN']
    MAPQ = th['mapq']
    P = lambda rel: os.path.join(data, rel)
    numbers = json.load(open(os.path.join(case, 'numbers.json'))) if os.path.exists(os.path.join(case, 'numbers.json')) else {}
    checks = []

    NOTES = {
        'sr_flank': 'raw depth: samtools bedcov -Q 20; the original counted a fraction of a percent differently',
        'sr_rel': 'per-bin ratio of copies, child against the parents\' mean, over the 107 bins fully inside; rests on the raw depth',
        'sr_father': 'rests on the raw depth', 'sr_mother': 'rests on the raw depth', 'hifi_in': 'rests on the raw depth',
        'hifi_mat': 'bins more than 10 kb inside with 10 or more reads of each parent phased at parent-opposite sites; the original\'s rule is unpublished',
        'hifi_pat': 'as hifi_mat', 'hifi_phased_span': 'as hifi_mat', 'hifi_flank_bins': 'as hifi_mat', 'hifi_flank_mat': 'as hifi_mat',
        'hifi_flank_pat': 'as hifi_mat', 'consistency': 'alleles read from the primary alignments',
        'ont.len': 'median primary alignment length at MAPQ 20', 'ont.blocks': 'the block spanning the stretch is identical; the last block differs by 3 sites',
        'ont.linked': 'reads carrying two or more of her heterozygous sites',
        'switch_left': 'rule: two or more phased sites each side, two in three agreeing; the father\'s = his own alleles and none of hers',
        'switch_right': 'as switch_left', 'realign.per_piece': 'rests on each read\'s parent, read from its primary alignment',
        'realign.one_copy': 'the left flank from 20 kb inside the long-read window to 10 kb before the segment',
        'realign.flank_pat': 'as realign.one_copy', 'realign.flank_mat': 'as realign.one_copy',
    }

    def check(key, expected, got, note=''):
        if not note:
            for k, v in NOTES.items():
                if key == k or key.startswith(k + '.') or (k.startswith('realign.joins') and key.startswith(k)):
                    note = v
            if key.startswith('realign.joins.') and ('edge_' in key):
                note = NOTES['realign.one_copy']
        ok = expected == got or (isinstance(expected, float) and isinstance(got, (int, float)) and abs(expected - got) <= 1e-9 * max(1, abs(expected)))
        if isinstance(expected, float) and isinstance(got, (int, float)) and not ok:
            ok = abs(expected - got) < 5e-4
        checks.append({'key': key, 'expected': expected, 'got': got, 'ok': ok, 'note': note})
        return ok

    log = lambda *a: print(*a, file=sys.stderr, flush=True)

    # ---------------------------------------------------------------------------------------- 1. sites
    log('sites')
    sites, n_snvs, med_depth = load_sites(cfg, data)
    site_by_pos = {s['pos']: s for s in sites}
    kinds = collections.Counter(s['kind'] for s in sites)
    kinds_in = collections.Counter(s['kind'] for s in sites if S <= s['pos'] <= E)
    check('snvs', numbers.get('snvs'), n_snvs)
    for k in ('OH', 'MH', 'FH', 'BH'):
        check(f'sites.{k}', numbers.get('sites', {}).get(k), kinds[k])
        check(f'sites_in.{k}', numbers.get('sites_in', {}).get(k), kinds_in[k])
    last_oh = max(s['pos'] for s in sites if s['kind'] == 'OH' and s['pos'] <= E)
    sites_bed = os.path.join(cache, 'sites.bed')
    with open(sites_bed, 'w') as f:
        for s in sites:
            f.write(f'{chrom}\t{s["pos"] - 1}\t{s["pos"]}\n')

    # ---------------------------------------------------------------------------------------- 2. HiFi reads on GRCh38
    log('hifi reads')
    hifi_bam = P(cfg['files']['hifi'])
    want_seq = set([cfg['dotplot']['read']] + cfg['dotplot']['join_reads'])
    hifi = group_alignments(stream(['samtools', 'view', '-F', '0x704', hifi_bam]), keep_seq_for=want_seq)
    # junction reads need their sequence too: find them first, then fetch their sequences
    slack, minclip, flank = th['junction_slack'], th['junction_clip'], th['through_flank']

    def junction_parts(alns, p_end, p_start, slack, minclip, mapq):
        ends = [a for a in alns if a.mapq >= mapq and abs(a.end - p_end) <= slack and a.rclip >= minclip]
        starts = [a for a in alns if a.mapq >= mapq and abs(a.pos - p_start) <= slack and a.lclip >= minclip]
        for a in ends:
            for b in starts:
                if a.strand == b.strand and a is not b:
                    return a, b
        return (ends[0] if ends else None), (starts[0] if starts else None)

    reads = {}
    jr = []
    for name, alns in hifi.items():
        prim = [a for a in alns if a.primary]
        a_end, b_start = junction_parts(alns, E, S, slack, minclip, MAPQ)
        is_jr = a_end is not None and b_start is not None and a_end.strand == b_start.strand and a_end is not b_start
        rec = {'name': name, 'cell': name.split('/')[0],
               'instrument': 'Revio' if name.startswith('m84') else ('Sequel II' if name.startswith('m64') else '?'),
               'strand': prim[0].strand if prim else alns[0].strand,
               'length': max(a.qlen for a in alns), 'mapq': prim[0].mapq if prim else alns[0].mapq,
               'start': min(a.pos for a in alns), 'end': max(a.end for a in alns), 'n_aln': len(alns),
               'primary': bool(prim), 'junction': is_jr,
               'stop': (not is_jr) and ((a_end is not None) != (b_start is not None)),
               'through_S': any(a.mapq >= MAPQ and a.through(S, flank) for a in alns),
               'through_E': any(a.mapq >= MAPQ and a.through(E, flank) for a in alns)}
        if is_jr:
            rec.update({'side1': (a_end.pos, a_end.end, a_end.mapq), 'side2': (b_start.pos, b_start.end, b_start.mapq)})
            jr.append(name)
        reads[name] = rec
    log(f'  {len(reads)} reads, {sum(1 for r in reads.values() if r["primary"])} primary, {len(jr)} junction reads')
    check('jr', numbers.get('jr'), len(jr))
    check('jr_mapq60', numbers.get('jr_mapq60'), sum(1 for n in jr if reads[n]['side1'][2] == 60 and reads[n]['side2'][2] == 60))
    check('jr_starts', numbers.get('jr_starts'), len({reads[n]['side1'][0] for n in jr}))
    check('jr_ends', numbers.get('jr_ends'), len({reads[n]['side2'][1] for n in jr}))
    side1 = [reads[n]['side1'][1] - reads[n]['side1'][0] + 1 for n in jr]
    side2 = [reads[n]['side2'][1] - reads[n]['side2'][0] + 1 for n in jr]
    check('jr_left', numbers.get('jr_left'), [min(side1), max(side1)] if jr else None)
    check('jr_right', numbers.get('jr_right'), [min(side2), max(side2)] if jr else None)
    check('jr_runs', numbers.get('jr_runs'), dict(collections.Counter(reads[n]['cell'] for n in jr)))
    check('through_e', numbers.get('through_e'), sum(1 for r in reads.values() if r['through_E']))
    check('hifi_stop', numbers.get('hifi_stop'), sum(1 for r in reads.values() if r['stop']))
    # sequences of the junction reads and the dot-plot reads (primary record = whole read, reference orientation)
    seq_names = set(jr) | want_seq
    seqs_reads = {}
    for line in stream(['samtools', 'view', '-F', '0x904', hifi_bam]):
        f = line.split('\t', 10)
        if f[0] in seq_names:
            seqs_reads[f[0]] = f[9]

    # ---------------------------------------------------------------------------------------- 3. HiFi pileup, read phase, site counts
    log('hifi pileup at the informative sites')
    region = f'{chrom}:{LWIN[0]}-{LWIN[1]}'
    # every read's parent rests on its primary alignment's alleles; a junction read's on both its sides, within 30 kb
    pile = pileup(hifi_bam, sites_bed, region, MAPQ, os.path.join(cache, 'hifi_primary.pileup'), primary_only=True)
    pile_all = pileup(hifi_bam, sites_bed, region, MAPQ, os.path.join(cache, 'hifi_all.pileup'))
    read_alleles = collections.defaultdict(dict)
    read_alleles_all = collections.defaultdict(dict)
    for src, dst in ((pile, read_alleles), (pile_all, read_alleles_all)):
        for pos, rows in src.items():
            s = site_by_pos[pos]
            for name, base in rows:
                if base in (s['ref'], s['alt']) and pos not in dst[name]:
                    dst[name][pos] = base
    jwin = th['junction_phase_window']
    near_bp = lambda pos: abs(pos - S) <= jwin or abs(pos - E) <= jwin
    for name, rec in reads.items():
        al = read_alleles.get(name, {})
        parent, how, oh_m, oh_f, own_m, own_f = phase_read(al, site_by_pos, th)
        rec.update({'parent': parent, 'phased_by': how, 'oh_mat': oh_m, 'oh_pat': oh_f, 'own_mat': own_m, 'own_pat': own_f,
                    'alleles': al, 'jr_parent': None})
        if rec['junction']:
            al2 = read_alleles_all.get(name, {})
            jp, jhow, j_oh_m, j_oh_f, j_own_m, j_own_f = phase_read(al2, site_by_pos, th, restrict=near_bp)
            rec.update({'jr_parent': jp, 'jr_mat_alleles': j_oh_m + j_own_m, 'jr_pat_alleles': j_oh_f + j_own_f,
                        'alleles_all': {p: b for p, b in al2.items() if near_bp(p)}})
    jr_parent = collections.Counter(reads[n]['jr_parent'] for n in jr)
    check('jr_mat', numbers.get('jr_mat'), jr_parent['mother'])
    check('jr_pat', numbers.get('jr_pat'), jr_parent['father'])
    check('jr_unphased', numbers.get('jr_unphased'), jr_parent['unphased'])
    mom_counts = [reads[n]['jr_mat_alleles'] for n in jr if reads[n]['jr_parent'] == 'mother']
    check('jr_mom_own', numbers.get('jr_mom_own'), [min(mom_counts), max(mom_counts)] if mom_counts else None)
    check('jr_oh_m', numbers.get('jr_oh_m'), sum(mom_counts))
    check('jr_dad_own', numbers.get('jr_dad_own'), sum(reads[n]['jr_pat_alleles'] for n in jr))
    ts = collections.Counter(reads[n]['parent'] for n in reads if reads[n]['through_S'])
    check('through_s.mat', numbers.get('through_s', {}).get('mat'), ts['mother'])
    check('through_s.pat', numbers.get('through_s', {}).get('pat'), ts['father'])
    # the phasing checks itself
    cons = collections.Counter()
    for rec in reads.values():
        if rec['phased_by'] != 'opposite':
            continue
        key = 'mat' if rec['parent'] == 'mother' else 'pat'
        cons[f'{key}|reads'] += 1
        if rec['own_mat']:
            cons[f'{key}|with a mother-only allele'] += 1
        if rec['own_pat']:
            cons[f'{key}|with a father-only allele'] += 1
        cons[f'{key}|mother-only alleles'] += rec['own_mat']
        cons[f'{key}|father-only alleles'] += rec['own_pat']
    for k, v in numbers.get('consistency', {}).items():
        check(f'consistency.{k}', v, cons[k])
    # site counts in HiFi
    for s in sites:
        rows = pile.get(s['pos'], [])
        s['hifi_ref'] = sum(1 for _, b in rows if b == s['ref'])
        s['hifi_alt'] = sum(1 for _, b in rows if b == s['alt'])
        s['hifi_dp'] = len(rows)

    def share(s, which):
        """The child's share of the given allele: (short reads, HiFi)."""
        c_ref, c_alt = s['child_ad']
        if which is None:
            return None, None
        is_alt = which == s['alt']
        sr = (c_alt if is_alt else c_ref) / (c_ref + c_alt) if c_ref + c_alt else None
        hf = (s['hifi_alt'] if is_alt else s['hifi_ref']) / (s['hifi_ref'] + s['hifi_alt']) if s['hifi_ref'] + s['hifi_alt'] else None
        return sr, hf

    for s in sites:
        s['share_sr'], s['share_hifi'] = share(s, s['mother_allele']) if s['kind'] in ('OH', 'MH') else (None, None)
        s['in_segment'] = S <= s['pos'] <= E
    oh_in = [s for s in sites if s['kind'] == 'OH' and s['in_segment']]
    oh_out = [s for s in sites if s['kind'] == 'OH' and not s['in_segment']]
    mean = lambda xs: sum(xs) / len(xs) if xs else None
    median = lambda xs: statistics.median(xs) if xs else None
    minh = th['site_min_hifi']
    shares = {}
    for label, group in (('oh_in', oh_in), ('oh_out', oh_out)):
        sr = [s['share_sr'] for s in group if s['share_sr'] is not None]
        hf = [s['share_hifi'] for s in group if s['share_hifi'] is not None]
        shares[label] = [median(sr), len(sr), median(hf), len(hf)]
        for i, k in enumerate(('sr', 'sr_n', 'hifi', 'hifi_n')):
            check(f'{label}.{k}', numbers.get(label, [None] * 4)[i], shares[label][i])
    check('last_oh', numbers.get('last_oh'), last_oh)

    # ---------------------------------------------------------------------------------------- 4. depth in copies
    log('depth')

    def bedcov(bed, bams):
        rows = []
        for line in run(['samtools', 'bedcov', '-Q', str(MAPQ), bed] + bams).splitlines():
            f = line.split('\t')
            rows.append((int(f[1]), int(f[2]), [int(x) / (int(f[2]) - int(f[1])) for x in f[3:]]))
        return rows

    sr_bins = bedcov(P(cfg['files']['sr_bins']), [P(cfg['files']['sr'][m]) for m in ('child', 'father', 'mother')])
    hifi_bins = bedcov(P(cfg['files']['hifi_bins']), [hifi_bam])
    is_flank = lambda b0, b1: b1 <= S - 10000 or b0 + 1 >= E + 10000
    is_inside = lambda b0, b1: b0 + 1 >= S and b1 <= E
    overlaps = lambda b0, b1: b1 >= S and b0 + 1 <= E
    sr_flank = {}
    for i, m in enumerate(('child', 'father', 'mother')):
        sr_flank[m] = mean([d[i] for b0, b1, d in sr_bins if is_flank(b0, b1)])
        check(f'sr_flank.{m}', numbers.get('sr_flank', {}).get(m), round(sr_flank[m], 4))
    depth_rows = []
    rel, own = [], {'child': [], 'father': [], 'mother': []}
    for b0, b1, d in sr_bins:
        cp = {m: 2 * d[i] / sr_flank[m] for i, m in enumerate(('child', 'father', 'mother'))}
        row = {'bin0': b0, 'bin1': b1, 'sr_child': d[0], 'sr_father': d[1], 'sr_mother': d[2],
               'cp_child': cp['child'], 'cp_father': cp['father'], 'cp_mother': cp['mother'],
               'cp_rel': cp['child'] / ((cp['father'] + cp['mother']) / 2), 'inside': is_inside(b0, b1), 'flank': is_flank(b0, b1)}
        depth_rows.append(row)
        if is_inside(b0, b1):
            rel.append(row['cp_rel'])
            for m in own:
                own[m].append(cp[m] / 2)
    check('sr_rel', numbers.get('sr_rel'), mean(rel))
    check('sr_child', numbers.get('sr_child'), mean(own['child']))
    check('sr_father', numbers.get('sr_father'), mean(own['father']))
    check('sr_mother', numbers.get('sr_mother'), mean(own['mother']))
    hifi_flank = mean([d[0] for b0, b1, d in hifi_bins if is_flank(b0, b1)])
    hrow = {b0: d[0] for b0, b1, d in hifi_bins}
    # the reads phased at parent-opposite sites, by their primary alignments: bases and reads per bin and parent
    phased_bases = collections.defaultdict(collections.Counter)
    phased_reads = collections.defaultdict(collections.Counter)
    for name, rec in reads.items():
        if rec['phased_by'] != 'opposite':
            continue
        for a in hifi[name]:
            if not a.primary or a.mapq < MAPQ:
                continue
            for b in range((a.pos - 1) // 10000, (a.end - 1) // 10000 + 1):
                lo, hi = max(a.pos, b * 10000 + 1), min(a.end, (b + 1) * 10000)
                phased_bases[b * 10000][rec['parent']] += hi - lo + 1
                phased_reads[b * 10000][rec['parent']] += 1
    hifi_in, mat_in, pat_in, mat_fl, pat_fl = [], [], [], [], []
    phased_span = []
    deep_inside = lambda b0, b1: b0 + 1 >= S + 10000 and b1 <= E - 10000
    for row in depth_rows:
        b0, b1 = row['bin0'], row['bin1']
        if b0 in hrow:
            row['hifi'] = hrow[b0]
            row['cp_hifi'] = 2 * hrow[b0] / hifi_flank
            pr, pb = phased_reads.get(b0, {}), phased_bases.get(b0, {})
            row['phased_reads_mat'], row['phased_reads_pat'] = pr.get('mother', 0), pr.get('father', 0)
            if pr.get('mother', 0) >= 10 and pr.get('father', 0) >= 10:
                m, p = pb['mother'], pb['father']
                row['cp_hifi_mat'] = row['cp_hifi'] * m / (m + p)
                row['cp_hifi_pat'] = row['cp_hifi'] * p / (m + p)
                if deep_inside(b0, b1):
                    mat_in.append(row['cp_hifi_mat'])
                    pat_in.append(row['cp_hifi_pat'])
                    phased_span.append(b0)
                elif row['flank']:
                    mat_fl.append(row['cp_hifi_mat'])
                    pat_fl.append(row['cp_hifi_pat'])
            if row['inside']:
                hifi_in.append(row['cp_hifi'])
    check('hifi_in', numbers.get('hifi_in'), mean(hifi_in))
    check('hifi_mat', numbers.get('hifi_mat'), mean(mat_in))
    check('hifi_pat', numbers.get('hifi_pat'), mean(pat_in))
    check('hifi_phased_bins', numbers.get('hifi_phased_bins'), len(mat_in))
    check('hifi_phased_span', numbers.get('hifi_phased_span'), [min(phased_span) + 1, max(phased_span) + 10000] if phased_span else None)
    check('hifi_flank_bins', numbers.get('hifi_flank_bins'), len(mat_fl))
    check('hifi_flank_mat', numbers.get('hifi_flank_mat'), mean(mat_fl))
    check('hifi_flank_pat', numbers.get('hifi_flank_pat'), mean(pat_fl))
    check('hifi_len', numbers.get('hifi_len'), int(statistics.median(r['length'] for r in reads.values() if r['primary'])))

    # ---------------------------------------------------------------------------------------- 5. the trio's short reads at the junction
    log('short reads at the junction')
    sr_evidence = {}
    sr_pairs_rows = []
    for m in ('child', 'father', 'mother'):
        bam = P(cfg['files']['sr'][m])
        regions = [f'{chrom}:{S - 2000}-{S + 2000}', f'{chrom}:{E - 2000}-{E + 2000}']
        near = group_alignments(stream(['samtools', 'view', '-F', '0x704', bam] + regions), by_mate=True)
        split, split_mapq, clip_s, clip_e = [], [], 0, 0
        for name, alns in near.items():
            a_end, b_start = junction_parts(alns, E, S, 5, 10, 0)
            if a_end is not None and b_start is not None and a_end is not b_start and a_end.strand == b_start.strand:
                split.append(name)
                split_mapq.append(a_end.mapq if a_end.primary else b_start.mapq)
            for a in alns:
                if a.primary and a.pos == S and a.lclip >= 10:
                    clip_s += 1
                if a.primary and a.end == E and a.rclip >= 10:
                    clip_e += 1
        # pairs in a tandem duplication's orientation: the forward mate ends within 1 kb before E, its reverse mate starts within 1 kb after S
        by_pair = group_alignments(stream(['samtools', 'view', '-F', '0x704', bam] + regions))
        pairs = []
        for name, alns in by_pair.items():
            fw = [a for a in alns if a.primary and a.strand == '+' and a.flag & 1 and E - 1000 <= a.end <= E + 50]
            rv = [a for a in alns if a.primary and a.strand == '-' and a.flag & 1 and S - 50 <= a.pos <= S + 1000]
            if fw and rv and max(fw[0].mapq, rv[0].mapq) >= MAPQ:
                pairs.append((name, fw[0], rv[0]))
        sr_evidence[m] = {'split': len(split), 'split_mapq': sorted(split_mapq), 'pairs': len(pairs),
                          'pair_mapq_med': statistics.median([max(a.mapq, b.mapq) for _, a, b in pairs]) if pairs else None,
                          'clip_s': clip_s, 'clip_e': clip_e, 'split_names': split}
        for name, a, b in pairs:
            sr_pairs_rows.append({'member': m, 'name': name, 'fwd_start': a.pos, 'fwd_end': a.end, 'fwd_mapq': a.mapq,
                                  'rev_start': b.pos, 'rev_end': b.end, 'rev_mapq': b.mapq})
    check('split', numbers.get('split'), sr_evidence['child']['split'])
    check('split_mapq', numbers.get('split_mapq'), sr_evidence['child']['split_mapq'])
    check('pairs', numbers.get('pairs'), sr_evidence['child']['pairs'])
    check('pair_mapq_med', numbers.get('pair_mapq_med'), sr_evidence['child']['pair_mapq_med'])
    check('clip_s', numbers.get('clip_s'), sr_evidence['child']['clip_s'])
    check('clip_e', numbers.get('clip_e'), sr_evidence['child']['clip_e'])
    for m in ('father', 'mother'):
        check(f'parents_sr.{m}', numbers.get('parents_sr', {}).get(m), [sr_evidence[m]['split'], sr_evidence[m]['pairs'], sr_evidence[m]['clip_s'] + sr_evidence[m]['clip_e']])

    # ---------------------------------------------------------------------------------------- 6. the mother's nanopore reads
    log("the mother's nanopore reads")
    ont_bam = P(cfg['files']['ont_mother'])
    ont = group_alignments(stream(['samtools', 'view', '-F', '0x704', ont_bam]))
    ont_rows = []
    ont_j = ont_ts = ont_te = ont_reach = 0
    for name, alns in ont.items():
        a_end, b_start = junction_parts(alns, E, S, slack, minclip, MAPQ)
        is_j = a_end is not None and b_start is not None and a_end is not b_start and a_end.strand == b_start.strand
        t_s = any(a.mapq >= MAPQ and a.through(S, flank) for a in alns)
        t_e = any(a.mapq >= MAPQ and a.through(E, flank) for a in alns)
        reach = any(a.mapq >= MAPQ and (a.pos - 200 <= S <= a.end + 200 or a.pos - 200 <= E <= a.end + 200) for a in alns)
        ont_j += is_j
        ont_ts += t_s
        ont_te += t_e
        ont_reach += reach
        prim = [a for a in alns if a.primary] or alns
        if max(a.mapq for a in alns) < MAPQ:
            continue
        ont_rows.append({'name': name, 'length': max(a.qlen for a in alns), 'mapq': prim[0].mapq, 'start': min(a.pos for a in alns),
                         'end': max(a.end for a in alns), 'junction': is_j, 'through_S': t_s, 'through_E': t_e, 'reach': reach})
    check('ont.junction', numbers.get('ont', {}).get('junction'), ont_j)
    check('ont.through_s', numbers.get('ont', {}).get('through_s'), ont_ts)
    check('ont.through_e', numbers.get('ont', {}).get('through_e'), ont_te)
    check('ont.reads', numbers.get('ont', {}).get('reads'), len(ont_rows))
    check('ont.reach', numbers.get('ont', {}).get('reach'), ont_reach)
    check('ont.len', numbers.get('ont', {}).get('len'), int(statistics.median(a.qlen for alns in ont.values() for a in alns if a.primary and a.mapq >= MAPQ)))
    # her heterozygous sites, phased from her reads
    het_sites = [s for s in sites if s['kind'] in ('MH', 'BH') and LWIN[0] <= s['pos'] <= LWIN[1]]
    het_bed = os.path.join(cache, 'mother_het.bed')
    with open(het_bed, 'w') as f:
        for s in het_sites:
            f.write(f'{chrom}\t{s["pos"] - 1}\t{s["pos"]}\n')
    opile = pileup(ont_bam, het_bed, region, MAPQ, os.path.join(cache, 'ont.pileup'))
    ont_alleles = collections.defaultdict(dict)
    for pos, rows in opile.items():
        s = site_by_pos[pos]
        for name, base in rows:
            if base in (s['ref'], s['alt']) and pos not in ont_alleles[name]:
                ont_alleles[name][pos] = base
    votes = collections.defaultdict(lambda: [0, 0])  # (i,j) -> [same, different]
    linked_reads = set()
    for name, al in ont_alleles.items():
        ps = sorted(al)
        if len(ps) < 2:
            continue
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                a, b = ps[i], ps[j]
                same = (al[a] == site_by_pos[a]['ref']) == (al[b] == site_by_pos[b]['ref'])
                votes[(a, b)][0 if same else 1] += 1
    links = []
    for (a, b), (same, diff) in votes.items():
        n = same + diff
        if n >= th['ont_link_min_reads'] and max(same, diff) >= th['ont_link_agree'] * n:
            links.append((n, a, b, same >= diff))
    links.sort(key=lambda x: (-x[0], x[1], x[2]))
    parent = {}
    parity = {}  # parity relative to root: 0 = ref alleles on the same haplotype as the root's ref

    def find(x):
        px = parent.get(x, x)
        if px == x:
            return x, 0
        r, p = find(px)
        parent[x] = r
        parity[x] = (parity.get(x, 0) + p) % 2
        return r, parity[x]

    conflicts = 0
    for n, a, b, same in links:
        ra, pa = find(a)
        rb, pb = find(b)
        want = 0 if same else 1
        if ra == rb:
            if (pa + pb) % 2 != want:
                conflicts += 1
            continue
        parent[rb] = ra
        parity[rb] = (pa + pb + want) % 2
    link_pairs = {(a, b) for _, a, b, _ in links}
    linked_reads = {name for name, al in ont_alleles.items() if len(al) >= 2}
    comps = collections.defaultdict(list)
    for s in het_sites:
        if s['pos'] in parent or any(s['pos'] in (a, b) for (a, b) in link_pairs):
            r, p = find(s['pos'])
            comps[r].append((s['pos'], p))
    blocks = []
    for r, members in comps.items():
        members.sort()
        if len(members) < 2:
            continue
        # haplotype 1 carries the reference allele at the block's first site
        p0 = members[0][1]
        blocks.append({'start': members[0][0], 'end': members[-1][0], 'n': len(members),
                       'sites': {pos: (1 if (p - p0) % 2 == 0 else 2) for pos, p in members}})
    blocks.sort(key=lambda b: b['start'])
    big = [b for b in blocks if b['n'] >= 10]
    check('ont.blocks', numbers.get('ont', {}).get('blocks'), [[b['start'], b['end'], b['n']] for b in big],
          f'{len(blocks)} blocks of 2 or more sites; those of 10 or more listed')
    check('ont.conflicts', numbers.get('ont', {}).get('conflicts'), conflicts)
    check('ont.linked', numbers.get('ont', {}).get('linked'), len(linked_reads))
    for bi, b in enumerate(blocks):
        for pos, hap in b['sites'].items():
            site_by_pos[pos]['block'] = bi
            site_by_pos[pos]['ref_on_hap'] = hap  # the haplotype carrying the reference allele
    for s in sites:
        if 'ref_on_hap' in s and s['mother_allele'] is not None:
            s['own_on_hap'] = s['ref_on_hap'] if s['mother_allele'] == s['ref'] else 3 - s['ref_on_hap']

    # ---------------------------------------------------------------------------------------- 7. where the copy was made: the stretch
    log('the stretch')
    for s in sites:
        s['copies'] = None
        if s['kind'] == 'MH' and s['share_hifi'] is not None and s['hifi_dp'] >= minh:
            x = s['share_hifi']
            s['copies'] = 0 if x < 1 / 6 else (1 if x <= 0.42 else (None if x <= 0.58 else (2 if x <= 0.83 else 3)))
            s['agree'] = None if s['copies'] is None else (s['copies'] in (0, 2))
    mh_seg = [s for s in sites if s['kind'] == 'MH' and s['in_segment'] and s.get('copies') is not None]
    # the stretch: the longest run of 'differ' sites, merged across at most two 'agree' sites
    best, cur, gap = [], [], 0
    for s in mh_seg:
        if s['copies'] == 1:
            cur.append(s)
            gap = 0
        else:
            if cur:
                gap += 1
                if gap > 2:
                    if len([x for x in cur if x['copies'] == 1]) > len([x for x in best if x['copies'] == 1]):
                        best = cur
                    cur, gap = [], 0
                else:
                    cur.append(s)
    if cur and len([x for x in cur if x['copies'] == 1]) > len([x for x in best if x['copies'] == 1]):
        best = cur
    while best and best[-1]['copies'] != 1:
        best.pop()
    stretch = [best[0]['pos'], best[-1]['pos']] if best else None
    in_stretch = [s for s in mh_seg if stretch and stretch[0] <= s['pos'] <= stretch[1]]
    check('stretch', numbers.get('stretch'), stretch)
    check('stretch_differ', numbers.get('stretch_differ'), sum(1 for s in in_stretch if s['copies'] == 1))
    check('stretch_n', numbers.get('stretch_n'), len(in_stretch))
    check('agree_before', numbers.get('agree_before'), sum(1 for s in mh_seg if stretch and s['pos'] < stretch[0] and s['copies'] in (0, 2)))
    check('agree_after', numbers.get('agree_after'), sum(1 for s in mh_seg if stretch and s['pos'] > stretch[1] and s['copies'] in (0, 2)))
    differ_outside = [s['pos'] for s in mh_seg if stretch and not (stretch[0] <= s['pos'] <= stretch[1]) and s['copies'] == 1]
    check('differ_outside', numbers.get('differ_outside'), differ_outside)
    before = [s for s in mh_seg if stretch and s['pos'] < stretch[0]]
    after = [s for s in mh_seg if stretch and s['pos'] > stretch[1]]
    bounds = [before[-1]['pos'] if before else None, stretch[0] if stretch else None, stretch[1] if stretch else None, after[0]['pos'] if after else None]
    check('bounds', numbers.get('bounds'), bounds)
    patchy_lo = differ_outside[0] if differ_outside else (stretch[0] if stretch else None)
    for s in sites:
        s['in_stretch'] = bool(stretch and stretch[0] <= s['pos'] <= stretch[1])
    # composition: at her phased mother-only sites, how many of the child's maternal copies carry her haplotype 1
    compo = {'left': collections.Counter(), 'inside': collections.Counter(), 'right': collections.Counter()}
    pw = th['pair_window']
    block_sp = None
    if stretch:
        for b in blocks:
            if b['start'] <= stretch[0] and b['end'] >= stretch[1]:
                block_sp = b
        for s in mh_seg:
            if 'own_on_hap' not in s or s['copies'] is None:
                continue
            if block_sp and not (block_sp['start'] <= s['pos'] <= block_sp['end']):
                continue
            if not block_sp and not (stretch[0] - pw <= s['pos'] <= stretch[1] + pw):
                continue
            k = s['copies'] if s['own_on_hap'] == 1 else 2 - s['copies']
            s['hap1_copies'] = k
            where = 'left' if s['pos'] < stretch[0] else ('right' if s['pos'] > stretch[1] else 'inside')
            compo[where][str(k)] += 1
    for w in ('left', 'inside', 'right'):
        check(f'compo.{w}', numbers.get('compo', {}).get(w), dict(compo[w]))
    check('compo.block', numbers.get('compo', {}).get('block'), [block_sp['start'], block_sp['end'], block_sp['n']] if block_sp else None)

    # reads across each end of the stretch: which of her haplotypes they carry before and after
    def hap_of(rec, lo, hi):
        """The haplotype a read carries at her phased sites in [lo, hi] (the block of the stretch only):
        two or more sites, and two in three agreeing."""
        c = collections.Counter()
        for pos, base in rec['alleles'].items():
            s = site_by_pos[pos]
            if lo <= pos <= hi and 'ref_on_hap' in s and (block_sp is None or s.get('block') == blocks.index(block_sp)):
                hap = s['ref_on_hap'] if base == s['ref'] else 3 - s['ref_on_hap']
                c[hap] += 1
        n = sum(c.values())
        if n < 2:
            return None
        (h1, n1), = c.most_common(1)
        return h1 if n1 >= 2 * n / 3 else None

    switch = {'left': collections.Counter(), 'right': collections.Counter()}
    switch_reads = {'left': [], 'right': []}
    if stretch and all(bounds):
        for end, (lo, hi) in (('left', (bounds[0], bounds[1])), ('right', (bounds[2], bounds[3]))):
            for name, rec in reads.items():
                if not (rec['start'] <= lo and rec['end'] >= hi):
                    continue
                if not any(a.mapq >= MAPQ and a.pos <= lo and a.end >= hi for a in hifi[name]):
                    continue
                hb, ha = hap_of(rec, rec['start'], lo - 1), hap_of(rec, hi + 1, rec['end'])
                if rec['parent'] == 'father' or (rec['own_pat'] >= 1 and rec['own_mat'] == 0 and rec['oh_mat'] == 0):
                    key = 'dad'
                elif hb is None or ha is None:
                    continue
                else:
                    key = f'{hb}>{ha}'
                switch[end][key] += 1
                switch_reads[end].append((name, key))
    check('switch_left', numbers.get('switch_left'), dict(switch['left']))
    check('switch_right', numbers.get('switch_right'), dict(switch['right']))

    # ---------------------------------------------------------------------------------------- 8. the assembly's pieces
    log('the assembly')
    ref_seq = read_fasta(P(cfg['files']['ref']))[cfg['files']['ref_name']]
    ref_offset = WIN[0]
    asm = cfg['assembly']
    pieces = {}
    extra_records = []
    for hap, bam in cfg['files']['asm_aln'].items():
        recs = [PieceAln(Aln(line.split('\t'), keep_seq=True), hap) for line in stream(['samtools', 'view', P(bam)])]
        recs = [r for r in recs if r.re >= S and r.rs <= E]
        recs.sort(key=lambda r: (r.contig, r.qs))
        names = asm['piece_names'][hap]
        if isinstance(names, str):
            for i, r in enumerate(recs):
                pieces[names[i]] = r
        else:
            for r in recs:
                pieces[names[r.contig]] = r
    for letter, r in pieces.items():
        r.letter = letter
        seqinfo = next(v for v in asm['sequences'].values() if v['contig'] == r.contig)
        r.parent_hap = seqinfo['parent']
        r.local_name = next(k for k, v in asm['sequences'].items() if v['contig'] == r.contig)
        r.local_offset = seqinfo['local_offset']
        r.variants_snv, r.variants_ins, r.variants_del = r.variants(ref_seq, ref_offset, max(r.rs, WIN[0]), min(r.re, WIN[1]))
        # the piece's parent by the trio's alleles inside the segment
        mom = dad = 0
        for s in sites:
            if not s['in_segment'] or not (r.rs <= s['pos'] <= r.re):
                continue
            base = r.variants_snv[s['pos']][1] if s['pos'] in r.variants_snv else s['ref']
            if s['pos'] in r.variants_del:
                continue
            if s['kind'] == 'OH':
                if base == s['mother_allele']:
                    mom += 1
                elif base == s['father_allele']:
                    dad += 1
            elif s['kind'] == 'MH' and base == s['mother_allele']:
                mom += 1
            elif s['kind'] == 'FH' and base == s['father_allele']:
                dad += 1
        r.mom_own, r.dad_own = mom, dad
        mn, oth = th['parent_piece_min'], th['parent_piece_other']
        r.parent = 'maternal' if mom >= mn and dad <= oth * (mom + dad) else ('paternal' if dad >= mn and mom <= oth * (mom + dad) else 'unknown')
        exp = numbers.get('pieces', {}).get(letter, {})
        for k, v in (('rs', r.rs), ('re', r.re), ('qs', r.qs), ('qe', r.qe), ('mom_own', mom), ('dad_own', dad), ('nm', r.nm),
                     ('aligned', r.aligned), ('parent', r.parent), ('kind', r.kind), ('contig', r.contig)):
            check(f'pieces.{letter}.{k}', exp.get(k), v)
    M = pieces[asm['maternal_chromosome']]
    m_segment = [M.ref_to_contig(S) - M.local_offset, M.ref_to_contig(E) - M.local_offset]
    check('realign.m_segment', numbers.get('realign', {}).get('m_segment'), m_segment)
    # the second copy's pieces and the missing copy
    second = [r for r in pieces.values() if r.parent == 'maternal' and r.letter != M.letter]
    cov = sorted((max(r.rs, S), min(r.re, E)) for r in second)
    holes, cur = [], S
    for a, b in cov:
        if a > cur:
            holes.append([cur, a - 1])
        cur = max(cur, b + 1)
    if cur <= E:
        holes.append([cur, E])
    check('holes', numbers.get('holes'), holes)
    check('holes_total', sum(b - a + 1 for a, b in numbers.get('holes', [])) or None, sum(b - a + 1 for a, b in holes))
    # copies the assembly holds per 10-kb bin over the segment
    asm_copies = collections.Counter()
    for row in depth_rows:
        b0, b1 = row['bin0'], row['bin1']
        n = sum(1 for r in pieces.values() if r.rs <= b0 + 5000 <= r.re)
        row['asm_copies'] = n
        if overlaps(b0, b1):
            asm_copies[str(n)] += 1
    check('copies', numbers.get('copies'), dict(asm_copies))
    # joins: consecutive pieces of a scaffold across its gaps
    gap_rows = [row for path in cfg['files']['qc']['gaps'] for row in load_bed_rows(P(path))]
    joins = []
    scaffold = sorted([r for r in pieces.values() if r.local_name == 'hap1'], key=lambda r: r.qs)
    for left, right in zip(scaffold, scaffold[1:]):
        gaps_between = [(g[1] + 1, g[2]) for g in gap_rows if g[0] == left.contig and left.qe < g[1] + 1 and g[2] < right.qs]
        gap = gaps_between[0] if gaps_between else None
        p1, p2 = left.re, right.rs
        through1 = sum(1 for r in reads.values() if any(a.mapq >= MAPQ and a.through(p1, flank) for a in hifi[r['name']]))
        through2 = sum(1 for r in reads.values() if any(a.mapq >= MAPQ and a.through(p2, flank) for a in hifi[r['name']]))
        carry = stop = 0
        gs = th['gap_slack']
        for name, alns in hifi.items():
            ends = [a for a in alns if a.mapq >= MAPQ and abs(a.end - p1) <= gs and a.rclip >= minclip]
            starts = [a for a in alns if a.mapq >= MAPQ and abs(a.pos - p2) <= gs and a.lclip >= minclip]
            if ends and starts and any(a.strand == b.strand and a is not b for a in ends for b in starts):
                carry += 1
            elif any(a.mapq >= MAPQ and ((a.rclip >= minclip and min(abs(a.end - p1), abs(a.end - p2)) <= gs) or
                                          (a.lclip >= minclip and min(abs(a.pos - p1), abs(a.pos - p2)) <= gs)) for a in alns):
                stop += 1
        joins.append({'left': left.letter, 'right': right.letter, 'p1': p1, 'p2': p2, 'gap': list(gap) if gap else None,
                      'unaligned': (right.qs - left.qe - 1), 'back': p1 - p2, 'carry': carry, 'through1': through1,
                      'through2': through2, 'stop': stop})
    for i, j in enumerate(joins):
        exp = numbers.get('joins', [{}] * len(joins))[i] if i < len(numbers.get('joins', [])) else {}
        for k in ('left', 'right', 'p1', 'p2', 'gap', 'unaligned', 'back', 'carry', 'through1', 'through2', 'stop'):
            check(f'joins.{i}.{k}', exp.get(k), j[k])
    check('gap_count', numbers.get('gap_count'), sum(1 for g in gap_rows if g[0] == scaffold[0].contig and scaffold[0].qs <= g[1] <= scaffold[-1].qe) if scaffold else None)
    # pieces against each other where both run, over the segment and 20 kb either side
    lo_w, hi_w = S - pw, E + pw
    where_of = lambda p: 'stretch' if stretch and stretch[0] <= p <= stretch[1] else ('patchy' if patchy_lo and stretch and patchy_lo <= p < stretch[0] else 'elsewhere')
    piece_pairs = []
    order = list(pieces)
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            x, y = pieces[order[i]], pieces[order[j]]
            lo, hi = max(x.rs, y.rs, lo_w), min(x.re, y.re, hi_w)
            if lo > hi:
                continue
            recs_x = {('snv', p, v) for p, v in x.variants_snv.items()} | {('ins', p, v) for p, v in x.variants_ins.items()} | {('del', p, v) for p, v in x.variants_del.items()}
            recs_y = {('snv', p, v) for p, v in y.variants_snv.items()} | {('ins', p, v) for p, v in y.variants_ins.items()} | {('del', p, v) for p, v in y.variants_del.items()}
            diff = [d for d in recs_x ^ recs_y if lo <= d[1] <= hi]
            n_snv = sum(1 for d in diff if d[0] == 'snv')
            n_indel = len(diff) - n_snv
            allp = sorted(d[1] for d in diff)
            piece_pairs.append({'x': x.letter, 'y': y.letter, 'lo': lo, 'hi': hi, 'len': hi - lo + 1, 'snv': n_snv,
                                'indel': n_indel, 'per_kb': len(diff) / ((hi - lo + 1) / 1000),
                                'where': dict(collections.Counter(where_of(p) for p in allp)), 'first': allp[0] if allp else None,
                                'last': allp[-1] if allp else None, 'same_scaffold': x.contig == y.contig,
                                'parents': [x.parent, y.parent]})
    exp_pairs = {(p['x'], p['y']): p for p in numbers.get('piece_pairs', [])}
    for p in piece_pairs:
        e = exp_pairs.get((p['x'], p['y'])) or exp_pairs.get((p['y'], p['x']))
        if e is None:
            continue
        for k in ('lo', 'hi', 'len', 'snv', 'indel', 'first', 'last', 'same_scaffold'):
            check(f'piece_pairs.{p["x"]}{p["y"]}.{k}', e.get(k), p[k])
        check(f'piece_pairs.{p["x"]}{p["y"]}.where', e.get('where'), p['where'])
    # the missing copy's reads: which piece they match, over each hole inside the segment
    missing = []
    for lo, hi in holes:
        if hi - lo < 10000:
            continue
        tally = collections.Counter()
        for s in sites:
            if not (lo <= s['pos'] <= hi) or s['kind'] not in ('OH', 'MH', 'FH'):
                continue
            # the missing copy is maternal: at a parent-opposite site it carries her allele; at a father-only site it lacks his;
            # at a mother-only site her own allele is on k of the child's three copies, and M holds one maternal copy
            m_base = M.variants_snv[s['pos']][1] if s['pos'] in M.variants_snv else s['ref']
            others = [r for r in pieces.values() if r.letter != M.letter and r.rs <= s['pos'] <= r.re and r.parent == 'paternal']
            if s['kind'] == 'OH':
                allele, m_has, missing_has = s['mother_allele'], m_base == s['mother_allele'], True
            elif s['kind'] == 'FH':
                allele, m_has, missing_has = s['father_allele'], m_base == s['father_allele'], False
            else:
                k = s.get('copies')
                if k not in (0, 1, 2):
                    continue
                allele = s['mother_allele']
                m_has = m_base == allele
                missing_has = (k - (1 if m_has else 0)) >= 1
            tally['M same' if missing_has == m_has else 'M differs'] += 1
            for r in others:
                base = r.variants_snv[s['pos']][1] if s['pos'] in r.variants_snv else s['ref']
                tally[f'{r.letter} same' if (base == allele) == missing_has else f'{r.letter} differs'] += 1
        missing.append({'lo': lo, 'hi': hi, 'tally': dict(tally)})
    for i, mrow in enumerate(missing):
        exp = numbers.get('missing', [{}] * len(missing))[i] if i < len(numbers.get('missing', [])) else {}
        check(f'missing.{i}.tally', exp.get('tally'), mrow['tally'])
    # the QC calls, placed on GRCh38 through the piece that holds their start
    flags_placed = []
    qc_rows = []
    for track in ('nucflag', 'flagger_hifi', 'flagger_ont'):
        for path in cfg['files']['qc'][track]:
            for contig, b0, b1, label in load_bed_rows(P(path)):
                qc_rows.append((track, contig, b0, b1, label))
    for track, contig, b0, b1, label in qc_rows:
        holders = [r for r in pieces.values() if r.contig == contig and r.qs <= b0 + 1 <= r.qe]
        if not holders:
            holders = [r for r in pieces.values() if r.contig == contig and r.qs <= b1 and b0 + 1 <= r.qe]
            if not holders:
                continue
        r = holders[0]
        if track == 'nucflag' and label != 'HET':
            flags_placed.append({'track': track, 'contig': contig, 'qs': b0 + 1, 'qe': b1, 'label': label, 'piece': r.letter,
                                 'rs': r.contig_to_ref(b0 + 1), 're': r.contig_to_ref(min(b1, r.qe))})
    exp_fl = numbers.get('flags_placed', [])
    for i, fp in enumerate(flags_placed):
        e = exp_fl[i] if i < len(exp_fl) else {}
        for k in ('contig', 'qs', 'qe', 'label', 'piece', 'rs', 're'):
            check(f'flags_placed.{i}.{k}', e.get(k), fp[k])
    # uniqueness: the 240 pieces of the segment against all of GRCh38
    paf_best, paf_other, paf_bp = {}, collections.Counter(), {'before the start': [], 'after the start': [], 'before the end': [], 'after the end': []}
    bp_pieces = {}
    with open(P(cfg['files']['pieces_paf'])) as f:
        for line in f:
            fl = line.rstrip('\n').split('\t')
            q, qlen, t, ts, te, matches, blen, mq = fl[0], int(fl[1]), fl[5], int(fl[7]), int(fl[8]), int(fl[9]), int(fl[10]), int(fl[11])
            qspan = int(fl[3]) - int(fl[2])
            m = re.match(r'(\w+):(\d+)-(\d+)', q)
            qc, qs_, qe_ = m.group(1), int(m.group(2)), int(m.group(3))
            own = t == qc and ts + 1 <= qe_ and te >= qs_
            entry = paf_best.setdefault(q, {'own_first': None, 'other': []})
            if entry['own_first'] is None:
                entry['own_first'] = (own, mq)
            if not own and blen >= 1000 and matches / blen >= 0.9:
                entry['other'].append((t, ts + 1, te, qspan, round(matches / blen, 3)))
    self_pieces = len(paf_best)
    self_first = sum(1 for e in paf_best.values() if e['own_first'] and e['own_first'][0] and e['own_first'][1] == 60)
    self_other = sum(1 for e in paf_best.values() if e['other'])
    check('self_pieces', numbers.get('self_pieces'), self_pieces)
    check('self_first', numbers.get('self_first'), self_first)
    check('self_other', numbers.get('self_other'), self_other)
    for q, e in paf_best.items():
        m = re.match(r'(\w+):(\d+)-(\d+)', q)
        qs_, qe_ = int(m.group(2)), int(m.group(3))
        key = None
        if qe_ == S - 1:
            key = 'before the start'
        elif qs_ == S and qe_ == S + 4999:
            key = 'after the start'
        elif qe_ == E and qs_ == E - 4999:
            key = 'before the end'
        elif qs_ == E + 1:
            key = 'after the end'
        if key and e['other']:
            best_o = max(e['other'], key=lambda o: o[3])
            paf_bp[key].append([best_o[0], best_o[1], best_o[2], best_o[3], best_o[4]])
    check('self_bp', numbers.get('self_bp'), paf_bp)
    self_rows = []
    for q, e in paf_best.items():
        if e['other']:
            best_o = max(e['other'], key=lambda o: o[3])
            self_rows.append({'piece': q, 'matches': len(e['other']), 'longest': f'{best_o[0]}:{best_o[1]}-{best_o[2]}', 'bp': best_o[3], 'identity': best_o[4]})
    # microhomology at the junction
    mh = 0
    while mh < 50 and ref_seq[E - ref_offset - mh].upper() == ref_seq[S - ref_offset - 1 - mh].upper():
        mh += 1
    check('mh', numbers.get('mh'), mh)

    # ---------------------------------------------------------------------------------------- 9. the reads given back to the assembly
    log('the reads given back to the assembly')
    realign_bam = P(cfg['files']['realign'])
    local = read_fasta(P(cfg['files']['asm_local']))
    local_len = {k: len(v) for k, v in local.items()}
    rl = group_alignments(stream(['samtools', 'view', '-F', '0x704', realign_bam]))
    # pieces in local coordinates
    piece_local = {}
    for r in pieces.values():
        # clipped to the local slice: A and E run past it
        piece_local[r.letter] = (r.local_name, max(1, r.qs - r.local_offset), min(local_len[r.local_name], r.qe - r.local_offset))
    gaps_local = []
    for contig, b0, b1, _ in gap_rows:
        for name, info in asm['sequences'].items():
            if info['contig'] == contig:
                a, b = b0 + 1 - info['local_offset'], b1 - info['local_offset']
                if 1 <= a <= local_len[name]:
                    gaps_local.append((name, a, b))
    gaps_local.sort()

    def piece_at(seq, pos):
        best = None
        for letter, (sname, a, b) in piece_local.items():
            if sname == seq and a <= pos <= b:
                return letter
        # in a gap or unaligned bases: the next piece to the right
        cands = [(a, letter) for letter, (sname, a, b) in piece_local.items() if sname == seq and a > pos]
        if cands:
            return min(cands)[1]
        cands = [(b, letter) for letter, (sname, a, b) in piece_local.items() if sname == seq]
        return max(cands)[1] if cands else None

    m_local = M.local_name
    ms, me = m_segment
    per_piece = collections.defaultdict(collections.Counter)
    jr_at_start = jr_at_end = jr_split = jr_through = 0
    for name, alns in rl.items():
        rec = reads.get(name)
        if rec is None:
            continue
        prim = [a for a in alns if a.primary]
        if not prim:
            continue
        p = prim[0]
        letter = piece_at(p.chrom, p.pos)
        rec['asm_piece'] = letter
        rec['asm_mapq'] = p.mapq
        rec['asm_seq'] = p.chrom
        rec['asm_start'] = p.pos
        rec['asm_end'] = p.end
        lab = {'mother': 'mat', 'father': 'pat', 'unphased': '?'}[rec['parent']]
        per_piece[letter]['reads'] += 1
        per_piece[letter][lab] += 1
        if p.mapq == 0:
            per_piece[letter]['q0'] += 1
        at_end = any(a.chrom == m_local and abs(a.end - me) <= slack and a.rclip >= minclip for a in alns)
        at_start = any(a.chrom == m_local and abs(a.pos - ms) <= slack and a.lclip >= minclip for a in alns)
        thr = any(a.chrom == m_local and (a.through(ms, flank) or a.through(me, flank)) for a in alns)
        rec['asm_split_end'] = at_end
        rec['asm_split_start'] = at_start
        rec['asm_through_M'] = thr
        if rec['junction']:
            jr_at_end += at_end
            jr_at_start += at_start
            jr_split += at_end and at_start
            jr_through += thr
    for letter, exp in numbers.get('realign', {}).get('per_piece', {}).items():
        for k, v in exp.items():
            check(f'realign.per_piece.{letter}.{k}', v, per_piece[letter][k])
    check('realign.reads', numbers.get('realign', {}).get('reads'), sum(c['reads'] for c in per_piece.values()))
    check('realign.jr_split_at_m', numbers.get('realign', {}).get('jr_split_at_m'), jr_split)
    check('realign.jr_at_end', numbers.get('realign', {}).get('jr_at_end'), jr_at_end)
    check('realign.jr_at_start', numbers.get('realign', {}).get('jr_at_start'), jr_at_start)
    check('realign.jr_through', numbers.get('realign', {}).get('jr_through'), jr_through)
    # depth on the assembly in 5-kb bins by parent, and at MAPQ 0
    BIN = 5000
    adepth = {seq: collections.defaultdict(lambda: collections.Counter()) for seq in local_len}
    base_depth = {seq: collections.Counter() for seq in local_len}  # per-base over 2-kb edges only (filled below)
    edges = []
    for seq, a, b in gaps_local:
        edges.append((seq, a - 2000, a - 1, 'left'))
        edges.append((seq, b + 1, b + 2000, 'right'))
    for name, alns in rl.items():
        rec = reads.get(name, {'parent': 'unphased'})
        lab = {'mother': 'mat', 'father': 'pat', 'unphased': '?'}[rec['parent']]
        for a in alns:
            d = adepth[a.chrom]
            b0 = (a.pos - 1) // BIN
            b1 = (a.end - 1) // BIN
            for b in range(b0, b1 + 1):
                lo, hi = max(a.pos, b * BIN + 1), min(a.end, (b + 1) * BIN)
                n = hi - lo + 1
                d[b][lab] += n
                if a.mapq == 0:
                    d[b]['q0'] += n
            for seq, elo, ehi, side in edges:
                if a.chrom == seq and a.end >= elo and a.pos <= ehi:
                    base_depth[seq][(elo, ehi)] += min(a.end, ehi) - max(a.pos, elo) + 1
    # one haplotype's depth: the left flank beyond 10 kb of the segment, on each chromosome
    flank_depth = {}
    for r in (pieces[asm['piece_names']['hap1'][0]], M):
        q0 = r.ref_to_contig(LWIN[0] + 20000) - r.local_offset
        q1 = r.ref_to_contig(S - 10000) - r.local_offset
        b0, b1 = (q0 - 1) // BIN + 1, (q1 - 1) // BIN
        tot = sum(sum(v for k, v in adepth[r.local_name][b].items() if k != 'q0') for b in range(b0, b1))
        flank_depth[r.local_name] = tot / ((b1 - b0) * BIN) if b1 > b0 else None
    one_copy = mean([v for v in flank_depth.values() if v])
    check('realign.one_copy', numbers.get('realign', {}).get('one_copy'), round(one_copy, 2))
    check('realign.flank_pat', numbers.get('realign', {}).get('flank_pat'), round(flank_depth.get('hap1', 0), 2))
    check('realign.flank_mat', numbers.get('realign', {}).get('flank_mat'), round(flank_depth.get('hap2', 0), 2))
    asm_depth_rows = []
    for seq in ('hap1', 'hap2', 'U'):
        if seq not in local_len:
            continue
        for b in range((local_len[seq] - 1) // BIN + 1):
            d = adepth[seq][b]
            n = min(BIN, local_len[seq] - b * BIN)
            asm_depth_rows.append({'seq': seq, 'bin0': b * BIN, 'bin1': b * BIN + n, 'mat': d['mat'] / n / one_copy,
                                   'pat': d['pat'] / n / one_copy, 'unphased': d['?'] / n / one_copy, 'q0': d['q0'] / n / one_copy,
                                   'piece': piece_at(seq, b * BIN + n // 2)})
    # the joins on the assembly: where the reads straight through each point align, crossing, edge depth
    rjoins = []
    for i, j in enumerate(joins):
        seq, a, b = gaps_local[i] if i < len(gaps_local) else (None, None, None)
        entry = {'left': j['left'], 'right': j['right'], 'lend': a - 1 if a else None, 'rbeg': None, 'cross': 0, 'stop_left': 0, 'stop_right': 0}
        if a:
            rb = piece_local[j['right']][1]
            entry['rbeg'] = rb
            for name, alns in rl.items():
                ends = [x for x in alns if x.chrom == seq and abs(x.end - (a - 1)) <= th['gap_slack'] and x.rclip >= minclip]
                starts = [x for x in alns if x.chrom == seq and abs(x.pos - rb) <= th['gap_slack'] and x.lclip >= minclip]
                if ends and starts:
                    entry['cross'] += 1
                elif ends:
                    entry['stop_left'] += 1
                elif starts:
                    entry['stop_right'] += 1
            entry['edge_left'] = base_depth[seq][(a - 2000, a - 1)] / 2000 / one_copy
            entry['edge_right'] = base_depth[seq][(b + 1, b + 2000)] / 2000 / one_copy
        for key, point in (('p1', j['p1']), ('p2', j['p2'])):
            where = collections.Counter()
            n = 0
            for name, rec in reads.items():
                if any(x.mapq >= MAPQ and x.through(point, flank) for x in hifi[name]) and rec.get('asm_piece'):
                    n += 1
                    lab = {'mother': 'mat', 'father': 'pat', 'unphased': '?'}[rec['parent']]
                    where[f'{rec["asm_piece"]}|{lab}'] += 1
            entry[key] = {'n': n, 'where': dict(where)}
        rjoins.append(entry)
    for i, e in enumerate(rjoins):
        exp = numbers.get('realign', {}).get('joins', [])
        exp = exp[i] if i < len(exp) else {}
        for k in ('lend', 'rbeg', 'cross', 'stop_left', 'stop_right', 'p1', 'p2'):
            check(f'realign.joins.{i}.{k}', exp.get(k), e[k])
        for k in ('edge_left', 'edge_right'):
            check(f'realign.joins.{i}.{k}', exp.get(k), e.get(k))
    # the local sequences: N count and record agreement
    for seq, info in asm['sequences'].items():
        if seq in local:
            check(f'realign.seqs.{seq}.len', numbers.get('realign', {}).get('seqs', {}).get(seq, {}).get('len'), local_len[seq])
            check(f'realign.seqs.{seq}.N', numbers.get('realign', {}).get('seqs', {}).get(seq, {}).get('N'), local[seq].upper().count('N'))

    # ---------------------------------------------------------------------------------------- 10. the call sets
    log('the call sets')
    calls = []
    yale_summary = {}
    for caller, (path, made_from) in cfg['files']['calls_yale'].items():
        samples, recs = parse_vcf(P(path))
        found, bnd_near, n50, longest = [], 0, 0, 0
        for rec in recs:
            svtype, pos, end, svlen, mate = sv_of(rec)
            gt = rec['gts'][0] if rec['gts'] else '.'
            if svtype not in ('BND', 'SNV') and svlen is not None and LWIN[0] <= pos <= LWIN[1]:
                n50 += svlen >= 50
                longest = max(longest, svlen)
            if svtype == 'BND':
                near_s, near_e = abs(pos - S) <= 5000, abs(pos - E) <= 5000
                same = mate is not None and mate[0] == chrom
                mate_near = same and (abs(mate[1] - S) <= 5000 or abs(mate[1] - E) <= 5000)
                inside = same and S - 5000 <= pos <= E + 5000 and S - 5000 <= mate[1] <= E + 5000 and abs(mate[1] - pos) >= 10000
                if (near_s or near_e) and mate_near:
                    reflects = 'the junction' if (near_s and abs(mate[1] - E) <= 5000) or (near_e and abs(mate[1] - S) <= 5000) else 'other'
                elif inside:
                    # a breakend pair within the segment, 10 kb or more apart: from the assembly, a step of the scaffold
                    reflects = 'the scaffold' if made_from == 'the assembly' else 'other'
                else:
                    bnd_near += near_s or near_e
                    continue
                found.append({'set': caller, 'made_from': made_from, 'type': 'BND', 'alt': rec['alt'], 'pos': pos, 'end': mate[1], 'len': abs(mate[1] - pos),
                              'gt': gt, 'id': rec['id'], 'filter': rec['filter'], 'qual': rec['qual'], 'support': support_of(rec), 'reflects': reflects})
                continue
            if svlen is None or svlen < 10000:
                continue
            if end < S - 5000 or pos > E + 5000 or pos < LWIN[0] or end > LWIN[1]:
                continue
            if abs(pos - S) <= 100 and abs(end - E) <= 100 and svtype == 'DUP':
                reflects = 'the junction'
            else:
                reflects = 'the scaffold' if made_from == 'the assembly' else 'other'
            found.append({'set': caller, 'made_from': made_from, 'type': svtype, 'alt': rec['alt'] if rec['alt'].startswith('<') else f'<{svtype}>',
                          'pos': pos, 'end': end, 'len': svlen, 'gt': gt, 'id': rec['id'], 'filter': rec['filter'], 'qual': rec['qual'],
                          'support': support_of(rec), 'reflects': reflects})
        calls.extend(found)
        yale_summary[caller] = {'made_from': made_from, 'n': len(found), 'bnd_near': bnd_near, 'n50': n50, 'longest': longest, 'records': len(recs),
                                'verdict': 'silent' if not found else ('calls it' if any(f['reflects'] == 'the junction' for f in found) else 'calls the scaffold')}
    published = []
    for setname, path, jpath in cfg['files']['calls_published']:
        meta = json.load(open(P(jpath)))
        samples, recs = parse_vcf(P(path))
        held = meta.get('held', [])
        n = longest = n_big = 0
        big = []
        for rec in recs:
            svtype, pos, end, svlen, mate = sv_of(rec)
            if svtype in ('SNV', None) or svlen is None:
                continue
            svlen = carried_len(rec, svlen)
            carried = any(a.isdigit() and int(a) > 0 for g in rec['gts'] for a in g.replace('|', '/').split('/')) if rec['gts'] else False
            if not held:
                carried = False
            if carried and svlen >= 50 and end >= S and pos <= E:
                n += 1
                longest = max(longest, svlen)
                if svlen >= 10000:
                    big.append(rec)
            if svlen >= 10000 and end >= S - 5000 and pos <= E + 5000 and pos >= LWIN[0] and end <= LWIN[1]:
                n_big += 1
                gts = ';'.join(f'{s}={g}' for s, g in zip(samples, rec['gts'])) if rec['gts'] else ''
                calls.append({'set': setname, 'made_from': 'published', 'type': svtype, 'alt': rec['alt'] if rec['alt'].startswith('<') else f'<{svtype}>',
                              'pos': pos, 'end': end, 'len': svlen, 'gt': gts, 'id': rec['id'], 'filter': rec['filter'], 'qual': rec['qual'],
                              'support': support_of(rec), 'reflects': 'the junction' if abs(pos - S) <= 100 and abs(end - E) <= 100 else 'other'})
        published.append({'set': setname, 'file': os.path.basename(meta.get('url', path)), 'held': held, 'n_sv50': n, 'longest': longest,
                          'n_big': n_big, 'records': len(recs), 'columns': meta.get('columns')})
    samples, recs = parse_vcf(P(cfg['files']['event_all_samples']))
    carriers = []
    for rec in recs:
        for s_, g in zip(samples, rec['gts']):
            if '1' in g.replace('|', '/').split('/'):
                carriers.append(s_)
    check('carriers', numbers.get('carriers'), sorted(set(carriers)))

    # ---------------------------------------------------------------------------------------- 11. the junction as called (optional alignment)
    junction_ref_rel = cfg['files'].get('junction_ref')
    if args.align and junction_ref_rel:
        log('the junction as called')
        L = 30000
        left = ref_seq[E - L + 1 - ref_offset:E - ref_offset + 1]
        right = ref_seq[S - ref_offset:S + L - ref_offset]
        name = f'junction_{chrom}_{E - L + 1}-{E}_{S}-{S + L - 1}'
        jpath = P(junction_ref_rel)
        with open(jpath, 'w') as f:
            f.write(f'>{name}\n')
            seq = left + right
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + '\n')
        run(['samtools', 'faidx', jpath])
        reads_fa = os.path.join(cache, 'hifi_reads.fa')
        with open(reads_fa, 'w') as f:
            f.write(run(['samtools', 'fasta', '-F', '0x904', hifi_bam], stderr=subprocess.DEVNULL))
        sam = os.path.join(cache, 'on_junction.sam')
        with open(sam, 'w') as f:
            subprocess.run(['minimap2', '-ax', 'map-hifi', '--eqx', '-Y', '-t', '8', jpath, reads_fa], check=True, stdout=f, stderr=subprocess.DEVNULL)
        tmp_bam = os.path.join(cache, 'on_junction.tmp.bam')
        run(['samtools', 'sort', '-o', tmp_bam, sam])
        # a clean header: no @PG, no UR
        hdr = [l for l in run(['samtools', 'view', '-H', tmp_bam]).splitlines() if not l.startswith('@PG')]
        hdr = ['\t'.join(t for t in l.split('\t') if not t.startswith('UR:')) for l in hdr]
        hdr_path = os.path.join(cache, 'on_junction.header.sam')
        with open(hdr_path, 'w') as f:
            f.write('\n'.join(hdr) + '\n')
        out_bam = P(cfg['files']['hifi_on_junction'])
        mapped_bam = os.path.join(cache, 'on_junction.mapped.bam')
        with open(mapped_bam, 'wb') as f:
            subprocess.run(['samtools', 'view', '--no-PG', '-b', '-F', '4', tmp_bam], stdout=f, check=True)
        with open(out_bam, 'wb') as f:
            subprocess.run(['samtools', 'reheader', '-P', hdr_path, mapped_bam], stdout=f, check=True)
        run(['samtools', 'index', out_bam])
        for pth in (sam, tmp_bam, mapped_bam):
            os.remove(pth)
        n_junction_full = 0
        for line in stream(['samtools', 'view', '-F', '0x904', out_bam]):
            f_ = line.split('\t')
            a = Aln(f_)
            if a.name in jr and a.pos <= L - 500 and a.end >= L + 500:
                n_junction_full += 1
        log(f'  {n_junction_full} of the {len(jr)} junction reads span the join point unbroken on the junction as called')
        checks.append({'key': 'junction_as_called.jr_unbroken', 'expected': len(jr), 'got': n_junction_full, 'ok': n_junction_full == len(jr), 'note': 'the junction reads on the junction as called, 500 bp either side'})

    # ---------------------------------------------------------------------------------------- 11b. the assembly against itself
    # hap1 and U aligned to hap2 directly (minimap2 asm5), no reference involved: where the paternal scaffold holds the
    # maternal chromosome's sequence, base for base
    # each gap-free piece of the paternal scaffold, and U, is aligned on its own, so that no chain runs across a gap
    log('the assembly against itself')
    q_fa, t_fa = os.path.join(cache, 'asm_query.fa'), os.path.join(cache, 'asm_target.fa')
    piece_offset = {}
    with open(q_fa, 'w') as f:
        for letter, (sname, a, b) in piece_local.items():
            if sname == m_local:
                continue
            f.write(f'>{letter}\n{local[sname][a - 1:b]}\n')
            piece_offset[letter] = (sname, a - 1)
    with open(t_fa, 'w') as f:
        f.write(f'>{m_local}\n{local[m_local]}\n')
    asm_vs_asm = []
    paf_text = run(['minimap2', '-x', 'asm20', '-c', '--secondary=no', '-t', '8', t_fa, q_fa], stderr=subprocess.DEVNULL)
    with open(os.path.join(out, 'asm_vs_asm.paf'), 'w') as f:
        f.write('# each piece of hap1 (its name) and U against hap2, minimap2 -x asm20 -c; query coordinates are the piece\'s own\n')
        f.write(paf_text)
    for line in paf_text.splitlines():
        fl = line.split('\t')
        q, qs_, qe_, strand, t, ts_, te_, matches, blen, mq = fl[0], int(fl[2]), int(fl[3]), fl[4], fl[5], int(fl[7]), int(fl[8]), int(fl[9]), int(fl[10]), int(fl[11])
        if blen < 5000:
            continue
        sname, off = piece_offset[q]
        asm_vs_asm.append({'piece': q, 'q': sname, 'qs': off + qs_ + 1, 'qe': off + qe_, 'strand': strand, 't': t, 'ts': ts_ + 1, 'te': te_,
                           'matches': matches, 'blen': blen, 'identity': round(matches / blen, 5), 'mapq': mq})

    # ---------------------------------------------------------------------------------------- 12. the reference for the read panels
    # chr2 as GRCh38 has it, N outside the window, bgzipped and indexed: igv.js can then show the reads without any
    # server but this repository's (GitHub Pages or serve.py)
    ref_gz_rel = cfg['files'].get('ref_padded')
    if ref_gz_rel:
        log('the padded reference')
        chrom_len = None
        for line in run(['samtools', 'view', '-H', hifi_bam]).splitlines():
            if line.startswith('@SQ') and f'SN:{chrom}\t' in line + '\t':
                chrom_len = int([t for t in line.split('\t') if t.startswith('LN:')][0][3:])
        if chrom_len is None:
            raise RuntimeError(f'{chrom} not in the HiFi BAM header')
        gz = P(ref_gz_rel)
        os.makedirs(os.path.dirname(gz), exist_ok=True)
        with open(gz, 'wb') as fout:
            bg = subprocess.Popen(['bgzip', '-c'], stdin=subprocess.PIPE, stdout=fout)
            w = bg.stdin
            w.write(f'>{chrom}\n'.encode())
            width = 60
            pieces_iter = [('N', WIN[0] - 1), (ref_seq.upper(), None), ('N', chrom_len - WIN[1])]
            carry = ''
            for base, count in pieces_iter:
                chunk_src = (base * count) if count is not None else base
                data = carry + chunk_src
                nlines = len(data) // width
                for i in range(0, nlines * width, 6000000):
                    block = data[i:min(nlines * width, i + 6000000)]
                    w.write(('\n'.join(block[j:j + width] for j in range(0, len(block), width)) + '\n').encode())
                carry = data[nlines * width:]
            if carry:
                w.write((carry + '\n').encode())
            w.close()
            bg.wait()
        for ext in ('.fai', '.gzi'):
            if os.path.exists(gz + ext):
                os.remove(gz + ext)
        run(['samtools', 'faidx', gz])
        # the window reads back identical
        got = run(['samtools', 'faidx', gz, f'{chrom}:{WIN[0]}-{WIN[1]}']).split('\n', 1)[1].replace('\n', '')
        checks.append({'key': 'ref_padded.window', 'expected': True, 'got': got == ref_seq.upper(), 'ok': got == ref_seq.upper(),
                       'note': f'{ref_gz_rel}: the window reads back base for base from the padded chromosome'})

    # ---------------------------------------------------------------------------------------- 13. outputs
    log('writing')

    def tsv(path, rows, cols):
        with open(path, 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(cols)
            for r in rows:
                w.writerow(['' if r.get(c) is None else (fmt(r[c]) if isinstance(r[c], float) else (int(r[c]) if isinstance(r[c], bool) else r[c])) for c in cols])

    # reads.tsv
    read_rows = []
    through_points = {f'through_j{i + 1}p{k}': p for i, j in enumerate(joins) for k, p in ((1, j['p1']), (2, j['p2']))}
    for name, rec in reads.items():
        row = dict(rec)
        row.pop('alleles', None)
        for k, p in through_points.items():
            row[k] = any(a.mapq >= MAPQ and a.through(p, flank) for a in hifi[name])
        if rec['junction']:
            row['side1_start'], row['side1_end'], row['side1_mapq'] = rec['side1']
            row['side2_start'], row['side2_end'], row['side2_mapq'] = rec['side2']
        row['mat_alleles'] = rec['oh_mat'] + rec['own_mat']
        row['pat_alleles'] = rec['oh_pat'] + rec['own_pat']
        row.pop('alleles_all', None)
        read_rows.append(row)
    read_cols = ['name', 'cell', 'instrument', 'strand', 'length', 'mapq', 'start', 'end', 'n_aln', 'primary', 'parent', 'phased_by',
                 'oh_mat', 'oh_pat', 'own_mat', 'own_pat', 'mat_alleles', 'pat_alleles', 'junction', 'jr_parent', 'jr_mat_alleles', 'jr_pat_alleles',
                 'stop', 'through_S', 'through_E',
                 'side1_start', 'side1_end', 'side1_mapq', 'side2_start', 'side2_end', 'side2_mapq'] + list(through_points) + \
                ['asm_seq', 'asm_piece', 'asm_start', 'asm_end', 'asm_mapq', 'asm_split_start', 'asm_split_end', 'asm_through_M']
    tsv(os.path.join(out, 'reads.tsv'), read_rows, read_cols)
    tsv(os.path.join(out, 'ont_reads.tsv'), ont_rows, ['name', 'length', 'mapq', 'start', 'end', 'junction', 'through_S', 'through_E', 'reach'])
    # sites.tsv
    site_rows = []
    for s in sites:
        row = {'pos': s['pos'], 'ref': s['ref'], 'alt': s['alt'], 'kind': s['kind'], 'qual': s['qual'],
               'child_gt': s['child_gt'], 'child_ref': s['child_ad'][0], 'child_alt': s['child_ad'][1],
               'father_gt': s['father_gt'], 'father_ref': s['father_ad'][0], 'father_alt': s['father_ad'][1],
               'mother_gt': s['mother_gt'], 'mother_ref': s['mother_ad'][0], 'mother_alt': s['mother_ad'][1],
               'mother_allele': s['mother_allele'], 'father_allele': s['father_allele'],
               'hifi_ref': s['hifi_ref'], 'hifi_alt': s['hifi_alt'], 'hifi_dp': s['hifi_dp'],
               'share_sr': s['share_sr'], 'share_hifi': s['share_hifi'], 'in_segment': s['in_segment'], 'in_stretch': s['in_stretch'],
               'copies': s.get('copies'), 'block': s.get('block'), 'ref_on_hap': s.get('ref_on_hap'), 'own_on_hap': s.get('own_on_hap'),
               'hap1_copies': s.get('hap1_copies')}
        for letter, r in pieces.items():
            if r.rs <= s['pos'] <= r.re:
                row[f'base_{letter}'] = '-' if s['pos'] in r.variants_del else (r.variants_snv[s['pos']][1] if s['pos'] in r.variants_snv else s['ref'])
        site_rows.append(row)
    site_cols = ['pos', 'ref', 'alt', 'kind', 'qual', 'child_gt', 'child_ref', 'child_alt', 'father_gt', 'father_ref', 'father_alt',
                 'mother_gt', 'mother_ref', 'mother_alt', 'mother_allele', 'father_allele', 'hifi_ref', 'hifi_alt', 'hifi_dp', 'share_sr',
                 'share_hifi', 'in_segment', 'in_stretch', 'copies', 'block', 'ref_on_hap', 'own_on_hap', 'hap1_copies'] + [f'base_{l}' for l in pieces]
    tsv(os.path.join(out, 'sites.tsv'), site_rows, site_cols)
    tsv(os.path.join(out, 'mother_phase.tsv'), [{'block': i, 'start': b['start'], 'end': b['end'], 'n_sites': b['n']} for i, b in enumerate(blocks)],
        ['block', 'start', 'end', 'n_sites'])
    tsv(os.path.join(out, 'depth_grch38.tsv'), depth_rows, ['bin0', 'bin1', 'sr_child', 'sr_father', 'sr_mother', 'cp_child', 'cp_father', 'cp_mother',
                                                           'cp_rel', 'hifi', 'cp_hifi', 'cp_hifi_mat', 'cp_hifi_pat', 'asm_copies', 'inside', 'flank'])
    tsv(os.path.join(out, 'depth_assembly.tsv'), asm_depth_rows, ['seq', 'bin0', 'bin1', 'mat', 'pat', 'unphased', 'q0', 'piece'])
    piece_rows = [{'piece': r.letter, 'hap': r.local_name, 'contig': r.contig, 'parent_hap': r.parent_hap, 'qs': r.qs, 'qe': r.qe, 'rs': r.rs, 're': r.re,
                   'local_start': r.qs - r.local_offset, 'local_end': r.qe - r.local_offset, 'mapq': r.mapq, 'kind': r.kind, 'nm': r.nm,
                   'aligned': r.aligned, 'mom_own': r.mom_own, 'dad_own': r.dad_own, 'parent': r.parent} for r in pieces.values()]
    tsv(os.path.join(out, 'pieces.tsv'), piece_rows, list(piece_rows[0]))
    tsv(os.path.join(out, 'joins.tsv'), [dict(j, gap_start=(j['gap'] or [None, None])[0], gap_end=(j['gap'] or [None, None])[1],
                                              local_end=rjoins[i]['lend'], local_start=rjoins[i]['rbeg'], cross=rjoins[i]['cross'],
                                              edge_left=rjoins[i].get('edge_left'), edge_right=rjoins[i].get('edge_right')) for i, j in enumerate(joins)],
        ['left', 'right', 'p1', 'p2', 'gap_start', 'gap_end', 'unaligned', 'back', 'carry', 'through1', 'through2', 'stop', 'local_end', 'local_start',
         'cross', 'edge_left', 'edge_right'])
    tsv(os.path.join(out, 'piece_pairs.tsv'), [dict(p, where=json.dumps(p['where']), parents='/'.join(p['parents'])) for p in piece_pairs],
        ['x', 'y', 'parents', 'lo', 'hi', 'len', 'snv', 'indel', 'per_kb', 'where', 'first', 'last', 'same_scaffold'])
    tsv(os.path.join(out, 'qc_placed.tsv'), flags_placed, ['track', 'contig', 'qs', 'qe', 'label', 'piece', 'rs', 're'])
    tsv(os.path.join(out, 'calls.tsv'), calls, ['set', 'made_from', 'type', 'alt', 'pos', 'end', 'len', 'gt', 'id', 'filter', 'qual', 'support', 'reflects'])
    tsv(os.path.join(out, 'sr_pairs.tsv'), sr_pairs_rows, ['member', 'name', 'fwd_start', 'fwd_end', 'fwd_mapq', 'rev_start', 'rev_end', 'rev_mapq'])
    tsv(os.path.join(out, 'self_matches.tsv'), self_rows, ['piece', 'matches', 'longest', 'bp', 'identity'])

    # igv.js tracks (BED, 0-based)
    def bed(path, rows):
        # igv.js splits BED lines on any whitespace, so no field may hold a space
        with open(path, 'w') as f:
            f.write('track itemRgb="On"\n')
            for r in rows:
                f.write('\t'.join(re.sub(r'\s+', '_', str(x)) for x in r) + '\n')

    col = {'paternal': '47,111,219', 'maternal': '184,134,43', 'unknown': '143,150,163'}
    bed(os.path.join(out, 'tracks', 'grch38_segment.bed'), [[chrom, S - 1, E, 'the segment (DUP, 1 copy more on the maternal chromosome)', 0, '+', S - 1, E, '209,51,110']])
    bed(os.path.join(out, 'tracks', 'grch38_pieces.bed'), [[chrom, r.rs - 1, r.re, f'{r.letter} {r.contig}:{r.qs}-{r.qe} ({r.parent} by alleles)', 0, '+', r.rs - 1, r.re, col[r.parent]]
                                                          for r in pieces.values()])
    bed(os.path.join(out, 'tracks', 'grch38_missing.bed'), [[chrom, a - 1, b, f'second maternal copy missing: {b - a + 1} bp', 0, '+', a - 1, b, '209,51,110'] for a, b in holes])
    bed(os.path.join(out, 'tracks', 'grch38_nucflag.bed'), [[chrom, fp['rs'] - 1, fp['re'], f'NucFlag {fp["label"]} on {fp["piece"]}', 0, '+', fp['rs'] - 1, fp['re'], '255,165,0' if fp['label'] == 'MISJOIN' else '200,30,30']
                                                           for fp in flags_placed])
    bed(os.path.join(out, 'tracks', 'grch38_junction_reads.bed'),
        [[chrom, reads[n][side][0] - 1, reads[n][side][1], f'{n} {"side 1" if side == "side1" else "side 2"} ({reads[n]["parent"]})', 0, reads[n]['strand'],
          reads[n][side][0] - 1, reads[n][side][1], '184,134,43' if reads[n]['parent'] == 'mother' else '143,150,163'] for n in jr for side in ('side1', 'side2')])
    kind_col = {'OH': '209,51,110', 'MH': '184,134,43', 'FH': '47,111,219', 'BH': '143,150,163'}
    kind_name = {'OH': 'parent-opposite', 'MH': 'mother-only', 'FH': 'father-only', 'BH': 'both heterozygous'}
    bed(os.path.join(out, 'tracks', 'grch38_sites.bed'), [[chrom, s['pos'] - 1, s['pos'], f'{kind_name[s["kind"]]} {s["ref"]}/{s["alt"]}', 0, '+', s['pos'] - 1, s['pos'], kind_col[s['kind']]] for s in sites])
    bed(os.path.join(out, 'tracks', 'grch38_calls.bed'), [[chrom, c['pos'] - 1, max(c['end'], c['pos']), f'{c["set"]} {c["type"]} {c["len"]:,} bp {c["gt"]} ({c["reflects"]})', 0, '+',
                                                          c['pos'] - 1, max(c['end'], c['pos']), '209,51,110' if c['reflects'] == 'the junction' else ('47,111,219' if c['reflects'] == 'the scaffold' else '143,150,163')]
                                                         for c in calls])
    bed(os.path.join(out, 'tracks', 'asm_pieces.bed'), [[sname, a - 1, b, f'{letter} ({pieces[letter].parent} by alleles; {pieces[letter].contig}:{pieces[letter].qs}-{pieces[letter].qe} = GRCh38 {pieces[letter].rs}-{pieces[letter].re})', 0, '+', a - 1, b, col[pieces[letter].parent]]
                                                       for letter, (sname, a, b) in piece_local.items()])
    bed(os.path.join(out, 'tracks', 'asm_gaps.bed'), [[seq, a - 1, b, f'gap {i + 1}: {b - a + 1} N', 0, '+', a - 1, b, '0,0,0'] for i, (seq, a, b) in enumerate(gaps_local)])
    bed(os.path.join(out, 'tracks', 'asm_m_segment.bed'), [[m_local, ms - 1, me, f'M: the segment once ({M.contig}:{ms + M.local_offset}-{me + M.local_offset})', 0, '+', ms - 1, me, '184,134,43']])
    qc_local = []
    for track, contig, b0, b1, label in qc_rows:
        for name, info in asm['sequences'].items():
            if info['contig'] == contig:
                a, b = b0 + 1 - info['local_offset'], b1 - info['local_offset']
                if b >= 1 and a <= local_len.get(name, 0):
                    qc_local.append([name, max(a, 1) - 1, min(b, local_len[name]), f'{track} {label}', 0, '+', max(a, 1) - 1, min(b, local_len[name]),
                                     '255,165,0' if label == 'MISJOIN' else ('200,30,30' if label == 'COLLAPSE' else ('0,0,0' if label == 'NNN' else '120,120,220'))])
    bed(os.path.join(out, 'tracks', 'asm_nucflag.bed'), [r for r in qc_local if r[3].startswith('nucflag')])
    bed(os.path.join(out, 'tracks', 'asm_flagger.bed'), [r for r in qc_local if r[3].startswith('flagger')])

    # checks
    tsv(os.path.join(out, 'checks.tsv'), [dict(c, expected=json.dumps(c['expected']), got=json.dumps(c['got']), ok=int(c['ok'])) for c in checks],
        ['key', 'expected', 'got', 'ok', 'note'])

    # the bundle for the page
    def table(rows, cols):
        """A compact table: {'cols': [...], 'rows': [[...], ...]}; booleans as 0/1, floats to 4 decimals."""
        def cell(v):
            if v is None or v == '':
                return None
            if isinstance(v, bool):
                return 1 if v else 0
            if isinstance(v, float):
                return round(v, 4)
            return v
        return {'cols': cols, 'rows': [[cell(r.get(c)) for c in cols] for r in rows]}

    bundle_read_cols = [c for c in read_cols if c not in ('n_aln', 'primary', 'mat_alleles', 'pat_alleles')]
    bundle = {
        'case': cfg, 'numbers': numbers,
        'stretch': stretch, 'bounds': bounds, 'patchy_lo': patchy_lo, 'last_oh': last_oh,
        'reads': table([r for r in read_rows if r['primary']], bundle_read_cols),
        'ont_reads': table([r for r in ont_rows if r['reach'] or r['through_S'] or r['through_E'] or r['junction']],
                           ['name', 'length', 'mapq', 'start', 'end', 'junction', 'through_S', 'through_E', 'reach']),
        'sites': table(site_rows, site_cols),
        'blocks': [{'id': i, 'start': b['start'], 'end': b['end'], 'n': b['n']} for i, b in enumerate(blocks)],
        'depth_grch38': table(depth_rows, ['bin0', 'bin1', 'sr_child', 'sr_father', 'sr_mother', 'cp_child', 'cp_father', 'cp_mother', 'cp_rel',
                                           'hifi', 'cp_hifi', 'cp_hifi_mat', 'cp_hifi_pat', 'phased_reads_mat', 'phased_reads_pat', 'asm_copies', 'inside', 'flank']),
        'depth_assembly': table(asm_depth_rows, ['seq', 'bin0', 'bin1', 'mat', 'pat', 'unphased', 'q0', 'piece']),
        'pieces': piece_rows, 'piece_local': piece_local, 'local_len': local_len, 'm_segment': m_segment, 'gaps_local': gaps_local,
        'joins': [dict(j, realign=rjoins[i]) for i, j in enumerate(joins)], 'holes': holes, 'piece_pairs': piece_pairs,
        'missing': missing, 'qc_placed': flags_placed, 'qc_rows': qc_rows, 'asm_vs_asm': asm_vs_asm,
        'calls': calls, 'yale_summary': yale_summary, 'published': published, 'carriers': sorted(set(carriers)),
        'sr': {m: {k: v for k, v in e.items() if k != 'split_names'} for m, e in sr_evidence.items()}, 'sr_pairs': sr_pairs_rows,
        'sr_flank': sr_flank, 'hifi_flank': hifi_flank, 'one_copy': one_copy, 'flank_depth': flank_depth,
        'self': {'pieces': self_pieces, 'first': self_first, 'other': self_other, 'bp': paf_bp, 'rows': self_rows}, 'mh': mh,
        'ont': {'junction': ont_j, 'through_s': ont_ts, 'through_e': ont_te, 'reads': len(ont_rows), 'reach': ont_reach, 'conflicts': conflicts, 'linked': len(linked_reads)},
        'compo': {k: dict(v) for k, v in compo.items()}, 'switch': {k: dict(v) for k, v in switch.items()},
        'switch_reads': {k: v for k, v in switch_reads.items()},
        'read_alleles': {name: sorted((rec['alleles_all'] if rec['junction'] else rec['alleles']).items()) for name, rec in reads.items()
                         if rec['junction'] or any(name == n for end in switch_reads.values() for n, _ in end)},
        'checks': checks,
    }
    with open(os.path.join(out, 'bundle.js'), 'w') as f:
        f.write('window.ASSEMBLYREVIEW_CASE = ')
        json.dump(bundle, f, separators=(',', ':'))
        f.write(';\n')
    # the sequences for the dot plots
    def ref_window(lo, hi):
        return {'name': f'GRCh38 {chrom}:{lo:,}-{hi:,}', 'coord': [chrom, lo], 'seq': ref_seq[lo - ref_offset:hi - ref_offset + 1]}

    def local_window(seq, lo, hi, label):
        lo, hi = max(1, lo), min(local_len[seq], hi)
        return {'name': f'{label} {seq}:{lo:,}-{hi:,}', 'coord': [seq, lo], 'seq': local[seq][lo - 1:hi]}

    g1 = gaps_local[0] if gaps_local else None
    g4 = gaps_local[-1] if gaps_local else None
    windows = {
        'grch38_start': ref_window(S - 5000, S + 30000), 'grch38_end': ref_window(E - 30000, E + 5000),
        'M_start': local_window(m_local, ms - 5000, ms + 30000, 'M'), 'M_end': local_window(m_local, me - 30000, me + 5000, 'M'),
    }
    if g1:
        windows['hap1_gap1'] = local_window(g1[0], g1[1] - 15000, g1[2] + 15000, 'paternal scaffold, gap 1')
    if g4:
        windows['hap1_gap4'] = local_window(g4[0], g4[1] - 15000, g4[2] + 15000, 'paternal scaffold, gap 4')
    seqs = {'windows': windows, 'reads': {n: seqs_reads[n] for n in seqs_reads}, 'segment': [S, E], 'm_segment': m_segment,
            'gaps_local': gaps_local, 'junction_reads': jr}
    with open(os.path.join(out, 'sequences.js'), 'w') as f:
        f.write('window.ASSEMBLYREVIEW_SEQ = ')
        json.dump(seqs, f, separators=(',', ':'))
        f.write(';\n')

    bad = [c for c in checks if not c['ok']]
    log(f'{len(checks)} checks, {len(checks) - len(bad)} match numbers.json')
    for c in bad:
        log(f'  MISMATCH {c["key"]}: expected {c["expected"]!r}, got {c["got"]!r} {c["note"]}')
    return 0 if not bad or args.no_check else 1


if __name__ == '__main__':
    sys.exit(main())
