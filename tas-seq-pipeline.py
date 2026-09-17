#!/usr/bin/env python3
"""
tas-seq-pipeline : Integrated TAS-Seq pipeline (single-file, self-contained)

Everything is in this one file:
  - bd-extract core (BD Rhapsody Enhanced-beads CB+UMI extraction & correction)
  - sam demultiplexing  (former sam_demutiplexing.rb)
  - parallel SAM->BAM   (former sam2bam.rb)
  - parallel salmon quant (former salmon_HT.rb, salmon -p 2)

External tools required on PATH: fastp, bwa, samtools, salmon
  Python deps: pysam (for dedup). umi_tools no longer required.

Flow:
  extract (CLS-list correction) -> fastp -> bwa mem -> sorted BAM
    -> directional dedup (self) -> demultiplex SAM -> SAM->BAM -> salmon quant
    -> summary TSV (UMI conversion rate etc.)
    -> cleanup (keep: salmon/, fastp, reads, summary)

CPU threads auto-detected (nproc). salmon runs 2 threads per job.

Usage:
  tas-seq-pipeline \
    -i sample_R1.fastq.gz -f reference.fa \
    --cls1 BD1_CB.txt --cls2 BD2_CB.txt --cls3 BD3_CB.txt \
    [--adapter TAS-adapter.fa] [--threads N] [--keep-intermediates]
"""

import argparse
import json
import os
import sys
import glob
import gzip
import shutil
import time
import itertools
import subprocess
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from collections import defaultdict
# pysam is imported lazily inside dedup functions (only needed for Step 4)

# ==============================================================================
# PART 1: bd-extract core (embedded)
# ==============================================================================
# ── Linker sequences (Enhanced beads) ────────────────────────────────────────
L1_SEQ = "GTGA"
L2_SEQ = "GACA"
CLS_LEN = 9
UMI_LEN = 8

# ── CB extraction ─────────────────────────────────────────────────────────────
def extract_enhanced(r1seq: str):
    for vb_len in (0, 1, 2, 3):
        l1_pos = vb_len + CLS_LEN
        if r1seq[l1_pos:l1_pos+4] == L1_SEQ:
            cls2_pos = l1_pos + 4
            l2_pos   = cls2_pos + CLS_LEN
            if r1seq[l2_pos:l2_pos+4] != L2_SEQ:
                continue
            cls3_pos = l2_pos + 4
            umi_pos  = cls3_pos + CLS_LEN
            return (r1seq[vb_len:l1_pos]
                    + r1seq[cls2_pos:cls2_pos+CLS_LEN]
                    + r1seq[cls3_pos:cls3_pos+CLS_LEN],
                    r1seq[umi_pos:umi_pos+UMI_LEN])
    return None, None

def extract_enhanced_1mm(r1seq: str):
    cb, umi = extract_enhanced(r1seq)
    if cb is not None:
        return cb, umi
    for l1_pos in range(9, 13):
        obs_l1 = r1seq[l1_pos:l1_pos+4]
        if sum(a != b for a, b in zip(obs_l1, L1_SEQ)) <= 1:
            cls2_pos = l1_pos + 4
            l2_pos   = cls2_pos + CLS_LEN
            if sum(a != b for a, b in zip(r1seq[l2_pos:l2_pos+4], L2_SEQ)) <= 1:
                cls3_pos = l2_pos + 4
                umi_pos  = cls3_pos + CLS_LEN
                return (r1seq[l1_pos-CLS_LEN:l1_pos]
                        + r1seq[cls2_pos:cls2_pos+CLS_LEN]
                        + r1seq[cls3_pos:cls3_pos+CLS_LEN],
                        r1seq[umi_pos:umi_pos+UMI_LEN])
    return None, None

def extract_v1(r1seq: str):
    return (r1seq[0:9] + r1seq[21:30] + r1seq[43:52],
            r1seq[52:60])

# ── Per-part CLS correction (BD official CLS lists) ───────────────────────────
def load_cls_list(path: str):
    """Load a CLS sequence list file (one sequence per line, ignores blanks/headers)."""
    seqs = []
    with open(path) as f:
        for line in f:
            s = line.strip().split()[0] if line.strip() else ""
            # accept lines that look like DNA only
            if s and all(c in "ACGTN" for c in s.upper()):
                seqs.append(s.upper())
    return seqs

def build_part_correction_map(cls_seqs, mismatch: int):
    """
    Build {observed_part -> correct_part} for ONE CLS position.
    Ambiguous 1MM (near 2+ valid parts) -> discarded.
    """
    cmap = {}
    ambiguous = set()
    for s in cls_seqs:
        cmap[s] = s  # exact
    if mismatch >= 1:
        for s in cls_seqs:
            for variant in _build_1mm_variants(s):
                if variant in cmap:
                    if cmap[variant] != s:
                        ambiguous.add(variant)
                else:
                    cmap[variant] = s
        for v in ambiguous:
            cmap.pop(v, None)
    return cmap, len(ambiguous)

class PartCorrector:
    """
    Holds three per-part correction maps (CLS1, CLS2, CLS3).
    correct(cb27) -> corrected 27nt CB or None (if any part fails/ambiguous).
    """
    def __init__(self, cls1_path, cls2_path, cls3_path, mismatch):
        self.cls1 = load_cls_list(cls1_path)
        self.cls2 = load_cls_list(cls2_path)
        self.cls3 = load_cls_list(cls3_path)
        self.m1, a1 = build_part_correction_map(self.cls1, mismatch)
        self.m2, a2 = build_part_correction_map(self.cls2, mismatch)
        self.m3, a3 = build_part_correction_map(self.cls3, mismatch)
        print(f"[bd-extract] CLS lists: {len(self.cls1)}/{len(self.cls2)}/"
              f"{len(self.cls3)} seqs", file=sys.stderr)
        print(f"[bd-extract] Per-part 1MM maps: "
              f"{len(self.m1)}/{len(self.m2)}/{len(self.m3)} entries "
              f"(ambiguous removed: {a1}/{a2}/{a3})", file=sys.stderr)
        combos = len(self.cls1) * len(self.cls2) * len(self.cls3)
        print(f"[bd-extract] Theoretical CB space: {combos:,}", file=sys.stderr)

    def correct(self, cb27: str):
        p1 = self.m1.get(cb27[0:9])
        if p1 is None:
            return None
        p2 = self.m2.get(cb27[9:18])
        if p2 is None:
            return None
        p3 = self.m3.get(cb27[18:27])
        if p3 is None:
            return None
        return p1 + p2 + p3

# ── Correction map builders ───────────────────────────────────────────────────
def _build_1mm_variants(cb: str):
    """Generate all 1MM variants of cb."""
    for pos in range(len(cb)):
        orig = cb[pos]
        for base in "ACGTN":
            if base != orig:
                yield cb[:pos] + base + cb[pos+1:]

def build_correction_map_from_whitelist(whitelist_path: str, mismatch: int):
    """
    From a whitelist TSV (umi-tools format):
      col1: correct CB
      col2: comma-separated error CBs (optional)
    Returns correction_map: {observed_cb -> correct_cb}
    Ambiguous (1MM to multiple correct CBs) -> not added (will be discarded)
    """
    correct_cbs = []
    correction_map = {}   # observed -> correct
    ambiguous = set()

    with open(whitelist_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            cb = parts[0]
            correct_cbs.append(cb)
            correction_map[cb] = cb   # exact
            if len(parts) >= 2 and parts[1]:
                for err in parts[1].split(","):
                    err = err.strip()
                    if err:
                        if err in correction_map and correction_map[err] != cb:
                            ambiguous.add(err)
                        else:
                            correction_map[err] = cb

    if mismatch >= 1:
        print(f"[bd-extract] Building 1MM map ({len(correct_cbs)} CBs)...",
              file=sys.stderr)
        for cb in correct_cbs:
            for variant in _build_1mm_variants(cb):
                if variant in correction_map:
                    if correction_map[variant] != cb:
                        ambiguous.add(variant)  # ambiguous
                else:
                    correction_map[variant] = cb
        # Remove ambiguous
        for v in ambiguous:
            correction_map.pop(v, None)
        print(f"[bd-extract] 1MM map: {len(correction_map)} entries "
              f"({len(ambiguous)} ambiguous removed)", file=sys.stderr)

    return correction_map

def build_correction_map_from_counts(cb_counts: dict, knee_threshold: int, mismatch: int):
    """
    From CB frequency dict, detect valid CBs by knee threshold,
    build correction map toward valid CB set.
    Ambiguous -> discarded.
    """
    valid_cbs = [cb for cb, n in cb_counts.items() if n >= knee_threshold]
    print(f"[bd-extract] Knee threshold: {knee_threshold}  "
          f"Valid CBs: {len(valid_cbs)}", file=sys.stderr)

    correction_map = {}
    ambiguous = set()

    for cb in valid_cbs:
        correction_map[cb] = cb   # exact

    if mismatch >= 1:
        print(f"[bd-extract] Building 1MM map ({len(valid_cbs)} valid CBs)...",
              file=sys.stderr)
        for cb in valid_cbs:
            for variant in _build_1mm_variants(cb):
                if variant in correction_map:
                    if correction_map[variant] != cb:
                        ambiguous.add(variant)
                else:
                    correction_map[variant] = cb
        for v in ambiguous:
            correction_map.pop(v, None)
        print(f"[bd-extract] 1MM map: {len(correction_map)} entries "
              f"({len(ambiguous)} ambiguous removed)", file=sys.stderr)

    return correction_map

# ── Knee detection ────────────────────────────────────────────────────────────
def detect_knee(cb_counts: dict, expected_cells: int = None):
    """
    Simple knee detection by sorted CB count drop.
    Returns threshold count.
    """
    counts = sorted(cb_counts.values(), reverse=True)
    if not counts:
        return 1

    if expected_cells and expected_cells < len(counts):
        # Use expected cell count as a guide (pick inflection near it)
        window = max(1, expected_cells // 10)
        lo = max(0, expected_cells - window)
        hi = min(len(counts)-1, expected_cells + window)
        # Find steepest drop in window
        best_drop = 0
        best_idx = expected_cells
        for i in range(lo, hi):
            drop = counts[i] - counts[i+1] if i+1 < len(counts) else 0
            if drop > best_drop:
                best_drop = drop
                best_idx = i
        return counts[best_idx]

    # Fallback: find global steepest drop in log space
    import math
    log_counts = [math.log10(c+1) for c in counts]
    best_drop = 0
    best_idx = 0
    for i in range(len(log_counts)-1):
        drop = log_counts[i] - log_counts[i+1]
        if drop > best_drop:
            best_drop = drop
            best_idx = i
    return counts[best_idx]

# ── FASTQ I/O ─────────────────────────────────────────────────────────────────
# Prefer external multi-threaded (de)compressors to keep the main process from
# becoming a serial gzip bottleneck. Fallback chain: pigz -> igzip(isal) -> gzip.
import shutil as _shutil

def _find_decompressor():
    for tool in ("pigz", "unpigz", "igzip"):
        if _shutil.which(tool):
            return [tool, "-dc"]
    # isal python module CLI (no separate binary needed)
    try:
        import isal.igzip  # noqa
        return [sys.executable, "-m", "isal.igzip", "-dc"]
    except ImportError:
        return None

def _find_compressor():
    if _shutil.which("pigz"):
        return ["pigz", "-1", "-c"]
    if _shutil.which("igzip"):
        return ["igzip", "-1", "-c"]
    try:
        import isal.igzip  # noqa
        return [sys.executable, "-m", "isal.igzip", "-1", "-c"]
    except ImportError:
        return None

_DECOMP = _find_decompressor()
_COMP = _find_compressor()

def open_fastq_read(path):
    """Return (file_obj, proc_or_None). Streams via pigz/igzip if available."""
    if path.endswith(".gz") and _DECOMP:
        cmd = list(_DECOMP) + [path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                                bufsize=1 << 20)
        return proc.stdout, proc
    if path.endswith(".gz"):
        try:
            from isal import igzip as _igzip
            return _igzip.open(path, "rt"), None
        except ImportError:
            return gzip.open(path, "rt", encoding="ascii"), None
    return open(path, "r"), None

def open_fastq_write(path, threads=4):
    """Return (file_obj, proc_or_None) for writing, streaming via pigz/igzip."""
    if path.endswith(".gz") and _COMP:
        fout = open(path, "wb")
        cmd = list(_COMP)
        # pigz supports -p for thread count
        if cmd[0] == "pigz":
            cmd = ["pigz", "-1", "-p", str(threads), "-c"]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=fout,
                                text=True, bufsize=1 << 20)
        return proc.stdin, (proc, fout)
    if path.endswith(".gz"):
        try:
            from isal import igzip as _igzip
            return _igzip.open(path, "wt", compresslevel=1), None
        except ImportError:
            return gzip.open(path, "wt", compresslevel=1), None
    return open(path, "w"), None

def close_write(handle, meta):
    handle.close()
    if meta is not None:
        proc, fout = meta
        proc.wait()
        fout.close()

def open_fastq(path):
    """Legacy single-handle reader (used by knee-count pass)."""
    f, _ = open_fastq_read(path)
    return f

def read_fastq_chunks(path, chunk_size):
    f, proc = open_fastq_read(path)
    chunk = []
    try:
        while True:
            name = f.readline()
            if not name:
                break
            seq  = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")
            qual = f.readline().rstrip("\n")
            chunk.append((name.rstrip("\n"), seq, plus, qual))
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk
    finally:
        f.close()
        if proc is not None:
            proc.wait()

def _read_fastq_chunks_orig(path, chunk_size):
    chunk = []
    with open_fastq(path) as f:
        while True:
            name = f.readline()
            if not name:
                break
            seq  = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")
            qual = f.readline().rstrip("\n")
            chunk.append((name.rstrip("\n"), seq, plus, qual))
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    if chunk:
        yield chunk

# ── Pass 1: count CB frequencies ─────────────────────────────────────────────
def count_cbs_chunk(args):
    r1_chunk, bead, linker_mm = args
    if bead == "enhanced":
        extractor = extract_enhanced_1mm if linker_mm else extract_enhanced
    else:
        extractor = extract_v1
    counts = defaultdict(int)
    n_nolinker = 0
    for (_, r1seq, _, _) in r1_chunk:
        cb, _ = extractor(r1seq)
        if cb is None:
            n_nolinker += 1
        else:
            counts[cb] += 1
    return dict(counts), n_nolinker

# ── Pass 2 / extract worker ───────────────────────────────────────────────────
def process_chunk(args):
    r1_chunk, r2_chunk, bead, correction_map, linker_mm, want_r1 = args

    if bead == "enhanced":
        extractor = extract_enhanced_1mm if linker_mm else extract_enhanced
    else:
        extractor = extract_v1

    # correction_map can be a dict {obs->correct} or a PartCorrector
    is_part = isinstance(correction_map, PartCorrector)

    r1_out = []
    r2_out = []
    n_total = n_pass = n_corrected = n_fail = n_nolinker = 0

    for (r1name, r1seq, r1plus, r1qual), (r2name, r2seq, r2plus, r2qual) in \
            zip(r1_chunk, r2_chunk):
        n_total += 1

        cb, umi = extractor(r1seq)
        if cb is None:
            n_nolinker += 1
            n_fail += 1
            continue

        if is_part:
            corrected = correction_map.correct(cb)
        else:
            corrected = correction_map.get(cb)
        if corrected is None:
            n_fail += 1
            continue
        if corrected != cb:
            n_corrected += 1
        cb = corrected

        n_pass += 1
        read_id = r1name.split()[0]
        tag = f"{read_id}_{cb}_{umi}"
        if want_r1:
            r1_out.append(f"{tag}\n{r1seq}\n{r1plus}\n{r1qual}\n")
        r2_out.append(f"{tag}\n{r2seq}\n{r2plus}\n{r2qual}\n")

    stats = dict(total=n_total, pass_=n_pass, corrected=n_corrected,
                 fail=n_fail, nolinker=n_nolinker)
    return ("".join(r1_out) if want_r1 else ""), "".join(r2_out), stats

# ── Writer ────────────────────────────────────────────────────────────────────
def write_output(out_path, text):
    if out_path.endswith(".gz"):
        with gzip.open(out_path, "at", compresslevel=1) as f:
            f.write(text)
    else:
        with open(out_path, "a") as f:
            f.write(text)

# ── Parallel runner ───────────────────────────────────────────────────────────
def run_extract(r1_path, r2_path, r1_out, r2_out,
                correction_map, bead, linker_mm,
                threads, chunk_size):
    """
    r1_out=None  -> skip writing R1 entirely (faster; use when only R2 is needed).
    Output streams through pigz/igzip when available so the main process is not
    a serial gzip bottleneck.
    """
    want_r1 = r1_out is not None
    for p in [r1_out, r2_out]:
        if p and os.path.exists(p):
            os.remove(p)

    t0 = time.time()
    total_stats = defaultdict(int)
    read_count = 0

    # Compression threads: give each writer a slice, leave the rest for workers
    comp_threads = max(2, threads // 4)

    # Open streaming writers (stay open for the whole run; no per-chunk reopen)
    w2, w2_meta = open_fastq_write(r2_out, comp_threads)
    w1, w1_meta = (open_fastq_write(r1_out, comp_threads)
                   if want_r1 else (None, None))

    r1_iter = read_fastq_chunks(r1_path, chunk_size)
    r2_iter = read_fastq_chunks(r2_path, chunk_size)
    chunk_args = (
        (r1c, r2c, bead, correction_map, linker_mm, want_r1)
        for r1c, r2c in zip(r1_iter, r2_iter)
    )

    try:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futures = []
            chunk_iter = iter(chunk_args)
            for ca in itertools.islice(chunk_iter, threads * 2):
                futures.append(pool.submit(process_chunk, ca))

            while futures:
                r1t, r2t, stats = futures.pop(0).result()
                if want_r1 and r1t:
                    w1.write(r1t)
                if r2t:
                    w2.write(r2t)
                for k, v in stats.items():
                    total_stats[k] += v
                read_count += stats["total"]
                elapsed = time.time() - t0
                speed = read_count / max(elapsed, 0.001) / 1e6
                print(f"\r[bd-extract] {read_count:,} reads | {speed:.2f}M/s | "
                      f"pass={total_stats['pass_']:,} "
                      f"nolinker={total_stats['nolinker']:,} "
                      f"fail={total_stats['fail']:,}",
                      end="", file=sys.stderr)
                try:
                    futures.append(pool.submit(process_chunk, next(chunk_iter)))
                except StopIteration:
                    pass
    finally:
        if want_r1:
            close_write(w1, w1_meta)
        close_write(w2, w2_meta)

    elapsed = time.time() - t0
    t = total_stats
    pct = lambda n: f"{100*n/max(t['total'],1):.1f}%"
    print(f"\n[bd-extract] Done in {elapsed:.1f}s", file=sys.stderr)
    print(f"  Total:        {t['total']:>10,}", file=sys.stderr)
    print(f"  Pass:         {t['pass_']:>10,}  ({pct(t['pass_'])})", file=sys.stderr)
    print(f"  CB corrected: {t['corrected']:>10,}", file=sys.stderr)
    print(f"  No linker:    {t['nolinker']:>10,}  ({pct(t['nolinker'])})", file=sys.stderr)
    print(f"  Discarded:    {t['fail']:>10,}  ({pct(t['fail'])})", file=sys.stderr)
    return dict(t)


# ==============================================================================
# PART 2: pipeline orchestration
# ==============================================================================
def plog(msg):
    print(f"[tas-seq] {msg}", file=sys.stderr, flush=True)

def run_cmd(cmd, **kwargs):
    plog(f"$ {cmd}")
    result = subprocess.run(cmd, shell=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {cmd}")
    return result

# ── Resume support: per-step marker files (workdir/.step_<name>.done) ────────
# Each marker holds a small JSON blob with whatever data the step needs to
# hand off to the summary section if it gets skipped on a later --resume run.
def _step_marker(workdir, name):
    return os.path.join(workdir, f".step_{name}.done")

def mark_step_done(workdir, name, data=None):
    with open(_step_marker(workdir, name), "w") as f:
        json.dump(data or {}, f)

def step_done_data(workdir, name, required_files=()):
    """Return the marker's JSON data if the step is done and all
    required_files still exist, else None."""
    path = _step_marker(workdir, name)
    if not os.path.exists(path):
        return None
    for rf in required_files:
        if not os.path.exists(rf):
            return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return {}

# ── Directional UMI dedup (self-implemented, replaces umi_tools dedup) ─────────
# Smith, Heger & Sudbery (2017). Per (contig, cell):
#   edge A->B if hamming(A,B)==1 and count_A >= 2*count_B - 1
#   connected components from high-count hubs = unique molecules.
# Equivalent to: umi_tools dedup --per-contig --per-gene --per-cell
# Representative read per molecule chosen deterministically.
def _hamming1(a, b):
    if len(a) != len(b):
        return False
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > 1:
                return False
    return d == 1

def _directional_clusters(umi_counts):
    """umi_counts: {umi->count} -> list of (hub_umi, {member umis})."""
    umis = sorted(umi_counts.keys(), key=lambda u: (-umi_counts[u], u))
    adj = defaultdict(list)
    n = len(umis)
    for i in range(n):
        ui = umis[i]; ci = umi_counts[ui]
        for j in range(n):
            if i == j:
                continue
            uj = umis[j]; cj = umi_counts[uj]
            if ci >= 2 * cj - 1 and _hamming1(ui, uj):
                adj[ui].append(uj)
    visited = set()
    clusters = []
    for u in umis:                      # high-count hubs first
        if u in visited:
            continue
        stack = [u]; comp = set()
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            for nb in adj[node]:
                if nb not in visited:
                    stack.append(nb)
        clusters.append((u, comp))       # u = hub (highest count) = representative UMI
    return clusters

def _dedup_one_group(args):
    """
    Dedup the contigs in one group. Reads the source BAM region, keeps one
    representative read per (contig, cell, molecule), writes a sub-BAM.
    Returns sub-BAM path.
    """
    import pysam
    in_bam, contigs, out_bam = args
    src = pysam.AlignmentFile(in_bam)

    # Gather reads per (contig, cell): {umi -> list of (sort_key, read)}
    # sort_key makes representative selection deterministic.
    from collections import defaultdict as dd
    groups = dd(lambda: dd(list))   # (contig,cell) -> umi -> [(key, aln)]
    for c in contigs:
        for a in src.fetch(c):
            if a.is_unmapped:
                continue
            parts = a.query_name.split("_")
            if len(parts) < 2:
                continue
            cell, umi = parts[-2], parts[-1]
            key = (a.reference_start, a.query_name)  # deterministic tie-break
            groups[(c, cell)][umi].append((key, a))

    dst = pysam.AlignmentFile(out_bam, "wb", template=src)
    for (contig, cell), umi_map in groups.items():
        umi_counts = {u: len(v) for u, v in umi_map.items()}
        for hub_umi, members in _directional_clusters(umi_counts):
            # collect all reads of the molecule (hub + absorbed satellites)
            mol_reads = []
            for u in members:
                mol_reads.extend(umi_map[u])
            # representative = smallest (reference_start, query_name)
            mol_reads.sort(key=lambda x: x[0])
            dst.write(mol_reads[0][1])
    dst.close()
    src.close()
    return out_bam

def parallel_dedup(in_bam, out_sam, threads, n_groups=None):
    """
    Parallel directional dedup by splitting contigs into balanced groups.
    Writes a single merged SAM (header from source). Result is deterministic.
    """
    import pysam
    src = pysam.AlignmentFile(in_bam)
    contigs = list(src.references)
    # read counts per contig for load balancing
    counts = {}
    for st in src.get_index_statistics():
        counts[st.contig] = st.mapped
    src.close()
    contigs = [c for c in contigs if counts.get(c, 0) > 0]
    if not contigs:
        plog("  WARNING: no mapped contigs found")
        # still produce an (almost empty) SAM with header
        s = pysam.AlignmentFile(in_bam)
        pysam.AlignmentFile(out_sam, "w", template=s).close()
        s.close()
        return (0, 0)

    n_in = sum(counts[c] for c in contigs)   # mapped reads entering dedup

    if n_groups is None:
        n_groups = threads
    n_groups = max(1, min(n_groups, len(contigs)))

    # greedy bin-packing by read count
    groups = [[] for _ in range(n_groups)]
    load = [0] * n_groups
    for c in sorted(contigs, key=lambda x: -counts[x]):
        g = load.index(min(load))
        groups[g].append(c)
        load[g] += counts[c]
    groups = [g for g in groups if g]

    tmpdir = out_sam + ".dd_tmp"
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    plog(f"  directional dedup: {len(contigs):,} contigs -> "
         f"{len(groups)} parallel groups")
    job_args = [(in_bam, groups[i], os.path.join(tmpdir, f"g{i}.bam"))
                for i in range(len(groups))]

    sub_bams = []
    with ProcessPoolExecutor(max_workers=threads) as pool:
        for ob in pool.map(_dedup_one_group, job_args):
            sub_bams.append(ob)

    # merge sub-BAMs into one SAM (header from source)
    src = pysam.AlignmentFile(in_bam)
    n_out = 0
    with pysam.AlignmentFile(out_sam, "w", template=src) as merged:
        for ob in sorted(sub_bams):
            sub = pysam.AlignmentFile(ob)
            for a in sub.fetch(until_eof=True):
                merged.write(a)
                n_out += 1
            sub.close()
    src.close()
    shutil.rmtree(tmpdir)
    conv = n_out / n_in if n_in else 0.0
    plog(f"  dedup done: {n_in:,} reads -> {n_out:,} molecules "
         f"(UMI conversion {100*conv:.1f}%)")
    return (n_in, n_out)

# ── demultiplex SAM by cell barcode (former sam_demutiplexing.rb) ─────────────
def demultiplex_sam(sam_path, out_dir):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    plog(f"Demultiplexing {sam_path} -> {out_dir}/")
    handles = {}
    MAX_OPEN = 500
    n_lines = 0
    try:
        with open(sam_path) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                tab = line.find("\t")
                if tab < 0:
                    continue
                fields = line[:tab].split("_")
                if len(fields) < 2:
                    continue
                bc = fields[1]
                fh = handles.get(bc)
                if fh is None:
                    if len(handles) >= MAX_OPEN:
                        for h in handles.values():
                            h.close()
                        handles.clear()
                    fh = open(os.path.join(out_dir, f"{bc}.sam"), "a")
                    handles[bc] = fh
                fh.write(line)
                n_lines += 1
    finally:
        for h in handles.values():
            h.close()
    n_bc = len(glob.glob(os.path.join(out_dir, "*.sam")))
    plog(f"  {n_lines:,} alignments -> {n_bc:,} barcodes")

# ── SAM -> BAM in parallel (former sam2bam.rb) ────────────────────────────────
def sam_to_bam_one(args):
    sam_path, refseq, out_dir = args
    base = os.path.splitext(os.path.basename(sam_path))[0]
    bam = os.path.join(out_dir, f"{base}.bam")
    rc = subprocess.run(
        f"samtools view -@ 2 -Shb {sam_path} -T {refseq} -o {bam} && "
        f"samtools index {bam}",
        shell=True, stderr=subprocess.DEVNULL).returncode
    return bam if rc == 0 else None

def sam_to_bam_all(sam_dir, refseq, out_dir, threads, resume=False):
    if os.path.exists(out_dir):
        if not resume:
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)
    sam_files = sorted(glob.glob(os.path.join(sam_dir, "*.sam")))
    if resume:
        todo = []
        n_skip = 0
        for s in sam_files:
            base = os.path.splitext(os.path.basename(s))[0]
            bam = os.path.join(out_dir, f"{base}.bam")
            if os.path.exists(bam) and os.path.exists(bam + ".bai") \
                    and os.path.getsize(bam) > 0:
                n_skip += 1
            else:
                todo.append(s)
        if n_skip:
            plog(f"  resume: {n_skip:,}/{len(sam_files):,} BAMs already present, skipping")
        sam_files = todo
    plog(f"SAM->BAM: {len(sam_files):,} files, {threads} threads")
    args = [(s, refseq, out_dir) for s in sam_files]
    done = 0
    with ThreadPoolExecutor(max_workers=threads) as pool:
        for _ in pool.map(sam_to_bam_one, args):
            done += 1
            if done % 500 == 0:
                plog(f"  {done:,}/{len(sam_files):,}")
    bams = sorted(glob.glob(os.path.join(out_dir, "*.bam")))
    plog(f"  -> {len(bams):,} BAMs")
    return bams

# ── salmon quant in parallel (former salmon_HT.rb, salmon -p 2) ───────────────
def salmon_one(args):
    bam_path, fasta, out_root = args
    base = os.path.splitext(os.path.basename(bam_path))[0]
    out_dir = os.path.join(out_root, base)
    rc = subprocess.run(
        f"salmon quant -t {fasta} -l IU -p 2 -a {bam_path} -o {out_dir}",
        shell=True, stderr=subprocess.DEVNULL).returncode
    logs = os.path.join(out_dir, "logs")
    if os.path.isdir(logs):
        shutil.rmtree(logs, ignore_errors=True)
    return base if rc == 0 else None

def salmon_all(bam_dir, fasta, out_root, jobs, resume=False):
    os.makedirs(out_root, exist_ok=True)
    bam_files = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))
    if resume:
        todo = []
        n_skip = 0
        for b in bam_files:
            base = os.path.splitext(os.path.basename(b))[0]
            qsf = os.path.join(out_root, base, "quant.sf")
            if os.path.exists(qsf) and os.path.getsize(qsf) > 0:
                n_skip += 1
            else:
                todo.append(b)
        if n_skip:
            plog(f"  resume: {n_skip:,}/{len(bam_files):,} quant dirs already present, skipping")
        bam_files = todo
    plog(f"salmon quant: {len(bam_files):,} BAMs, {jobs} concurrent jobs x 2 threads")
    args = [(b, fasta, out_root) for b in bam_files]
    done = 0
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for _ in pool.map(salmon_one, args):
            done += 1
            if done % 500 == 0:
                plog(f"  {done:,}/{len(bam_files):,}")
    quants = [d for d in glob.glob(os.path.join(out_root, "*")) if os.path.isdir(d)]
    plog(f"  -> {len(quants):,} quant dirs")

# ── R2 extracted read count per barcode ───────────────────────────────────────
def r2_read_count(trimmed_fastq, out_txt):
    plog(f"Counting R2 reads per barcode -> {out_txt}")
    counts = defaultdict(int)
    opener = gzip.open if trimmed_fastq.endswith(".gz") else open
    with opener(trimmed_fastq, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                parts = line.strip().split("_")
                if len(parts) >= 2:
                    counts[parts[1]] += 1
    with open(out_txt, "w") as out:
        for bc in sorted(counts):
            out.write(f"{bc}\t{counts[bc]}\n")
    plog(f"  {len(counts):,} barcodes")

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        prog="tas-seq-pipeline",
        description="Integrated single-file TAS-Seq pipeline "
                    "(embedded bd-extract + fastp + bwa + directional dedup + "
                    "demux + salmon)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True,
                   help="Input R1 FASTQ(.gz). R2 inferred by replacing R1->R2.")
    p.add_argument("-f", "--ref", required=True,
                   help="Reference FASTA (bwa index prefix & salmon -t target)")
    p.add_argument("--cls1", default="/mnt/hdd2/BD_CLS1_enh.txt",
                   help="CLS1 barcode list (default: /mnt/hdd2/BD_CLS1_enh.txt)")
    p.add_argument("--cls2", default="/mnt/hdd2/BD_CLS2_enh.txt",
                   help="CLS2 barcode list (default: /mnt/hdd2/BD_CLS2_enh.txt)")
    p.add_argument("--cls3", default="/mnt/hdd2/BD_CLS3_enh.txt",
                   help="CLS3 barcode list (default: /mnt/hdd2/BD_CLS3_enh.txt)")
    p.add_argument("--adapter", default="/mnt/hdd2/MK_bin/TAS-adapter.fa",
                   help="fastp adapter FASTA")
    p.add_argument("--mismatch", type=int, default=1,
                   help="CB correction mismatch (default: 1)")
    p.add_argument("--threads", type=int, default=multiprocessing.cpu_count(),
                   help="CPU threads (default: auto = nproc)")
    p.add_argument("--min-len", type=int, default=31,
                   help="fastp minimum length (default: 31)")
    p.add_argument("--chunk-size", type=int, default=100000,
                   help="bd-extract reads per chunk (default: 100000)")
    p.add_argument("--linker-mm", action="store_true",
                   help="Allow 1MM in linker sequences during extraction")
    p.add_argument("--keep-intermediates", action="store_true",
                   help="Do NOT delete intermediate bam/sam/gz at the end")
    p.add_argument("--resume", action="store_true",
                   help="Resume an interrupted run: skip steps whose outputs "
                        "already exist (per-step markers), and within the "
                        "SAM->BAM and salmon-quant steps, skip individual "
                        "files/cells already completed.")
    args = p.parse_args()

    # Validate CLS list files exist (defaults point to /mnt/hdd2)
    for name, path in [("--cls1", args.cls1), ("--cls2", args.cls2),
                       ("--cls3", args.cls3)]:
        if not os.path.exists(path):
            p.error(f"{name} file not found: {path}")

    R1 = args.input
    R2 = R1.replace("R1", "R2")
    if R2 == R1:
        plog("WARNING: R2 filename == R1 (no 'R1' substring). Check -i.")
    ncpu = args.threads
    plog(f"Threads: {ncpu}")
    plog(f"R1: {R1}")
    plog(f"R2: {R2}")

    sample = os.path.basename(R1).split(".")[0]
    workdir = f"TAS-Seq_{sample}"
    os.makedirs(workdir, exist_ok=True)

    tmp_R1       = os.path.join(workdir, "extracted_R1.fastq.gz")
    extracted_R2 = os.path.join(workdir, "extracted_R2.fastq.gz")
    trimmed_R2   = os.path.join(workdir, "trimmed_R2.fastq.gz")
    bam          = os.path.join(workdir, f"{sample}.bam")
    dedup_sam    = os.path.join(workdir, f"{sample}.dedup.sam")
    demux_dir    = os.path.join(workdir, "sam_demultiplexed")
    bam_split    = os.path.join(workdir, "bwa_mem_split")
    salmon_dir   = os.path.join(workdir, "salmon")
    reads_txt    = os.path.join(workdir, "R2-extracted-read.txt")
    fastp_html   = os.path.join(workdir, "fastp.html")
    fastp_json   = os.path.join(workdir, "fastp.json")
    summary_tsv  = os.path.join(workdir, f"summary_{sample}.tsv")

    if args.resume and os.path.exists(summary_tsv):
        plog(f"--resume: {summary_tsv} already exists -> pipeline already "
             f"completed for {sample}. Nothing to do.")
        return

    if args.resume:
        plog("--resume: will skip steps whose outputs are already present")

    # ── Step 1: extract (embedded bd-extract, CLS-list mode) ──────────────────
    _d = step_done_data(workdir, "extract", [extracted_R2]) if args.resume else None
    if _d is not None:
        plog("=== Step 1/7: extract - SKIPPED (resume, extracted_R2 present) ===")
        ext_stats = _d["ext_stats"]
    else:
        plog("=== Step 1/7: extract (CB+UMI extract & correct, embedded) ===")
        corrector = PartCorrector(args.cls1, args.cls2, args.cls3, args.mismatch)
        # R1 output is discarded downstream -> pass None to skip writing it entirely
        ext_stats = run_extract(R1, R2, None, extracted_R2,
                                corrector, "enhanced", args.linker_mm,
                                ncpu, args.chunk_size)
        mark_step_done(workdir, "extract", {"ext_stats": ext_stats})

    # ── Step 2: fastp poly-X trim ─────────────────────────────────────────────
    _d = step_done_data(workdir, "fastp", [trimmed_R2]) if args.resume else None
    if _d is not None:
        plog("=== Step 2/7: fastp - SKIPPED (resume, trimmed_R2 present) ===")
    else:
        plog("=== Step 2/7: fastp (poly-X trim) ===")
        run_cmd(f"fastp --trim_poly_x -w {min(ncpu, 16)} "
                f"-i {extracted_R2} -o {trimmed_R2} "
                f"--adapter_fasta={args.adapter} "
                f"-l {args.min_len} -h {fastp_html} -j {fastp_json}")
        mark_step_done(workdir, "fastp")

    # ── Step 3: bwa mem -> sorted BAM ─────────────────────────────────────────
    _d = step_done_data(workdir, "bwa", [bam, bam + ".bai"]) if args.resume else None
    if _d is not None:
        plog("=== Step 3/7: bwa mem -> BAM - SKIPPED (resume, BAM present) ===")
    else:
        plog("=== Step 3/7: bwa mem -> BAM ===")
        run_cmd(f"bwa mem -t {ncpu} {args.ref} {trimmed_R2} | "
                f"samtools view -@ {ncpu} -b | "
                f"samtools sort -@ {ncpu} -m 3G > {bam}")
        run_cmd(f"samtools index {bam}")
        mark_step_done(workdir, "bwa")

    # ── Step 4: directional UMI dedup (self-implemented, parallel) ────────────
    _d = step_done_data(workdir, "dedup", [dedup_sam]) if args.resume else None
    if _d is not None:
        plog("=== Step 4/7: dedup - SKIPPED (resume, dedup.sam present) ===")
        dedup_in, dedup_out = _d["dedup_in"], _d["dedup_out"]
    else:
        plog("=== Step 4/7: directional UMI dedup (parallel) ===")
        dedup_in, dedup_out = parallel_dedup(bam, dedup_sam, threads=ncpu)
        mark_step_done(workdir, "dedup", {"dedup_in": dedup_in, "dedup_out": dedup_out})

    # ── Step 5: demultiplex SAM by barcode ────────────────────────────────────
    _d = step_done_data(workdir, "demux", [demux_dir]) if args.resume else None
    if _d is not None:
        plog("=== Step 5/7: demultiplex - SKIPPED (resume, demux dir present) ===")
    else:
        plog("=== Step 5/7: demultiplex SAM by barcode ===")
        demultiplex_sam(dedup_sam, demux_dir)
        mark_step_done(workdir, "demux")

    # ── Step 6: SAM -> BAM (parallel; file-level resume) ──────────────────────
    plog("=== Step 6/7: SAM -> BAM (parallel) ===")
    sam_to_bam_all(demux_dir, args.ref, bam_split, ncpu, resume=args.resume)

    # ── Step 7: salmon quant (parallel, -p 2; file-level resume) ──────────────
    plog("=== Step 7/7: salmon quant (parallel) ===")
    salmon_jobs = max(1, ncpu // 2)   # 2 threads per job
    salmon_all(bam_split, args.ref, salmon_dir, salmon_jobs, resume=args.resume)

    # R2 extracted read count
    if not (args.resume and os.path.exists(reads_txt)):
        r2_read_count(trimmed_R2, reads_txt)

    # ── Summary (sample-level, machine-readable TSV) ──────────────────────────
    plog("=== Writing summary ===")
    n_cells = len([d for d in glob.glob(os.path.join(salmon_dir, "*"))
                   if os.path.isdir(d)])
    total_reads   = ext_stats.get("total", 0)
    pass_reads    = ext_stats.get("pass_", 0)
    corrected     = ext_stats.get("corrected", 0)
    nolinker      = ext_stats.get("nolinker", 0)
    discarded     = ext_stats.get("fail", 0)
    # UMI conversion rate = unique molecules (after dedup) / reads entering dedup
    umi_conv = (dedup_out / dedup_in) if dedup_in else 0.0
    pct = lambda a, b: (100.0 * a / b) if b else 0.0

    rows = [
        ("sample",                    sample),
        ("total_reads",              total_reads),
        ("cb_assigned_reads",        pass_reads),
        ("cb_assigned_pct",          f"{pct(pass_reads, total_reads):.2f}"),
        ("cb_corrected_reads",       corrected),
        ("no_linker_reads",          nolinker),
        ("discarded_reads",          discarded),
        ("mapped_reads_into_dedup",  dedup_in),
        ("unique_molecules",         dedup_out),
        ("umi_conversion_rate",      f"{umi_conv:.4f}"),
        ("umi_conversion_pct",       f"{pct(dedup_out, dedup_in):.2f}"),
        ("duplication_rate_pct",     f"{100 - pct(dedup_out, dedup_in):.2f}"),
        ("n_cells",                  n_cells),
        ("mean_molecules_per_cell",  f"{(dedup_out / n_cells):.1f}" if n_cells else "0"),
    ]
    with open(summary_tsv, "w") as out:
        out.write("metric\tvalue\n")
        for k, v in rows:
            out.write(f"{k}\t{v}\n")
    plog(f"Summary written: {summary_tsv}")
    plog(f"  UMI conversion rate: {100*umi_conv:.1f}%  "
         f"({dedup_in:,} reads -> {dedup_out:,} molecules)")

    # ── Cleanup ───────────────────────────────────────────────────────────────
    if not args.keep_intermediates:
        plog("=== Cleanup: removing intermediate bam/sam/gz "
             "(keeping salmon/, fastp, reads) ===")
        for f in [tmp_R1, extracted_R2, trimmed_R2,
                  bam, bam + ".bai", dedup_sam]:
            if os.path.exists(f):
                os.remove(f)
        for d in [demux_dir, bam_split]:
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)
        plog("Cleanup done. Kept: "
             f"{os.path.basename(salmon_dir)}/, fastp.html/json, "
             f"{os.path.basename(reads_txt)}")
    else:
        plog("Keeping all intermediates (--keep-intermediates)")

    plog(f"=== Pipeline finished. Output in {workdir}/ ===")

if __name__ == "__main__":
    main()

