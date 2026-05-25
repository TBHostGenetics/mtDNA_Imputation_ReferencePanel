#!/usr/bin/env python3
"""
msa_to_vcf.py — Memory-efficient MSA → VCF converter for large mtDNA panels.

Problem with the naive approach:
    58,808 sequences × 16,569 sites × 1 byte = ~975 MB just for the sequences,
    plus Python object overhead this easily exceeds 8–16 GB and gets OOM-killed.

Solution:
    1. Single pass through the FASTA to write a raw byte matrix to a memory-
       mapped file on disk (numpy memmap). Each row = one sequence, each column
       = one alignment site. Disk is used as scratch; only one row at a time
       is in RAM during the write pass.
    2. Iterate column-by-column over the memmap (columns are read in small
       chunks). RAM usage stays ~constant regardless of MSA size.
    3. Write VCF lines as we go — no VCF line accumulation in memory.

Memory use: O(n_sequences) for one column vector at a time (~60 KB for 58k seqs).
Disk use:   ~975 MB temporary memmap file (deleted on exit).

Usage:
    python3 msa_to_vcf.py \\
        --msa  filtered/filtered_msa.fasta \\
        --out  data/vcf/mtdna_raw.vcf \\
        --rcrs NC_012920.1 \\
        [--tmp  /path/to/scratch]   # default: same dir as --out
        [--chunk 500]               # columns processed per iteration
"""

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Base coding
# STANDARD bases stored as their ASCII byte values for fast numpy comparison
# ---------------------------------------------------------------------------
A = ord('A'); C = ord('C'); G = ord('G'); T = ord('T')
N = ord('N')   # missing / gap / ambiguous sentinel
GAP1 = ord('-'); GAP2 = ord('.'); STANDARD_SET = {A, C, G, T}


def code_byte(b: int) -> int:
    """Map a raw ASCII byte to a coded byte (standard base or N)."""
    if b == A or b == C or b == G or b == T:
        return b
    return N


# Vectorised version for a whole numpy array column
_LUT = np.full(256, N, dtype=np.uint8)
for _b in (A, C, G, T):
    _LUT[_b] = _b
# Also handle lowercase
for _b, _B in ((ord('a'), A), (ord('c'), C), (ord('g'), G), (ord('t'), T)):
    _LUT[_b] = _B


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--msa",   required=True,  help="Filtered MSA (FASTA)")
    p.add_argument("--out",   required=True,  help="Output VCF path")
    p.add_argument("--rcrs",  default="NC_012920.1", help="rCRS accession")
    p.add_argument("--tmp",   default=None,
                   help="Directory for temporary memmap file (default: output dir)")
    p.add_argument("--chunk", type=int, default=500,
                   help="Columns to process per iteration (default: 500)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Pass 1: stream FASTA → memmap matrix
# ---------------------------------------------------------------------------
def build_memmap(fasta_path: str, tmp_dir: str):
    """
    Stream the FASTA file and write a uint8 memmap of shape (n_seqs, aln_len).
    Returns (memmap_array, sample_ids, rcrs_row_index, aln_len).
    """
    print(f"[msa_to_vcf] Pass 1: scanning FASTA to determine dimensions...",
          flush=True)
    t0 = time.time()

    # --- Quick scan for dimensions (seq count + alignment length) ---
    n_seqs = 0
    aln_len = None
    with open(fasta_path) as fh:
        seq_buf = []
        current_id = None
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(seq_buf)
                    if aln_len is None:
                        aln_len = len(seq)
                    elif len(seq) != aln_len:
                        sys.exit(f"[ERROR] Sequence '{current_id}' has length "
                                 f"{len(seq)}, expected {aln_len}. MSA not aligned.")
                    n_seqs += 1
                    seq_buf = []
                current_id = line[1:].split()[0]
            else:
                seq_buf.append(line)
        # last record
        if current_id is not None:
            seq = "".join(seq_buf)
            if aln_len is None:
                aln_len = len(seq)
            n_seqs += 1

    if n_seqs == 0 or aln_len is None:
        sys.exit("[ERROR] No sequences found in FASTA.")

    print(f"[msa_to_vcf]   Sequences : {n_seqs:,}", flush=True)
    print(f"[msa_to_vcf]   Aln length: {aln_len:,} columns", flush=True)
    print(f"[msa_to_vcf]   Matrix    : {n_seqs * aln_len / 1e9:.2f} GB on disk",
          flush=True)

    # --- Allocate memmap ---
    mmap_path = os.path.join(tmp_dir, "msa_matrix.mmap")
    print(f"[msa_to_vcf]   Memmap at : {mmap_path}", flush=True)
    mat = np.lib.format.open_memmap(
        mmap_path, mode='w+', dtype=np.uint8, shape=(n_seqs, aln_len)
    )

    # --- Pass 2: fill memmap row by row ---
    print(f"[msa_to_vcf] Pass 2: writing sequences to memmap...", flush=True)
    sample_ids = []
    row = 0
    with open(fasta_path) as fh:
        seq_buf = []
        current_id = None
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seq_bytes = np.frombuffer(
                        "".join(seq_buf).upper().encode(), dtype=np.uint8
                    )
                    mat[row, :] = _LUT[seq_bytes]
                    row += 1
                    seq_buf = []
                    if row % 5000 == 0:
                        print(f"[msa_to_vcf]   Written {row:,}/{n_seqs:,} sequences...",
                              flush=True)
                current_id = line[1:].split()[0]
                sample_ids.append(current_id)
            else:
                seq_buf.append(line)
        # last record
        if current_id is not None:
            seq_bytes = np.frombuffer(
                "".join(seq_buf).upper().encode(), dtype=np.uint8
            )
            mat[row, :] = _LUT[seq_bytes]

    mat.flush()
    elapsed = time.time() - t0
    print(f"[msa_to_vcf] Memmap built in {elapsed:.1f}s", flush=True)
    return mat, sample_ids, mmap_path, aln_len


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    out_dir = str(Path(args.out).parent)
    tmp_dir = args.tmp if args.tmp else out_dir
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    mmap_path = None
    try:
        # ----------------------------------------------------------------
        # Build memmap
        # ----------------------------------------------------------------
        mat, sample_ids, mmap_path, aln_len = build_memmap(args.msa, tmp_dir)
        n_seqs = len(sample_ids)

        # ----------------------------------------------------------------
        # Find rCRS row
        # ----------------------------------------------------------------
        rcrs_row = None
        rcrs_actual_id = None
        rcrs_prefix = args.rcrs.split(".")[0]
        for i, sid in enumerate(sample_ids):
            if sid == args.rcrs or sid.split(".")[0] == rcrs_prefix:
                rcrs_row = i
                rcrs_actual_id = sid
                break
        if rcrs_row is None:
            sys.exit(f"[ERROR] rCRS '{args.rcrs}' not found in MSA.")
        print(f"[msa_to_vcf] rCRS found: {rcrs_actual_id} (row {rcrs_row})",
              flush=True)

        # Sample IDs exclude rCRS (rCRS defines REF, not a sample)
        sample_mask = np.ones(n_seqs, dtype=bool)
        sample_mask[rcrs_row] = False
        sample_ids_out = [sid for i, sid in enumerate(sample_ids)
                          if sample_mask[i]]
        n_samples = len(sample_ids_out)
        print(f"[msa_to_vcf] Samples (excl. rCRS): {n_samples:,}", flush=True)

        # ----------------------------------------------------------------
        # Write VCF
        # ----------------------------------------------------------------
        print(f"[msa_to_vcf] Writing VCF: {args.out}", flush=True)
        t0 = time.time()

        with open(args.out, "w") as vcf:
            # Header
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##reference={args.rcrs}\n")
            vcf.write('##contig=<ID=chrM,length=16569>\n')
            vcf.write('##INFO=<ID=AC,Number=A,Type=Integer,'
                      'Description="Allele count in samples">\n')
            vcf.write('##INFO=<ID=AN,Number=1,Type=Integer,'
                      'Description="Total allele count">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,'
                      'Description="Genotype">\n')
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                      + "\t".join(sample_ids_out) + "\n")

            rcrs_pos = 0        # 1-based physical position counter
            n_variant = 0
            n_invariant = 0
            n_rcrs_gap = 0

            # Process in column chunks for I/O efficiency
            chunk = args.chunk
            for col_start in range(0, aln_len, chunk):
                col_end = min(col_start + chunk, aln_len)
                block = mat[:, col_start:col_end]  # shape (n_seqs, chunk_width)

                for local_col in range(col_end - col_start):
                    col = block[:, local_col]       # shape (n_seqs,)
                    ref_base = col[rcrs_row]

                    # Always advance the rCRS position counter, even when the
                    # base is N. rCRS position 3107 is a known placeholder N
                    # (officially absent in the rCRS but retained in GenBank to
                    # preserve sequence length). Skipping the counter here caused
                    # every position after 3107 to be off by one, making real
                    # variant sites appear invariant in the VCF.
                    rcrs_pos += 1

                    if ref_base == N:
                        # rCRS base is N — cannot define REF allele, skip record.
                        # Counter already incremented so numbering stays correct.
                        n_rcrs_gap += 1
                        continue

                    # Sample bases (exclude rCRS row)
                    sample_col = col[sample_mask]   # shape (n_samples,)

                    # Unique non-REF, non-N bases = ALT alleles
                    unique = np.unique(sample_col)
                    alt_bases = sorted(
                        int(b) for b in unique
                        if b != ref_base and b != N
                    )

                    if not alt_bases:
                        n_invariant += 1
                        continue  # invariant — skip (BCFtools will also filter)

                    n_variant += 1
                    ref_char = chr(ref_base)
                    alt_chars = [chr(b) for b in alt_bases]
                    alt_str = ",".join(alt_chars)

                    # Allele index map: ref_base→0, alt[0]→1, alt[1]→2, N→"."
                    allele_map = np.full(256, 255, dtype=np.uint8)  # 255 = missing
                    allele_map[ref_base] = 0
                    for idx, ab in enumerate(alt_bases, 1):
                        allele_map[ab] = idx
                    allele_map[N] = 255

                    gt_indices = allele_map[sample_col]  # vectorised lookup

                    # Build GT strings
                    gt_strs = np.where(gt_indices == 255, ".", gt_indices.astype(str))

                    # AC and AN
                    an = int(np.sum(gt_indices != 255))
                    ac_list = [
                        str(int(np.sum(gt_indices == i)))
                        for i in range(1, len(alt_bases) + 1)
                    ]
                    info = f"AC={','.join(ac_list)};AN={an}"

                    # Assemble VCF line
                    gt_field = "\t".join(gt_strs.tolist())
                    vcf.write(
                        f"chrM\t{rcrs_pos}\t.\t{ref_char}\t{alt_str}\t"
                        f".\tPASS\t{info}\tGT\t{gt_field}\n"
                    )

                if (col_start // chunk) % 10 == 0:
                    pct = 100 * col_start / aln_len
                    print(f"[msa_to_vcf]   {pct:.0f}% ({col_start:,}/{aln_len:,} columns)...",
                          flush=True)

        elapsed = time.time() - t0
        print(f"\n[msa_to_vcf] VCF written in {elapsed:.1f}s", flush=True)
        print(f"[msa_to_vcf] rCRS positions  : {rcrs_pos:,}", flush=True)
        print(f"[msa_to_vcf] Variant sites   : {n_variant:,}", flush=True)
        print(f"[msa_to_vcf] Invariant sites : {n_invariant:,} (skipped)", flush=True)
        print(f"[msa_to_vcf] rCRS gap cols   : {n_rcrs_gap:,} (skipped)", flush=True)
        print(f"[msa_to_vcf] VCF output      : {args.out}", flush=True)

    finally:
        # Always clean up the temporary memmap file
        if mmap_path and os.path.exists(mmap_path):
            os.remove(mmap_path)
            print(f"[msa_to_vcf] Temporary memmap deleted: {mmap_path}", flush=True)


if __name__ == "__main__":
    main()
