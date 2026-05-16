#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
a3m_mask_query_protected_multimer.py

Mask specific query-referenced columns in a multimer A3M alignment
while NEVER modifying the FIRST sequence (query).

Multimer logic:
- Multimer FASTA has chains separated by ':' (e.g. A:B:C)
- A3M sequences are concatenated (no separators)
- Chain lengths are inferred from the FASTA
- Mask chain-local positions for a selected chain (A, B, C, ...)

Positions:
- 1-based
- Query numbering
- Gaps '-' and lowercase inserts are ignored for counting
"""

import argparse
import random
from typing import Dict, List, Optional, Tuple

# ---------------------------- A3M I/O ---------------------------- #

def parse_a3m(path: str) -> List[Tuple[str, str]]:
    records = []
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append((header, "".join(seq_chunks)))
    if not records:
        raise ValueError(f"No sequences found in {path}")
    return records


def write_a3m(records: List[Tuple[str, str]], path: str):
    with open(path, "w") as f:
        for h, s in records:
            f.write(h + "\n")
            for i in range(0, len(s), 120):
                f.write(s[i:i+120] + "\n")

# ---------------------------- FASTA ---------------------------- #

def read_first_fasta_sequence(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line)
    if not seq:
        raise ValueError(f"No FASTA sequence found in {path}")
    return "".join(seq)


def infer_chain_lengths_from_fasta(path: str, sep=":") -> List[int]:
    seq = read_first_fasta_sequence(path)
    parts = seq.split(sep)
    if len(parts) < 2:
        raise ValueError(
            f"Expected multimer FASTA with '{sep}' separators, got single chain."
        )
    return [len(p) for p in parts]


def chain_offsets(lengths: List[int]) -> List[Tuple[int, int]]:
    offsets = []
    s = 0
    for L in lengths:
        offsets.append((s, s + L))
        s += L
    return offsets

# ---------------------------- ranges ---------------------------- #

def parse_ranges(spec: str) -> List[int]:
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a, b = int(a), int(b)
            if a > b:
                a, b = b, a
            out.update(range(a, b + 1))
        else:
            out.add(int(part))
    return sorted(out)

# ---------------------------- query mapping ---------------------------- #

def build_query_map_chain(query_seq: str, start: int, end: int) -> Dict[int, List[int]]:
    """
    Map chain-local query positions (1-based, A–Z only)
    to GLOBAL alignment indices.
    """
    qmap = {}
    qpos = 0
    for i in range(start, end):
        ch = query_seq[i]
        if "A" <= ch <= "Z":
            qpos += 1
            qmap.setdefault(qpos, []).append(i)
    return qmap

# ---------------------------- masking ---------------------------- #

def mask_records(
    records: List[Tuple[str, str]],
    mutant_qpos: List[int],
    channel_qpos: List[int],
    channel_fraction: float,
    masking_char: str,
    chain_lengths: List[int],
    mask_chain_index: int,
    rng: Optional[random.Random] = None,
) -> Tuple[List[Tuple[str, str]], List[int], List[int]]:
    rng = rng or random

    query_header, query_seq = records[0]
    total_len = sum(chain_lengths)

    if len(query_seq) != total_len:
        raise ValueError(
            f"A3M query length ({len(query_seq)}) != sum of chain lengths ({total_len})."
        )

    offsets = chain_offsets(chain_lengths)
    if mask_chain_index >= len(offsets):
        raise ValueError("Requested chain index exceeds number of chains.")

    start, end = offsets[mask_chain_index]
    qmap = build_query_map_chain(query_seq, start, end)

    align_cols_mutant = []
    align_cols_channel = []
    applied = []
    missing = []

    for qp in sorted(set(mutant_qpos + channel_qpos)):
        if qp in qmap:
            applied.append(qp)
        else:
            missing.append(qp)

    for qp in mutant_qpos:
        if qp in qmap:
            align_cols_mutant.extend(qmap[qp])
    for qp in channel_qpos:
        if qp in qmap:
            align_cols_channel.extend(qmap[qp])

    align_cols_mutant = sorted(set(align_cols_mutant))
    align_cols_channel = sorted(set(align_cols_channel))

    def mask_char(ch):
        if "A" <= ch <= "Z":
            return masking_char
        return ch

    new_records = []
    for i, (hdr, seq) in enumerate(records):
        if i == 0:
            new_records.append((hdr, seq))
            continue

        s = list(seq)

        for col in align_cols_mutant:
            s[col] = mask_char(s[col])

        if align_cols_channel:
            if channel_fraction >= 1.0:
                for col in align_cols_channel:
                    s[col] = mask_char(s[col])
            elif channel_fraction > 0:
                for col in align_cols_channel:
                    if rng.random() < channel_fraction:
                        s[col] = mask_char(s[col])

        new_records.append((hdr, "".join(s)))

    return new_records, applied, missing

# ---------------------------- CLI ---------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Mask multimer A3M columns by chain-local query positions (query preserved)."
    )
    ap.add_argument("--input-a3m", required=True)
    ap.add_argument("--output-a3m", required=True)
    ap.add_argument("--multimer-fasta", required=True)
    ap.add_argument("--mask-chain", default="A")
    ap.add_argument("--chain-sep", default=":")
    ap.add_argument("--masking-char", default="X")

    ap.add_argument("--mutant-ranges", default="")
    ap.add_argument("--channel-masking", default="")
    ap.add_argument("--channel-mask-percent", type=float, default=1.0)
    ap.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible stochastic channel masking.",
    )

    args = ap.parse_args()

    if not (0.0 <= args.channel_mask_percent <= 1.0):
        raise SystemExit("--channel-mask-percent must be in [0,1]")

    if len(args.masking_char) != 1:
        raise SystemExit("--masking-char must be a single character")

    chain_letter = args.mask_chain.upper()
    if not ("A" <= chain_letter <= "Z"):
        raise SystemExit("--mask-chain must be A, B, C, ...")

    mask_chain_index = ord(chain_letter) - ord("A")

    records = parse_a3m(args.input_a3m)
    rng = random.Random(args.seed) if args.seed is not None else None

    chain_lengths = infer_chain_lengths_from_fasta(
        args.multimer_fasta, sep=args.chain_sep
    )

    mutant_qpos = parse_ranges(args.mutant_ranges) if args.mutant_ranges else []
    channel_qpos = parse_ranges(args.channel_masking) if args.channel_masking else []

    new_records, applied, missing = mask_records(
        records,
        mutant_qpos,
        channel_qpos,
        args.channel_mask_percent,
        args.masking_char,
        chain_lengths,
        mask_chain_index,
        rng=rng,
    )

    write_a3m(new_records, args.output_a3m)

    offsets = chain_offsets(chain_lengths)
    s, e = offsets[mask_chain_index]

    print(f"# Query preserved: {records[0][0]}")
    print(f"# Chains: {len(chain_lengths)} lengths={chain_lengths}")
    print(f"# Masking chain {chain_letter} region [{s}:{e}]")
    print(f"# Mutant positions (chain-local): {mutant_qpos}")
    print(f"# Channel positions (chain-local): {channel_qpos}")
    if args.seed is not None:
        print(f"# Random seed: {args.seed}")
    print(f"# Applied positions: {applied}")
    if missing:
        print(f"# WARNING missing positions (ignored): {missing}")
    print(f"# Output written: {args.output_a3m}")

if __name__ == "__main__":
    main()
