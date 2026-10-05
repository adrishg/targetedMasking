"""Validated masking of ColabFold A3M match columns (standard library only).

Stored blocks, not physical copies, define column offsets. Lowercase insertions
are retained verbatim and never advance a match-column counter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import random
import re
import string


_INSERTIONS = str.maketrans("", "", string.ascii_lowercase)
_SEQUENCE = re.compile(r"[A-Za-z-]+\Z")
_META = re.compile(r"#([1-9][0-9]*(?:,[1-9][0-9]*)*)\t([1-9][0-9]*(?:,[1-9][0-9]*)*)\Z")


def aligned(sequence):
    """Independent projection used by the verifier, not by the masking loop."""
    return sequence.translate(_INSERTIONS)


def unit_name(index):
    name = ""
    index += 1
    while index:
        index, digit = divmod(index - 1, 26)
        name = chr(65 + digit) + name
    return name


@dataclass(frozen=True)
class Document:
    metadata: str | None
    records: tuple[tuple[str, str], ...]

    def text(self):
        # One sequence per line works with old and current ColabFold readers.
        prefix = self.metadata + "\n" if self.metadata is not None else ""
        return prefix + "".join(h + "\n" + s + "\n" for h, s in self.records)


@dataclass(frozen=True)
class Layout:
    lengths: tuple[int, ...]
    copies: tuple[int, ...]
    queries: tuple[str, ...]

    @property
    def width(self):
        return sum(self.lengths)

    @property
    def kind(self):
        if sum(self.copies) == 1:
            return "monomer"
        return "homomer" if len(set(self.queries)) == 1 else "heteromer"

    @property
    def blocks(self):
        offset = 0
        result = {}
        for i, (length, copies) in enumerate(zip(self.lengths, self.copies)):
            result[unit_name(i)] = (offset, length, copies)
            offset += length
        return result


def parse_a3m_text(text):
    metadata = None
    records = []
    header = None
    chunks = []

    def finish():
        if header is not None:
            seq = "".join(chunks)
            if not seq:
                raise ValueError(f"Empty sequence under {header!r}")
            records.append((header, seq))

    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith("#"):
            if metadata is not None or header is not None or not _META.fullmatch(line):
                raise ValueError(
                    f"Line {number}: expected one initial ColabFold header "
                    "'#lengths<TAB>copy_counts'; extra comments/headers are unsupported."
                )
            metadata = line  # Preserve the exact header, including actual TAB.
        elif line.startswith(">"):
            finish()
            if not line[1:].strip():
                raise ValueError(f"Line {number}: empty sequence identifier")
            header, chunks = line, []
        else:
            if header is None:
                raise ValueError(f"Line {number}: sequence before FASTA header")
            seq = line.strip()
            if not _SEQUENCE.fullmatch(seq):
                raise ValueError(
                    f"Line {number}: unsupported sequence characters. Use A3M "
                    "letters and '-'; dots, separators, embedded whitespace and NULs "
                    "must not be silently interpreted as match columns."
                )
            chunks.append(seq)
    finish()
    if not records:
        raise ValueError("No A3M sequence records")
    doc = Document(metadata, tuple(records))
    validate_document(doc)
    return doc


def read_a3m(path):
    return parse_a3m_text(Path(path).read_text(encoding="utf-8"))


def validate_document(doc):
    if not doc.records:
        raise ValueError("No A3M sequence records")
    query = doc.records[0][1]
    if not query or any(not ("A" <= c <= "Z") for c in query):
        raise ValueError("The first query must be ungapped uppercase residues, without insertions")
    if doc.metadata is None:
        lengths, copies = (len(query),), (1,)
    else:
        match = _META.fullmatch(doc.metadata)
        if match is None:
            raise ValueError("Malformed ColabFold metadata header")
        lengths = tuple(map(int, match[1].split(",")))
        copies = tuple(map(int, match[2].split(",")))
        if len(lengths) != len(copies):
            raise ValueError("Header lengths and copy counts have different sizes")
    width = sum(lengths)  # Deliberately NOT sum(length * copies).
    if len(query) != width:
        raise ValueError(f"Query has {len(query)} columns; header declares {width}")
    for row, (header, seq) in enumerate(doc.records, 1):
        if not header.startswith(">") or "\n" in header or "\r" in header:
            raise ValueError(f"Row {row}: invalid FASTA header")
        if not _SEQUENCE.fullmatch(seq):
            raise ValueError(f"Row {row} {header}: unsupported A3M characters")
        observed = len(aligned(seq))
        if observed != width:
            raise ValueError(
                f"Row {row} {header}: observed {observed} A3M match states; expected {width}. "
                "Separate chain alignments must be supplied as separate files, "
                "not mixed into a combined A3M."
            )
    queries = []
    offset = 0
    for length in lengths:
        queries.append(query[offset:offset + length])
        offset += length
    return Layout(lengths, copies, tuple(queries))


def read_fasta_chains(path, sep=":"):
    sequences, chunks = [], []
    seen_header = False
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seen_header:
                if not chunks:
                    raise ValueError("Empty FASTA record")
                sequences.extend("".join(chunks).split(sep))
            chunks = []
            seen_header = True
        else:
            chunks.append(line)
    if chunks:
        sequences.extend("".join(chunks).split(sep))
    elif seen_header:
        raise ValueError("Empty FASTA record")
    if not sequences or any(not re.fullmatch(r"[A-Z]+", s) for s in sequences):
        raise ValueError("FASTA must contain nonempty uppercase ungapped protein chains")
    return sequences


def check_fasta(layout, chains):
    expected = [q for q, count in zip(layout.queries, layout.copies) for _ in range(count)]
    if chains != expected:
        raise ValueError(
            "FASTA sequences/order/copy counts disagree with the A3M query and header. "
            "Supply all physical copies in stored-block order; lengths alone are insufficient."
        )


def parse_ranges(spec, max_position=None):
    positions = set()
    if not spec.strip():
        return []
    for token in spec.split(","):
        match = re.fullmatch(r"\s*([1-9][0-9]*)(?:\s*-\s*([1-9][0-9]*))?\s*", token)
        if not match:
            raise ValueError(f"Invalid 1-based range: {token!r}")
        first, last = int(match[1]), int(match[2] or match[1])
        first, last = min(first, last), max(first, last)
        if max_position is not None and last > max_position:
            raise ValueError(f"Position {last} exceeds block length {max_position}")
        if last - first > 1_000_000:
            raise ValueError("Range is too large for a protein alignment")
        positions.update(range(first, last + 1))
    return sorted(positions)


def resolve_targets(layout, specs, all_copies=False):
    """Resolve repeated 'BLOCK:1,4-7' requests to zero-based match columns."""
    columns = set()
    for spec in specs:
        if ":" not in spec:
            raise ValueError(f"Expected BLOCK:POSITIONS, got {spec!r}")
        label, ranges = spec.split(":", 1)
        label = label.strip().upper()
        if label not in layout.blocks:
            raise ValueError(f"Unknown stored block {label!r}; available: {', '.join(layout.blocks)}")
        offset, length, copies = layout.blocks[label]
        positions = parse_ranges(ranges, length)
        if not positions:
            raise ValueError(f"No positions specified for block {label}")
        if copies > 1 and not all_copies:
            raise ValueError(
                f"Block {label} is shared by {copies} copies. Use --all-copies to mask "
                "the shared MSA. A compact A3M cannot target just one physical copy."
            )
        columns.update(offset + p - 1 for p in positions)
    return columns


def protected_rows(doc, layout):
    """Protect global query and identified per-block queries in heteromer MSAs.

    A block query is identified by BOTH its query ID and its exact padded query
    sequence. Identical homologs with other IDs remain eligible for masking.
    Duplicate monomer/homomer homolog rows retain the historical behavior.
    """
    protected = {0}
    if len(layout.lengths) < 2:
        return protected
    ids = doc.records[0][0][1:].split("\t")
    if len(ids) != len(layout.lengths):
        raise ValueError(
            "Combined multimer query header must provide one TAB-separated ID "
            "per stored block, so per-block query rows can be identified safely."
        )
    if len(set(ids)) != len(ids) or any(not x for x in ids):
        raise ValueError("Combined multimer query IDs must be nonempty and distinct")
    query_records = {doc.records[0]}
    offset = 0
    for identifier, length, query in zip(ids, layout.lengths, layout.queries):
        padded = "-" * offset + query + "-" * (layout.width - offset - length)
        query_records.add((">" + identifier, padded))
        offset += length
    for row, record in enumerate(doc.records):
        if record in query_records:
            protected.add(row)
    return protected


def validate_request(layout, deterministic, stochastic, fraction):
    if not 0 <= fraction <= 1:
        raise ValueError("Mask fraction must be finite and in [0, 1]")
    for col in deterministic | stochastic:
        if not isinstance(col, int) or not 0 <= col < layout.width:
            raise ValueError(f"Target column {col!r} is outside the alignment")


def mask_document(doc, deterministic, stochastic=(), fraction=1.0, seed=7):
    layout = validate_document(doc)
    deterministic, stochastic = set(deterministic), set(stochastic)
    validate_request(layout, deterministic, stochastic, fraction)
    protected = protected_rows(doc, layout)
    rng = random.Random(seed)
    output = []
    for row, (header, sequence) in enumerate(doc.records):
        if row in protected:
            output.append((header, sequence))
            continue
        chars = list(sequence)
        col = -1
        for raw, char in enumerate(sequence):
            if "a" <= char <= "z":
                continue
            col += 1  # Deletion gaps DO occupy match columns.
            if char == "-":
                continue
            selected = col in deterministic
            if not selected and col in stochastic:
                selected = fraction == 1 or (fraction > 0 and rng.random() < fraction)
            if selected:
                chars[raw] = "X"
        output.append((header, "".join(chars)))
    result = Document(doc.metadata, tuple(output))
    report = verify_masking(doc, result, deterministic, stochastic, fraction)
    report["seed"] = seed
    return result, report


def verify_masking(original, output, deterministic, stochastic=(), fraction=1.0):
    """Check coordinates by independently stripping insertions from both MSAs.

    Does not trust the masker's character-index map or only test equal lengths.
    At fractional masking, validates allowed coordinates and reports coverage;
    it does not claim an exact fraction or reproduce the RNG draw history.
    """
    layout = validate_document(original)
    if validate_document(output) != layout:
        raise ValueError("Output layout/query differs from the original")
    if original.metadata != output.metadata:
        raise ValueError("ColabFold metadata header changed")
    if len(original.records) != len(output.records):
        raise ValueError("Record count changed")
    deterministic, stochastic = set(deterministic), set(stochastic)
    validate_request(layout, deterministic, stochastic, fraction)
    allowed = deterministic | (stochastic if fraction > 0 else set())
    required = deterministic | (stochastic if fraction == 1 else set())
    protected = protected_rows(original, layout)
    coverage = {c: {"residues": 0, "gaps": 0, "preexisting_X": 0,
                    "newly_masked": 0, "unmasked_residues": 0}
                for c in sorted(deterministic | stochastic)}
    changed_rows = 0
    for row, ((ha, a), (hb, b)) in enumerate(zip(original.records, output.records)):
        prefix = f"Row {row + 1} {ha}"
        if ha != hb:
            raise ValueError(f"{prefix}: sequence header/order changed")
        if row in protected:
            if a != b:
                raise ValueError(f"{prefix}: protected query changed")
            continue
        if len(a) != len(b):
            raise ValueError(f"{prefix}: raw sequence length changed")
        # Position-by-position check also preserves the LOCATION of insertions.
        for x, y in zip(a, b):
            if x != y and (not "A" <= x <= "Z" or y != "X"):
                raise ValueError(f"{prefix}: insertion, gap or non-mask character changed")
        aa, bb = aligned(a), aligned(b)
        changed_rows += int(aa != bb)
        for col, (x, y) in enumerate(zip(aa, bb)):
            if x != y and col not in allowed:
                raise ValueError(f"{prefix}: off-target change at aligned column {col + 1}")
            if col in required and x != "-" and y != "X":
                raise ValueError(f"{prefix}: missed mask at aligned column {col + 1}")
            if col in coverage:
                stats = coverage[col]
                stats["gaps"] += int(x == "-")
                stats["residues"] += int(x != "-")
                stats["preexisting_X"] += int(x == "X")
                stats["newly_masked"] += int(x != y)
                stats["unmasked_residues"] += int(x != "-" and y != "X")
    blocks = []
    site_report = []
    for name, (offset, length, copies) in layout.blocks.items():
        blocks.append({"block": name, "length": length, "copies": copies,
                       "first_aligned_column": offset + 1})
        for col, stats in coverage.items():
            if offset <= col < offset + length:
                site_report.append({"block": name, "position": col - offset + 1,
                                    "aligned_column": col + 1, **stats})
    return {
        "verification": "passed", "kind": layout.kind,
        "metadata": original.metadata, "blocks": blocks,
        "records": len(original.records), "aligned_columns": layout.width,
        "protected_query_rows_1based": [i + 1 for i in sorted(protected)],
        "changed_rows": changed_rows,
        "newly_masked_residues": sum(s["newly_masked"] for s in coverage.values()),
        "outside_target_changes": 0, "insertions_and_gaps_preserved": True,
        "headers_and_record_order_preserved": True,
        "deterministic_columns_1based": [c + 1 for c in sorted(deterministic)],
        "stochastic_columns_1based": [c + 1 for c in sorted(stochastic)],
        "stochastic_fraction": fraction, "coverage": site_report,
    }
