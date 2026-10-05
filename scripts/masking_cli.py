"""Command-line interface for targetedMasking.py."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
import tempfile

if __package__:
    from .a3m_masking import (check_fasta, mask_document, read_a3m,
                             read_fasta_chains, resolve_targets,
                             validate_document, verify_masking)
else:
    from a3m_masking import (check_fasta, mask_document, read_a3m,
                            read_fasta_chains, resolve_targets,
                            validate_document, verify_masking)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_mode(layout, mode):
    if mode != "auto" and layout.kind != mode:
        raise ValueError(f"Expected a {mode} A3M, but this one describes a {layout.kind}")


def print_layout(layout, metadata, stream, label=None):
    """Explain stored blocks versus physical copies before selecting targets."""
    prefix = f"[{label}] " if label else ""
    print(f"{prefix}Identified {layout.kind}: {sum(layout.copies)} chain(s), "
          f"{len(layout.lengths)} stored alignment block(s), "
          f"{layout.width} aligned columns.", file=stream)
    if metadata is None:
        print(f"{prefix}No multimer metadata header: assuming one monomer; "
              "missing copy counts cannot be inferred from sequence alone.", file=stream)
    else:
        print(f"{prefix}A3M header: {metadata!r} (preserved)", file=stream)
    for name, (offset, length, copies) in layout.blocks.items():
        print(f"{prefix}Block {name}: {length} residues, {copies} copy/copies, "
              f"aligned columns {offset + 1}-{offset + length}.", file=stream)
        if copies > 1:
            print(f"{prefix}Block {name} is shared: masking it affects all {copies} copies "
                  "and requires --all-copies.", file=stream)


def targets(args, layout):
    full = list(args.mask)
    stochastic = list(args.stochastic_mask)
    if args.mutant_ranges:
        full.append(f"{args.mask_unit}:{args.mutant_ranges}")
    if args.channel_masking:
        stochastic.append(f"{args.mask_unit}:{args.channel_masking}")
    if not full and not stochastic:
        raise ValueError("Specify --mask BLOCK:POSITIONS or --mutant-ranges POSITIONS")
    return (resolve_targets(layout, full, args.all_copies),
            resolve_targets(layout, stochastic, args.all_copies))


def publish(jobs, report_path, report):
    """Stage and independently re-read every A3M before publishing new files.

    Never overwrite an input or existing output. Roll back files created by this
    invocation if publication fails. A machine crash is not a multi-file transaction.
    """
    paths = [Path(j["output"]).absolute() for j in jobs] + [Path(report_path).absolute()]
    resolved = [p.resolve() for p in paths]
    if len(set(resolved)) != len(paths):
        raise ValueError("Output/report paths must be distinct")
    for path in paths:
        if os.path.lexists(path):
            raise ValueError(f"Refusing to overwrite existing file: {path}")
        if not path.parent.is_dir():
            raise ValueError(f"Output directory does not exist: {path.parent}")
    staged, created = [], []
    try:
        for job, destination in zip(jobs, paths):
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                             dir=destination.parent, delete=False) as handle:
                stage = Path(handle.name)
                staged.append((stage, destination))
                handle.write(job["masked"].text())
            reread = read_a3m(stage)
            verify_masking(job["original"], reread, job["full"], job["stochastic"], job["fraction"])
            if reread != job["masked"]:
                raise ValueError("Serialized output differs from verified in-memory result")
            job["report"]["output_sha256"] = sha256(stage)
            job["report"]["serialized_output_verified"] = True
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=paths[-1].parent,
                                         delete=False) as handle:
            staged.append((Path(handle.name), paths[-1]))
            json.dump(report, handle, indent=2)
            handle.write("\n")
        for stage, destination in staged:
            os.link(stage, destination)  # Atomic creation; fails if destination now exists.
            created.append(destination)
    except BaseException:
        for path in reversed(created):
            path.unlink()
        raise
    finally:
        for path, _ in staged:
            path.unlink(missing_ok=True)


def prepare_job(args, source, full, stochastic, output=None, original=None):
    original = original or read_a3m(source)
    masked, report = mask_document(original, full, stochastic, args.channel_mask_percent, args.seed)
    report.update({"input": str(Path(source).resolve()), "input_sha256": sha256(source),
                   "output": str(Path(output).absolute()) if output else None})
    return dict(original=original, masked=masked, full=full, stochastic=stochastic,
                fraction=args.channel_mask_percent, output=output, report=report)


def separate_jobs(args):
    if args.input_a3m or args.output_a3m or args.verify_output or args.query_fasta:
        raise ValueError("Separate-file mode uses --chain-a3m and --output-dir; "
                         "do not mix combined-input/FASTA/verify-output options")
    sources = {}
    for spec in args.chain_a3m:
        name, sep, path = spec.partition("=")
        name = name.upper()
        if not sep or not name.isascii() or not name.isalpha() or not path or name in sources:
            raise ValueError(f"Expected distinct BLOCK=PATH assignments, got {spec!r}")
        sources[name] = path
    if len(sources) < 2:
        raise ValueError("Separate heteromer mode needs at least two subunit A3Ms")
    if not args.dry_run and not args.output_dir:
        raise ValueError("Separate-file mode requires --output-dir (an existing directory)")
    full_specs = list(args.mask)
    random_specs = list(args.stochastic_mask)
    if args.mutant_ranges:
        full_specs.append(f"{args.mask_unit}:{args.mutant_ranges}")
    if args.channel_masking:
        random_specs.append(f"{args.mask_unit}:{args.channel_masking}")
    if not full_specs and not random_specs:
        raise ValueError("Specify the subunit and positions with --mask BLOCK:POSITIONS")
    assigned = {name: [[], []] for name in sources}
    for group, specs in enumerate((full_specs, random_specs)):
        for spec in specs:
            name, sep, positions = spec.partition(":")
            name = name.upper()
            if not sep or name not in sources:
                raise ValueError(f"Unknown/missing subunit in {spec!r}")
            assigned[name][group].append("A:" + positions)
    docs = {name: read_a3m(path) for name, path in sources.items()}
    if len({doc.records[0][1] for doc in docs.values()}) < 2:
        raise ValueError("Separate heteromer inputs must contain at least two distinct queries")
    jobs = []
    stream = sys.stderr if args.dry_run else sys.stdout
    print(f"Identified {len(docs)} separate subunit A3Ms; cross-subunit pairing "
          "and complex copy counts are not encoded in these files.", file=stream)
    for name, source in sources.items():
        layout = validate_document(docs[name])
        check_mode(layout, "monomer")
        print_layout(layout, docs[name].metadata, stream, label=f"subunit {name}")
        full = resolve_targets(layout, assigned[name][0])
        stochastic = resolve_targets(layout, assigned[name][1])
        output = Path(args.output_dir) / f"{name}.masked.a3m" if args.output_dir else None
        job = prepare_job(args, source, full, stochastic, output, docs[name])
        job["report"]["subunit"] = name
        jobs.append(job)
    report = {"verification": "passed", "representation": "separate_subunit_a3ms",
              "pairing": "unchanged; no cross-subunit pairing is inferred or created",
              "files": [job["report"] for job in jobs]}
    return jobs, report


def main(mode="auto", argv=None):
    description = ("Auto-detect monomer/homomer/heteromer from ColabFold A3M metadata; "
                   "headerless inputs are treated as monomers. " if mode == "auto" else "")
    parser = argparse.ArgumentParser(description=description + f"Verified {mode} A3M masking; preserves query, "
                                     "headers, gaps, insertions and paired/unpaired row order.")
    parser.add_argument("--input-a3m")
    parser.add_argument("--output-a3m")
    parser.add_argument("--mask", action="append", default=[], metavar="BLOCK:POSITIONS",
                        help="Repeat for different stored blocks; 1-based local positions")
    parser.add_argument("--stochastic-mask", action="append", default=[], metavar="BLOCK:POSITIONS")
    parser.add_argument("--mask-unit", "--mask-chain", dest="mask_unit", default="A",
                        help="Stored block for legacy range options, NOT physical copy ID")
    parser.add_argument("--mutant-ranges", default="")
    parser.add_argument("--channel-masking", default="")
    parser.add_argument("--channel-mask-percent", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--masking-char", choices=["X"], default="X",
                        help="X preserves gap topology and paired/unpaired occupancy")
    parser.add_argument("--all-copies", action="store_true",
                        help="Acknowledge masking ALL copies sharing a compact alignment block")
    parser.add_argument("--query-fasta", "--multimer-fasta", dest="query_fasta")
    parser.add_argument("--chain-sep", default=":")
    parser.add_argument("--report-json", help="Default: OUTPUT.a3m.audit.json")
    parser.add_argument("--dry-run", action="store_true", help="Validate/mask in memory; write nothing")
    parser.add_argument("--verify-output", help="Audit an existing masked A3M against input; write nothing")
    if mode in ("heteromer", "auto"):
        parser.add_argument("--chain-a3m", action="append", default=[], metavar="BLOCK=PATH")
        parser.add_argument("--output-dir", help="Existing directory for separate subunit outputs")
    else:
        parser.set_defaults(chain_a3m=[], output_dir=None)
    args = parser.parse_args(argv)
    try:
        if not 0 <= args.channel_mask_percent <= 1:
            raise ValueError("--channel-mask-percent must be finite and in [0,1]")
        if (args.dry_run or args.verify_output) and args.report_json:
            raise ValueError("Read-only modes print JSON to stdout; omit --report-json")
        if args.verify_output and (args.dry_run or args.output_a3m):
            raise ValueError("--verify-output cannot be combined with --dry-run or --output-a3m")
        if args.chain_a3m:
            jobs, report = separate_jobs(args)
            report_path = args.report_json or str(Path(args.output_dir or ".") / "masking.audit.json")
        else:
            if args.output_dir:
                raise ValueError("--output-dir is for separate --chain-a3m inputs")
            if not args.input_a3m:
                raise ValueError("--input-a3m is required")
            doc = read_a3m(args.input_a3m)
            layout = validate_document(doc)
            check_mode(layout, mode)
            print_layout(layout, doc.metadata,
                         sys.stderr if args.dry_run or args.verify_output else sys.stdout)
            if args.query_fasta:
                if not args.chain_sep:
                    raise ValueError("--chain-sep must not be empty")
                check_fasta(layout, read_fasta_chains(args.query_fasta, args.chain_sep))
            full, stochastic = targets(args, layout)
            if args.verify_output:
                report = verify_masking(doc, read_a3m(args.verify_output), full,
                                        stochastic, args.channel_mask_percent)
                report.update({"input_sha256": sha256(args.input_a3m),
                               "output_sha256": sha256(args.verify_output)})
                print(json.dumps(report, indent=2))
                return
            if not args.dry_run and not args.output_a3m:
                raise ValueError("--output-a3m is required unless --dry-run is used")
            job = prepare_job(args, args.input_a3m, full, stochastic, args.output_a3m, doc)
            jobs, report = [job], job["report"]
            report_path = args.report_json or str(args.output_a3m) + ".audit.json"
        if args.dry_run:
            report["dry_run"] = True
            print(json.dumps(report, indent=2))
            return
        publish(jobs, report_path, report)
        for job in jobs:
            r = job["report"]
            print(f"Verified {r['kind']}: {job['output']} | {r['records']} rows | "
                  f"{r['newly_masked_residues']} new X | 0 off-target changes")
        print(f"Audit: {report_path}")
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
