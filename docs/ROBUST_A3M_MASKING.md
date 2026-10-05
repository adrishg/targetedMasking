# Verified A3M masking for LocalColabFold

These standalone commands preserve ColabFold's metadata and use aligned match
columns, not raw character indices. Python 3.10+ and the standard library are
sufficient. Run them from the repository root. Keep the scripts directory
together: the entry points import a shared engine and CLI.

Use `scripts/targetedMasking.py` for every layout:

```bash
python scripts/targetedMasking.py \
  --input-a3m original.a3m --output-a3m masked.a3m --mask A:20-25
```

It detects monomer, homomer or heteromer from stored query blocks and header
copy counts. The console and JSON audit identify the detected type. Compact shared
blocks still require `--all-copies`. It also accepts the separate-subunit
`--chain-a3m` options described below. A headerless A3M is treated as a monomer;
its sequence alone cannot reveal missing multimer metadata or copy counts.

The notebook's `make_masked_a3m` is a self-contained copy of the same rule, so it
runs in Colab without this repository. It masks one combined ColabFold A3M,
preserves the `#lengths<TAB>copy_counts` line, writes one line per sequence,
accepts only `X`, and re-reads the output to reject off-target or missed masks.
Its `mask_chain` names a unique chain (stored block), not a physical copy.

## Why the original multimer utility needed more changes

The previous CLI discarded the initial `#lengths<TAB>copy_counts` line, required
a colon-separated multimer FASTA, and assumed every physical chain had a stored
alignment block. The earlier match-index fix addressed lowercase insertions,
but not compact homomers or the metadata needed by ColabFold.

The updated CLI retains that header, derives stored block boundaries from it,
and independently verifies every resulting alignment before publishing it.
`targetedMasking.py` replaces the original `targetedMasking_multimer.py` and accepts
its options (`--multimer-fasta`, `--mask-chain`, `--mutant-ranges`, `--channel-masking`).

## Formats and coordinates

| Representation | Header example | Stored width | Meaning |
| --- | --- | --- | --- |
| Monomer | absent, or `#1572<TAB>1` | 1572 | One query sequence |
| Compact homotetramer | `#600<TAB>4` | 600 | One block shared by four copies |
| Heterodimer | `#600,150<TAB>1,1` | 750 | Two different stored blocks |
| Mixed stoichiometry | `#600,150<TAB>2,3` | 750 | Block A used twice; block B three times |

`<TAB>` represents a real tab, not a literal backslash followed by `t`.
All positions are **1-based within the stored query block**. `A`, `B`, etc.
refer to stored blocks in header order, not PDB chain IDs, physical copy IDs,
UniProt positions or full-length positions removed from a trimmed construct.

Uppercase residues and `-` each occupy a match column. Lowercase insertions do
not. Each homolog is mapped independently. Masking replaces uppercase residues
with `X`; gaps and insertions remain exactly where they were. No realignment is
performed. Only `X` is supported, to preserve gap topology and block occupancy.

For a compact block with more than one copy, `--all-copies` is required.
Masking that block affects every copy using its MSA. A single physical copy
cannot be selected from a compact shared block. This tool does not silently
expand or re-pair the complex to attempt copy-specific masking.

Combined heteromer A3Ms contain full-width paired rows and full-width unpaired
rows padded with `-` in the other blocks. Their first record must be the complete
ungapped concatenated query with one TAB-separated query ID per stored block.
The first query and exact, identified per-block query rows are protected. Other
homologs, including identical sequences with different IDs, remain eligible.
Duplicate monomer/homomer homolog rows after the first query remain eligible.

## Monomer: Nav1.5 or Cav1.2

```bash
python scripts/targetedMasking.py \
  --input-a3m original.a3m --output-a3m masked.a3m \
  --mask A:660
```

For an in-memory check that writes nothing, replace `--output-a3m masked.a3m`
with `--dry-run`. Examples here demonstrate syntax, not recommended biological
mask definitions. Recover the intended positions for your exact construct.

## Homomer: compact Kv2.1 tetramer

```bash
python scripts/targetedMasking.py \
  --input-a3m kv21_wt.a3m --output-a3m kv21_wt.masked.a3m \
  --mask A:288-328,370-384,401-417 --all-copies
```

The `#600<TAB>4` header remains unchanged.

## Heteromer: one combined ColabFold A3M

```bash
python scripts/targetedMasking.py \
  --input-a3m complex.a3m --output-a3m complex.masked.a3m \
  --mask A:20-25 --mask B:10,15
```

Each block has its own local numbering, even when lengths differ. Paired and
unpaired homologs are masked at the same aligned positions. All row identifiers,
row order, gaps in absent subunits and query records are preserved. If the header
declares repeated copies of a targeted block, add `--all-copies`.

## Heteromer: separate subunit A3Ms

```bash
mkdir -p masked_subunits
python scripts/targetedMasking.py \
  --chain-a3m A=subunit_a.a3m --chain-a3m B=subunit_b.a3m \
  --output-dir masked_subunits --mask A:20-25 --mask B:10
```

Each input must be a monomer A3M with its own query and homologs; the files may
have different lengths and numbers of sequences. Outputs are `A.masked.a3m`,
`B.masked.a3m`, and `masking.audit.json`. Subunits without requested masks are
preserved. Labels are explicit, so file order never establishes pairing.

**These outputs remain separate.** Passing their directory to `colabfold_batch`
does not turn them into one heteromer prediction. The tool does not fabricate
pairing or assemble a combined complex from independent homolog lists. For an
existing complex prediction, prefer masking its combined ColabFold A3M, which
already records the paired/unpaired relationships.

## Stochastic masking and existing options

```bash
python scripts/targetedMasking.py \
  --input-a3m complex.a3m --output-a3m complex.masked.a3m \
  --mask A:25 --stochastic-mask B:20-60 \
  --channel-mask-percent 0.35 --seed 7
```

The fraction is a Bernoulli probability per eligible homolog residue, not an
exact quota. Deterministic masks take precedence over stochastic masks at an
overlap. Zero applies no stochastic masks; one requires every eligible residue
in that region to be masked. Seeded runs are reproducible for the same inputs
and implementation; historical RNG draw sequences are not preserved.

The legacy CLI range spellings still work: `--mask-chain A --mutant-ranges 25`
and `--channel-masking 20-60`. `--mask-chain` now explicitly means a stored block;
`--mask-unit` is its clearer alias. Out-of-range positions cause an error.

An optional `--query-fasta` (alias `--multimer-fasta`) verifies actual sequences,
order and copy counts, rather than lengths alone. Use colon-separated chains
or one FASTA record per chain. Include all physical copies in stored-block
order: an A2B3 example requires `A:A:B:B:B`. Headers remain authoritative.
Multimer files missing their metadata are rejected rather than guessed from
total sequence length. A headerless monomer is accepted by the monomer command.

## What verification guarantees

The verifier independently removes lowercase insertions from original/output
rows and compares aligned positions. It does not reuse the masker's index map.
It checks:

- Same metadata, query sequences, record count, identifiers and record order.
- Correct match-state width in every row, before and after masking.
- All lowercase insertions and gaps unchanged at their original raw positions.
- Every substitution is an uppercase residue becoming `X` at an allowed column.
- Full coverage of deterministic targets and stochastic targets at fraction 1,
  except gaps and protected queries.

Every normal run writes `OUTPUT.a3m.audit.json` with per-site counts of residues,
gaps, pre-existing `X`, newly masked residues and unmasked residues; block widths
and copy counts; protected query rows; seed; and input/output SHA-256 hashes.
Stochastic audits verify coordinates and report coverage, not the RNG history.

Outputs are staged, re-read and checked before publication. Existing outputs,
inputs and audit files are never overwritten. For separate-file batches all
files are prepared before publication, with cleanup of this invocation's new
outputs on a caught failure; this is not a machine-crash-safe multi-file
transaction. Output directories must already exist.

Audit an existing output without writing anything:

```bash
python scripts/targetedMasking.py \
  --input-a3m original.a3m --verify-output candidate.masked.a3m \
  --mask A:660
```

This audit checks substitutions made relative to the supplied original. It
cannot recover amino acids from existing `X` or determine whether an original
was already incorrectly masked. Always regenerate from the unmasked MSA.

Malformed metadata, inconsistent row widths, lowercase/gapped first queries,
ambiguous combined query IDs, mixed-width per-chain records in one file, A2M
dots, internal sequence whitespace and unknown symbols are rejected with errors.
These are explicit supported-format boundaries, not automatic format repairs.

## Verification and authoritative format references

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest
```

An optional check executes the serializer/parser functions extracted from
trusted official ColabFold sources, with template generation stubbed out:

```bash
python -m tests.check_colabfold_compatibility /path/to/colabfold/input.py /path/to/colabfold/batch.py
# ColabFold 1.5.x keeps the serializer in batch.py, so one file is enough:
python -m tests.check_colabfold_compatibility /path/to/colabfold/batch.py
```

It exercises monomer, compact homotetramer, heterodimer and A2B3 serialization,
masking and decoding, checking query sequences and copy counts. No model
inference or GPU libraries are involved. Source hashes are printed for provenance.

- [LocalColabFold documentation](https://github.com/YoshitakaMo/localcolabfold):
  LocalColabFold installs the ColabFold pipeline.
- [ColabFold input serializer](https://github.com/sokrypton/ColabFold/blob/main/colabfold/input.py):
  `msa_to_str`, `pair_sequences`, `pad_sequences` define metadata, stored
  cardinalities, paired rows and gap-padded unpaired blocks.
- [ColabFold reader](https://github.com/sokrypton/ColabFold/blob/main/colabfold/batch.py):
  `unserialize_msa` (and `normalize_a3m` in newer releases) interpret that representation.
- [AlphaFold A3M parser](https://github.com/google-deepmind/alphafold/blob/main/alphafold/data/parsers.py):
  lowercase insertions are removed from aligned rows but contribute insertion
  counts to model features; they should not be discarded from the saved A3M.
