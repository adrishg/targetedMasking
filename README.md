# Targeted Masking

Targeted Masking is a framework built on top of AlphaFold2/ColabFold designed to reduce excessive conformational bias introduced by evolutionary information in multiple sequence alignments (MSAs). The goal is to improve exploration of alternative structural states, particularly for dynamic systems such as voltage-gated ion channels and membrane proteins.

In standard AlphaFold/ColabFold workflows, the MSA strongly constrains residue-residue relationships toward the dominant evolutionary state represented in sequence databases. While this is extremely powerful for structure prediction, it can limit sampling of biologically relevant alternative conformations, especially for proteins that undergo large-scale functional transitions.

Targeted Masking introduces controlled perturbations into the MSA by selectively masking:

- mutation sites,
- flexible loops,
- gating charge regions,
- experimentally relevant dynamic motifs,
- or user-defined residue ranges.

The framework currently supports:

- deterministic masking of specific residues ("mutant masking"),
- stochastic masking of selected regions ("channel masking" / conformational masking),
- multimer-aware masking,
- and customizable masking percentages to tune conformational variability.

The original motivation for this work came from studies of voltage-gated calcium ion channels (VGICs), particularly CaV1.2 voltage-sensing domains, where experimentally relevant conformational heterogeneity is difficult to capture using conventional AlphaFold pipelines alone. However, the approach is designed to be generalizable to many dynamic protein systems.

This repository aims to provide:

- reproducible masking utilities,
- ColabFold-compatible workflows,
- exploratory conformational sampling tools,
- and a lightweight framework for hypothesis-driven structural perturbation experiments.

The philosophy behind Targeted Masking is not to "force" a structure into a desired state, but rather to partially relax evolutionary constraints in carefully selected regions and allow the model to explore alternative energetic and conformational solutions.

## Repository Contents

- [targetedMasking_AF2_v1_Cav12_vsd2Test.ipynb](targetedMasking_AF2_v1_Cav12_vsd2Test.ipynb): ColabFold/AlphaFold2 notebook with targeted MSA masking controls and visualization.
- [scripts/targetedMasking_multimer.py](scripts/targetedMasking_multimer.py): standalone A3M masking utility for multimer-aware, query-referenced residue masking.
- [tests/test_targeted_masking_multimer.py](tests/test_targeted_masking_multimer.py): regression tests for the core masking behavior.

## Quick Start

### Notebook Workflow

Open the notebook in Google Colab, select a GPU runtime, and run the setup/input cells as usual. The notebook includes an optional targeted MSA masking section before prediction. Use that section to define mutation sites, dynamic regions, masking percentages, and chain-aware masking behavior.

### Command-Line A3M Masking

The CLI masks selected chain-local query positions in a multimer A3M while preserving the first sequence, which is treated as the query.

```bash
python scripts/targetedMasking_multimer.py \
  --input-a3m input.a3m \
  --output-a3m masked.a3m \
  --multimer-fasta multimer.fasta \
  --mask-chain A \
  --mutant-ranges 10,25,42-48 \
  --channel-masking 90-120 \
  --channel-mask-percent 0.35 \
  --seed 7
```

Input assumptions:

- multimer FASTA chains are separated with `:` by default, for example `CHAIN_A:CHAIN_B`;
- residue positions are 1-based and chain-local;
- lowercase A3M insertions and gaps do not count toward query residue numbering;
- the first A3M sequence is preserved to maintain the original query.

## Development

The standalone masking utility currently uses only the Python standard library. To run the regression tests:

```bash
python -m unittest
```

## Citation / Related Work

If you use this repository, please cite and/or refer to:

- Hernandez-Gonzalez et al., unpublished work on conformational state sampling in voltage-gated calcium channels.
- Hernandez-Gonzalez et al., unpublished VGIC mutants study.
- Related CaV1.2 conformational modeling repositories and associated future publications from the Yarov-Yarovoy Lab, UC Davis.

Additional citations and preprints will be added as manuscripts become publicly available.

## Disclaimer

This repository is intended primarily as a research and hypothesis-generation framework. Targeted masking modifies the information content of the MSA and therefore may produce structures that differ from the dominant evolutionary state. Interpretation of resulting models should be guided by experimental evidence, biophysical reasoning, and appropriate validation strategies.
