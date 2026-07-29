# Targeted Masking

[![Open the Targeted Masking notebook in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/targetedMasking/blob/main/targetedMasking_AF2_v1_Cav12_vsd2Test.ipynb)

<p align="right"><img src="./assets/targeted-masking-conformational-sampling_bgless.png" alt="Targeted masking and conformational sampling" width="300"></p>

Targeted Masking is an AlphaFold2/ColabFold workflow for hypothesis-driven conformational sampling.
It selectively masks user-defined columns of a multiple sequence alignment (MSA) while retaining the
query sequence and the evolutionary information outside the selected region.

The goal is not to force AlphaFold2 toward a chosen structure. The goal is to relax local evolutionary
constraints at mutation sites or conformationally variable regions while preserving the broader MSA
context that helps AlphaFold2 model large proteins accurately.

## Rationale

AlphaFold2 has learned general principles of protein structure, but predictions of large multidomain
proteins still benefit strongly from the residue-coupling information encoded in an MSA. For
voltage-gated ion channels, an intact MSA helps maintain the transmembrane fold, pore architecture,
domain packing, and long-range relationships across a very large sequence.

That evolutionary information can also favor the dominant sequence-compatible geometry. At mutation
sites or regions known or suspected to participate in alternative conformations, a strong local MSA
signal may reduce structural variability and bias repeated predictions toward one basin. Removing too
much MSA information is not an ideal solution because it can compromise the global accuracy gained
from homologous sequences.

Targeted masking is a middle ground:

- keep the query sequence unchanged;
- retain the MSA outside the selected region;
- mask selected alignment columns in homologous sequences;
- generate vanilla and targeted-masked ensembles under otherwise matched settings;
- test whether variability increases locally without disrupting the rest of the protein.

Candidate targets include mutation sites, voltage-sensor segments, intracellular gates, flexible
linkers, ligand-coupling motifs, and other regions supported by a structural or functional hypothesis.

This is a sampling strategy, not a state classifier. Increased variability does not by itself prove
that additional models represent functional conformations. Predictions require convergence checks,
structural-quality control, comparison across independent seeds, and validation against experimental
structures.

## Workflow

```text
sequence or multimer FASTA
          │
          ▼
ColabFold MSA generation
          │
          ├──────────────► vanilla AlphaFold2 prediction
          │
          ▼
select mutation sites or dynamic regions
          │
          ▼
mask those MSA columns in homologous sequences
keep the query and all other MSA columns intact
          │
          ▼
targeted-masked AlphaFold2 prediction
          │
          ▼
compare ensembles using structural and experimental coordinates
```

The companion [VGIC mutant ensemble analysis](https://github.com/adrishg/vgci_mutants) repository
implements the downstream evaluation for Cav1.2, Nav1.5, and Kv2.1. It compares vanilla and
targeted-masked ensembles using convergence-filtered distance distributions, mutation-site analyses,
structural integrity checks, and experimental PDB references.

## Repository contents

- [targetedMasking_AF2_v1_Cav12_vsd2Test.ipynb](targetedMasking_AF2_v1_Cav12_vsd2Test.ipynb):
  ColabFold/AlphaFold2 notebook that runs matched vanilla and targeted-masked predictions.
- [scripts/targetedMasking_multimer.py](scripts/targetedMasking_multimer.py):
  standalone A3M masking utility with multimer-aware, query-referenced residue numbering.
- [tests/test_targeted_masking_multimer.py](tests/test_targeted_masking_multimer.py):
  regression tests for range parsing, chain handling, query preservation, and masking behavior.

## Quick start

### Colab notebook

Open the notebook in Google Colab and select a GPU runtime. The default workflow:

1. accepts the input sequence or complex;
2. generates the original ColabFold MSA;
3. runs the vanilla prediction;
4. creates a targeted-masked A3M from the same MSA;
5. runs the targeted-masked prediction using matched model settings;
6. reports the applied positions and displays the resulting structures.

Define mutation sites and dynamic regions using 1-based, chain-local residue positions. Run the
vanilla and targeted-masked branches with the same seeds, model configuration, recycle settings, and
other sampling parameters whenever a controlled comparison is required.

### Command-line A3M masking

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
- the first A3M sequence is preserved as the query;
- targeted positions are masked in homologous rows, not deleted from the alignment.

## Choosing a masking experiment

A useful masking experiment should begin with a defined question:

- Which site or region is expected to be conformationally variable?
- What structural coordinate will report that variability?
- Which experimental structures define relevant reference geometries?
- What unmasked regions must remain stable for the model to be considered valid?

The masked region should be no broader than the hypothesis requires. Compare multiple mask definitions
when possible, and retain a vanilla control generated from the same starting MSA.

## Validation

At minimum, evaluate:

1. convergence across recycles and seeds;
2. global fold and domain integrity;
3. local distributions at the masked region;
4. control regions outside the mask;
5. overlap with coordinate-matched experimental structures;
6. whether apparent broadening represents reproducible clusters rather than isolated failures.

The [vgci_mutants](https://github.com/adrishg/vgci_mutants) analysis provides concrete examples of this
validation strategy, including extra pore-integrity checks for tetrameric Kv2.1 and IFM-latch
coordinates for Nav1.5.

## Development

The standalone masking utility uses the Python standard library. Run its regression tests with:

```bash
python -m unittest
```

## Citation and related work

If you use this repository, please cite or refer to:

- Hernandez-Gonzalez et al., unpublished work on conformational sampling in voltage-gated ion channels.
- Hernandez-Gonzalez et al., unpublished VGIC mutant ensemble study.
- Associated future publications from the Yarov-Yarovoy Lab, UC Davis.

Additional citations and preprints will be added as manuscripts become publicly available.

## Disclaimer

Targeted masking changes the information supplied to AlphaFold2. Resulting models are hypotheses that
must be evaluated with experimental evidence, biophysical reasoning, and appropriate structural
quality control. Masking should not be interpreted as evidence that a predicted state is populated in
cells or accessible under a particular physiological condition.
