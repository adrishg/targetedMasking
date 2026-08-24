# A3M match-state indexing correction

The earlier multimer utility derived raw character indices from the query and
applied those indices directly to every homolog row. In A3M format, lowercase
insertions are row-specific and do not occupy match states, so a homolog with
one or more insertions before a target could be masked at the wrong residue.

The corrected implementation maps each homolog row independently from A3M
match-state positions to raw-string indices. Uppercase residues and deletion
gaps occupy match states; lowercase insertions do not advance numbering; gaps
remain gaps; and the query row remains unchanged. Rows whose match-state count
does not match the multimer FASTA are rejected rather than partially masked.

Regression tests cover lowercase insertions, multiple and adjacent insertions,
deletion gaps, query preservation, match-state length validation, and multimer
chain boundaries.

The correction demonstrates a defect in the prior implementation. It does not
by itself establish that historical prediction inputs were affected. That
requires the exact production A3Ms and masking command records, which are not
available in the companion repository audit.
