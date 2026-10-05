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

Historical inputs were affected. The Kv2.1 masked MSAs (WT and F412L, all six
archived versions) match the raw-string-position rule in every row and the
query-column rule in only 7,337 of 20,478 homolog rows. Both the notebook's
`make_masked_a3m` and the pre-fix multimer script applied query-derived raw
indices to every row, which is that rule. Masked
conditions built with that code must be regenerated and re-predicted.
The notebook now uses the per-row query-column rule and verifies its output.
