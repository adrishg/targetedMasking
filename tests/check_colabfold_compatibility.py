"""Optional compatibility check against downloaded official ColabFold source.

Usage: python -m tests.check_colabfold_compatibility [INPUT_PY] BATCH_PY
ColabFold 1.5.x defines the serializer in batch.py, so BATCH_PY alone is enough;
newer releases move it to input.py and add normalize_a3m, which is loaded if present.
Only the serializer/parser function definitions are loaded; no GPU libraries or
model inference are needed. Supply trusted official source files.
"""
import ast
import hashlib
from pathlib import Path
import sys
import typing

from scripts.a3m_masking import mask_document, parse_a3m_text, resolve_targets, validate_document


REQUIRED = {"pair_sequences", "pad_sequences", "pair_msa", "msa_to_str", "unserialize_msa"}
OPTIONAL = {"normalize_a3m"}


def check(*paths):
    namespace = {**vars(typing), "mk_mock_template": lambda seq: {"query": seq}}
    found = set()
    for path in map(Path, paths):
        tree = ast.parse(path.read_text())
        nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name in (REQUIRED | OPTIONAL) - found]
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
        found.update(n.name for n in nodes)
        print(path.name, "sha256", hashlib.sha256(path.read_bytes()).hexdigest())
    if REQUIRED - found:
        raise ValueError(f"Serializer/parser functions not found: {sorted(REQUIRED - found)}")
    for queries, copies in [(["ACDEFG"], [1]), (["ACDEFG"], [4]),
                            (["ACDE", "FGH"], [1, 1]), (["ACDE", "FGH"], [2, 3])]:
        unpaired = [f">{101+i}\n{q}\n>hit{i}\n{q[:1]}gg{q[1:]}\n" for i, q in enumerate(queries)]
        paired = ([f">{101+i}\n{q}\n>pair{i}\n{q[:2]}t{q[2:]}\n" for i, q in enumerate(queries)]
                  if len(queries) > 1 else None)
        doc = parse_a3m_text(namespace["msa_to_str"](unpaired, paired, queries, copies))
        layout = validate_document(doc)
        result, report = mask_document(doc, resolve_targets(layout, ["A:2"], True))
        decoded = namespace["unserialize_msa"]([result.text()], doc.records[0][1])
        assert decoded[2] == queries, (decoded[2], queries)
        assert decoded[3] == copies, (decoded[3], copies)
        for index, msa in enumerate(decoded[0]):
            assert msa.splitlines()[1] == queries[index], msa
        print("PASS official serialize -> mask -> official parse:", layout.kind, copies,
              report["records"], "rows")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    check(*sys.argv[1:])
