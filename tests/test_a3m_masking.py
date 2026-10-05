import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from scripts.a3m_masking import (
    Document, check_fasta, mask_document, parse_a3m_text, protected_rows,
    read_fasta_chains, resolve_targets, validate_document, verify_masking,
)
from scripts.a3m_masking import read_a3m
from scripts.masking_cli import publish


ROOT = Path(__file__).resolve().parents[1]
MONOMER = ">query\nACDEFG\n>hit\nACggDEFG\n>gap\nAC-EFG\n"
# Same layout as ColabFold input.msa_to_str(pair_sequences + pad_sequences).
# Query, paired homolog, block A query/homolog, block B query/homolog.
HETEROMER = (
    "#4,3\t1,1\n>101\t102\nACDEFGH\n>pairA\tpairB\nACggDEtFGH\n"
    ">101\nACDE---\n>a_hit\nACggDE---\n>102\n----FGH\n>b_hit\n----FtG-\n"
)
HOMOMER = "#6\t4\n" + MONOMER


class CoreTests(unittest.TestCase):
    def test_shifted_raw_index_is_rejected_even_with_identical_lengths(self):
        original = parse_a3m_text(MONOMER)
        wrong = parse_a3m_text(MONOMER.replace("ACggDEFG", "ACggXEFG"))
        with self.assertRaisesRegex(ValueError, "off-target.*column 3"):
            verify_masking(original, wrong, {4})
        correct, report = mask_document(original, {4})
        self.assertEqual(correct.records[1][1], "ACggDEXG")
        self.assertEqual(report["outside_target_changes"], 0)

    def test_missing_full_mask_is_rejected(self):
        doc = parse_a3m_text(MONOMER)
        with self.assertRaisesRegex(ValueError, "missed mask"):
            verify_masking(doc, doc, {4})

    def test_gap_and_insertion_positions_are_preserved(self):
        doc = parse_a3m_text(MONOMER)
        masked, report = mask_document(doc, {2})
        self.assertEqual(masked.records[1][1], "ACggXEFG")
        self.assertEqual(masked.records[2], doc.records[2])
        self.assertEqual(report["coverage"][0]["gaps"], 1)
        wrong = Document(None, (doc.records[0], (">hit", "AgCgDEFG"), doc.records[2]))
        with self.assertRaisesRegex(ValueError, "insertion"):
            verify_masking(doc, wrong, set())

    def test_heteromer_boundaries_paired_unpaired_and_queries(self):
        doc = parse_a3m_text(HETEROMER)
        layout = validate_document(doc)
        targets = resolve_targets(layout, ["B:1-2"])
        masked, report = mask_document(doc, targets)
        self.assertEqual(masked.metadata, "#4,3\t1,1")
        self.assertEqual(masked.records[1][1], "ACggDEtXXH")
        self.assertEqual(masked.records[3], doc.records[3])
        self.assertEqual(masked.records[5][1], "----XtX-")
        for row in (0, 2, 4):
            self.assertEqual(masked.records[row], doc.records[row])
        self.assertEqual(report["protected_query_rows_1based"], [1, 3, 5])
        self.assertEqual([s["position"] for s in report["coverage"]], [1, 2])

    def test_homomer_width_is_not_multiplied_by_cardinality(self):
        doc = parse_a3m_text(HOMOMER)
        layout = validate_document(doc)
        self.assertEqual(layout.width, 6)
        self.assertEqual(layout.kind, "homomer")
        with self.assertRaisesRegex(ValueError, "all-copies"):
            resolve_targets(layout, ["A:5"])
        masked, _ = mask_document(doc, resolve_targets(layout, ["A:5"], True))
        self.assertEqual(masked.metadata, doc.metadata)
        self.assertEqual(masked.records[1][1], "ACggDEXG")

    def test_mixed_cardinality_uses_stored_block_offsets(self):
        doc = parse_a3m_text(HETEROMER.replace("\t1,1", "\t2,3"))
        layout = validate_document(doc)
        self.assertEqual(resolve_targets(layout, ["B:1"], True), {4})
        check_fasta(layout, ["ACDE", "ACDE", "FGH", "FGH", "FGH"])
        with self.assertRaisesRegex(ValueError, "FASTA"):
            check_fasta(layout, ["ACDE", "FGH"])

    def test_same_length_wrong_fasta_is_rejected(self):
        layout = validate_document(parse_a3m_text(MONOMER))
        with self.assertRaisesRegex(ValueError, "FASTA"):
            check_fasta(layout, ["ACDEFA"])

    def test_stochastic_is_reproducible_and_full_targets_take_priority(self):
        doc = parse_a3m_text(">q\nACDEFG\n" + "".join(f">h{i}\nACggDEFG\n" for i in range(30)))
        first, report = mask_document(doc, {0}, range(6), 0.35, 17)
        self.assertEqual(first, mask_document(doc, {0}, range(6), 0.35, 17)[0])
        self.assertNotEqual(first, mask_document(doc, {0}, range(6), 0.35, 18)[0])
        self.assertTrue(all(s[0] == "X" for _, s in first.records[1:]))
        self.assertGreater(sum(c["unmasked_residues"] for c in report["coverage"]), 0)
        self.assertEqual(mask_document(doc, set(), range(6), 0)[0], doc)

    def test_preexisting_x_does_not_count_as_new_mask(self):
        doc = parse_a3m_text(">q\nAC\n>hit\nAX\n")
        _, report = mask_document(doc, {1})
        self.assertEqual(report["newly_masked_residues"], 0)
        self.assertEqual(report["coverage"][0]["preexisting_X"], 1)

    def test_fraction_zero_verifier_rejects_new_random_masks(self):
        doc = parse_a3m_text(MONOMER)
        changed, _ = mask_document(doc, {4})
        with self.assertRaisesRegex(ValueError, "off-target"):
            verify_masking(doc, changed, set(), {4}, 0)

    def test_verifier_rejects_header_order_metadata_and_query_changes(self):
        doc = parse_a3m_text(HOMOMER)
        bad = [Document("#6\t3", doc.records),
               Document(doc.metadata, (doc.records[0], *reversed(doc.records[1:]))),
               Document(doc.metadata, ((">query", "XCDEFG"), *doc.records[1:])),
               Document(doc.metadata, doc.records[:-1])]
        for output in bad:
            with self.subTest(output=output), self.assertRaises(ValueError):
                verify_masking(doc, output, {0})

    def test_parser_rejects_ambiguous_or_invalid_formats(self):
        bad = ["#6\\t4\n" + MONOMER, "#6,3\t1\n" + MONOMER,
               "#7\t1\n" + MONOMER, ">q\nAC\n>h\nA\n",
               ">q\nAcC\n>h\nACC\n", ">q\nA-C\n", ">q\nAC\n>h\nA.C\n",
               ">q\nAC\n>h\nA C\n", ">q\nAC\n>empty\n", "AC\n>q\nAC\n",
               ">q\nAC\n#2\t1\n>h\nAC\n"]
        for text in bad:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_a3m_text(text)

    def test_invalid_targets_and_fractions_fail(self):
        doc = parse_a3m_text(MONOMER)
        layout = validate_document(doc)
        for spec in ("A:0", "A:7", "B:1", "A:-1", "A:1,", "A:"):
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                resolve_targets(layout, [spec])
        for fraction in (-1, 2, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                mask_document(doc, {1}, fraction=fraction)

    def test_wrapped_sequences_and_legacy_io_preserve_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "in.a3m", Path(tmp) / "out.a3m"
            source.write_text(HOMOMER.replace("ACggDEFG", "ACgg\nDEFG"))
            output.write_text(read_a3m(source).text())
            self.assertEqual(output.read_text(), HOMOMER)

    def test_multirecord_and_colon_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chains.fasta"
            for content in (">complex\nACDE:FGH\n", ">A\nACDE\n>B\nFGH\n"):
                path.write_text(content)
                self.assertEqual(read_fasta_chains(path), ["ACDE", "FGH"])


class CliTests(unittest.TestCase):
    def run_cli(self, mode, *args):
        assert mode == "auto"
        return subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / "targetedMasking.py"),
                               *map(str, args)], capture_output=True, text=True)

    def test_auto_detects_headerless_and_headered_monomer_homomer_and_heteromer(self):
        examples = [(MONOMER, "monomer", []), ("#6\t1\n" + MONOMER, "monomer", []),
                    (HOMOMER, "homomer", ["--all-copies"]), (HETEROMER, "heteromer", []),
                    (HETEROMER.replace("\t1,1", "\t2,3"), "heteromer", ["--all-copies"])]
        for content, kind, extra in examples:
            with self.subTest(kind=kind, extra=extra), tempfile.TemporaryDirectory() as tmp:
                source = Path(tmp) / "input.a3m"
                source.write_text(content)
                result = self.run_cli("auto", "--input-a3m", source, "--mask", "A:2", "--dry-run", *extra)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(result.stdout)["kind"], kind)
                self.assertIn(f"Identified {kind}", result.stderr)
                self.assertIn("Block A:", result.stderr)

    def test_layout_summary_distinguishes_blocks_and_physical_copies(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "input.a3m", Path(tmp) / "output.a3m"
            source.write_text(HETEROMER.replace("\t1,1", "\t2,3"))
            result = self.run_cli("auto", "--input-a3m", source, "--output-a3m", output,
                                  "--mask", "B:2", "--all-copies")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Identified heteromer: 5 chain(s), 2 stored alignment block(s)", result.stdout)
            self.assertIn("Block B: 3 residues, 3 copy/copies, aligned columns 5-7", result.stdout)
            self.assertIn("masking it affects all 3 copies", result.stdout)

    def test_auto_does_not_silently_mask_one_copy_of_a_compact_homomer(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.a3m"
            source.write_text(HOMOMER)
            result = self.run_cli("auto", "--input-a3m", source, "--mask", "A:2", "--dry-run")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--all-copies", result.stderr)

    def test_combined_entry_points_write_verified_audits(self):
        for mode, text, extra in (("auto", MONOMER, []),
                                  ("auto", HOMOMER, ["--all-copies"]),
                                  ("auto", HETEROMER, [])):
            with self.subTest(mode=mode, text=text[:6]), tempfile.TemporaryDirectory() as tmp:
                source, output = Path(tmp) / "in.a3m", Path(tmp) / "out.a3m"
                source.write_text(text)
                result = self.run_cli(mode, "--input-a3m", source, "--output-a3m", output,
                                      "--mask", "A:2", *extra)
                self.assertEqual(result.returncode, 0, result.stderr)
                report = json.loads(Path(str(output) + ".audit.json").read_text())
                self.assertTrue(report["serialized_output_verified"])
                self.assertEqual(report["outside_target_changes"], 0)
                self.assertEqual(source.read_text(), text)
                verify = self.run_cli(mode, "--input-a3m", source, "--verify-output", output,
                                      "--mask", "A:2", *extra)
                self.assertEqual(verify.returncode, 0, verify.stderr)

    def test_no_files_written_on_invalid_input_or_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / "in.a3m", Path(tmp) / "out.a3m"
            source.write_text(MONOMER)
            for extra in (["--mask", "A:7"], ["--mask", "A:2", "--report-json", source]):
                result = self.run_cli("auto", "--input-a3m", source, "--output-a3m", output, *extra)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(output.exists())
                self.assertEqual(source.read_text(), MONOMER)

    def test_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "in.a3m"
            source.write_text(MONOMER)
            args = ("--input-a3m", source, "--mask", "A:5", "--dry-run")
            result = self.run_cli("auto", *args)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(json.loads(result.stdout)["dry_run"])
            self.assertEqual(list(Path(tmp).iterdir()), [source])

    def test_separate_subunits_keep_independent_queries_and_do_not_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.a3m", Path(tmp) / "b.a3m"
            a.write_text(MONOMER)
            b.write_text(">queryB\nFGH\n>hitB\nFggGH\n")
            result = self.run_cli("auto", "--chain-a3m", f"A={a}", "--chain-a3m", f"B={b}",
                                  "--output-dir", tmp, "--mask", "B:2")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((Path(tmp) / "A.masked.a3m").read_text(), MONOMER)
            self.assertEqual((Path(tmp) / "B.masked.a3m").read_text(), ">queryB\nFGH\n>hitB\nFggXH\n")
            audit = json.loads((Path(tmp) / "masking.audit.json").read_text())
            self.assertEqual(len(audit["files"]), 2)
            self.assertTrue(all(f["serialized_output_verified"] for f in audit["files"]))

    def test_separate_subunit_failure_leaves_no_partial_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.a3m", Path(tmp) / "b.a3m"
            a.write_text(MONOMER)
            b.write_text(">q\nFGH\n>bad\nFG\n")
            result = self.run_cli("auto", "--chain-a3m", f"A={a}", "--chain-a3m", f"B={b}",
                                  "--output-dir", tmp, "--mask", "A:2")
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(set(Path(tmp).iterdir()), {a, b})

    def test_publication_rechecks_serialized_columns_before_saving(self):
        doc = parse_a3m_text(MONOMER)
        wrong = parse_a3m_text(MONOMER.replace("ACggDEFG", "ACggXEFG"))
        with tempfile.TemporaryDirectory() as tmp:
            report = {}
            job = dict(original=doc, masked=wrong, full={4}, stochastic=set(),
                       fraction=1, output=Path(tmp) / "out.a3m", report=report)
            with self.assertRaisesRegex(ValueError, "off-target"):
                publish([job], Path(tmp) / "audit.json", report)
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_publication_rolls_back_its_outputs_on_report_failure(self):
        import os
        doc = parse_a3m_text(MONOMER)
        masked, report = mask_document(doc, {4})
        real_link = os.link

        def fail_report(source, destination):
            if str(destination).endswith("audit.json"):
                raise OSError("Simulated publication failure")
            real_link(source, destination)

        with tempfile.TemporaryDirectory() as tmp:
            job = dict(original=doc, masked=masked, full={4}, stochastic=set(),
                       fraction=1, output=Path(tmp) / "out.a3m", report=report)
            with patch("scripts.masking_cli.os.link", side_effect=fail_report):
                with self.assertRaisesRegex(OSError, "Simulated"):
                    publish([job], Path(tmp) / "audit.json", report)
            self.assertEqual(list(Path(tmp).iterdir()), [])


class MatchStateRegressionTests(unittest.TestCase):
    """Cases carried over from the former targetedMasking_multimer tests."""

    def mask(self, text, columns):
        doc = parse_a3m_text(text)
        masked, _ = mask_document(doc, columns)
        self.assertEqual(masked.records[0], doc.records[0])
        return [s for _, s in masked.records[1:]]

    def test_lowercase_insertions_do_not_advance_match_state_numbering(self):
        self.assertEqual(self.mask(">q\nABCDE\n>hit\naAbbcBCdDE\n", {2}), ["aAbbcBXdDE"])

    def test_deletion_gap_remains_gap_at_masked_match_state(self):
        self.assertEqual(self.mask(">q\nABCDE\n>hit\nAB-DE\n", {2}), ["AB-DE"])

    def test_multiple_insertions_before_and_adjacent_to_target_are_preserved(self):
        self.assertEqual(self.mask(">q\nABCDE\n>hit\nAabcBdeCDE\n", {1, 2}), ["AabcXdeXDE"])

    def test_deterministic_and_stochastic_targets_together(self):
        doc = parse_a3m_text(">q\nABCDE\n>hit1\nABCDE\n>hit2\nABcCDE\n")
        masked, _ = mask_document(doc, {1}, {0, 2}, fraction=1.0)
        self.assertEqual([s for _, s in masked.records[1:]], ["XXXDE", "XXcXDE"])

    def test_heteromer_boundary_masks_only_selected_block(self):
        doc = parse_a3m_text("#3,2\t1,1\n>101\t102\nABCDE\n>hit\nABcCDE\n")
        columns = resolve_targets(validate_document(doc), ["B:1-2"])
        masked, _ = mask_document(doc, columns)
        self.assertEqual(masked.records[1][1], "ABcCXX")

    def test_inconsistent_homolog_match_state_count_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "observed 4 A3M match states"):
            parse_a3m_text(">q\nABCDE\n>bad_hit\nABcDE\n")

    def test_positions_beyond_block_are_rejected(self):
        doc = parse_a3m_text("#3,2\t1,1\n>101\t102\nABCDE\n>hit\nABCDE\n")
        with self.assertRaisesRegex(ValueError, "exceeds block length 3"):
            resolve_targets(validate_document(doc), ["A:4"])


if __name__ == "__main__":
    unittest.main()
