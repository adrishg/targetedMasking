import os
import tempfile
import unittest

from scripts.targetedMasking_multimer import (
    infer_chain_lengths_from_fasta,
    match_state_indices,
    mask_records,
    parse_ranges,
)


class TargetedMaskingMultimerTests(unittest.TestCase):
    def test_parse_ranges_accepts_ranges_and_singletons(self):
        self.assertEqual(parse_ranges("5, 2-4, 10"), [2, 3, 4, 5, 10])
        self.assertEqual(parse_ranges("4-2"), [2, 3, 4])

    def test_infer_chain_lengths_from_colon_separated_fasta(self):
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write(">example\n")
            handle.write("AAAA:BBB:CC\n")
            path = handle.name

        try:
            self.assertEqual(infer_chain_lengths_from_fasta(path), [4, 3, 2])
        finally:
            os.unlink(path)

    def test_mask_records_preserves_query_and_masks_selected_chain(self):
        records = [
            (">query", "ABCDE"),
            (">hit1", "ABCDE"),
            (">hit2", "ABcCDE"),
        ]

        masked, applied, missing = mask_records(
            records=records,
            mutant_qpos=[2],
            channel_qpos=[1, 3],
            channel_fraction=1.0,
            masking_char="X",
            chain_lengths=[3, 2],
            mask_chain_index=0,
        )

        self.assertEqual(masked[0], records[0])
        self.assertEqual(masked[1], (">hit1", "XXXDE"))
        self.assertEqual(masked[2], (">hit2", "XXcXDE"))
        self.assertEqual(applied, [1, 2, 3])
        self.assertEqual(missing, [])

    def test_lowercase_insertions_do_not_advance_match_state_numbering(self):
        records = [
            (">query", "ABCDE"),
            (">hit", "aAbbcBCdDE"),
        ]
        masked, _, _ = mask_records(
            records, mutant_qpos=[3], channel_qpos=[], channel_fraction=1.0,
            masking_char="X", chain_lengths=[3, 2], mask_chain_index=0,
        )
        self.assertEqual(masked[1], (">hit", "aAbbcBXdDE"))
        self.assertEqual(match_state_indices(records[1][1]), [1, 5, 6, 8, 9])

    def test_deletion_gap_remains_gap_at_masked_match_state(self):
        records = [(">query", "ABCDE"), (">hit", "AB-DE")]
        masked, _, _ = mask_records(
            records, mutant_qpos=[3], channel_qpos=[], channel_fraction=1.0,
            masking_char="X", chain_lengths=[3, 2], mask_chain_index=0,
        )
        self.assertEqual(masked[1], records[1])

    def test_multiple_insertions_before_and_adjacent_to_target_are_preserved(self):
        records = [(">query", "ABCDE"), (">hit", "AabcBdeCDE")]
        masked, _, _ = mask_records(
            records, mutant_qpos=[2, 3], channel_qpos=[], channel_fraction=1.0,
            masking_char="X", chain_lengths=[3, 2], mask_chain_index=0,
        )
        self.assertEqual(masked[1], (">hit", "AabcXdeXDE"))

    def test_multimer_boundary_masks_only_selected_chain(self):
        records = [(">query", "ABCDE"), (">hit", "ABcCDE")]
        masked, applied, missing = mask_records(
            records, mutant_qpos=[1, 2], channel_qpos=[], channel_fraction=1.0,
            masking_char="X", chain_lengths=[3, 2], mask_chain_index=1,
        )
        self.assertEqual(masked[0], records[0])
        self.assertEqual(masked[1], (">hit", "ABcCXX"))
        self.assertEqual(applied, [1, 2])
        self.assertEqual(missing, [])

    def test_inconsistent_homolog_match_state_count_is_rejected(self):
        records = [(">query", "ABCDE"), (">bad_hit", "ABcDE")]
        with self.assertRaisesRegex(ValueError, "observed 4 A3M match states"):
            mask_records(
                records, mutant_qpos=[2], channel_qpos=[], channel_fraction=1.0,
                masking_char="X", chain_lengths=[3, 2], mask_chain_index=0,
            )

    def test_missing_positions_are_reported(self):
        records = [
            (">query", "ABCDE"),
            (">hit1", "ABCDE"),
        ]

        _masked, applied, missing = mask_records(
            records=records,
            mutant_qpos=[4],
            channel_qpos=[],
            channel_fraction=1.0,
            masking_char="X",
            chain_lengths=[3, 2],
            mask_chain_index=0,
        )

        self.assertEqual(applied, [])
        self.assertEqual(missing, [4])


if __name__ == "__main__":
    unittest.main()
