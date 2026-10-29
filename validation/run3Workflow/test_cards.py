"""Preflight regression for Run-3 cards with explicit veto-map settings."""
from pathlib import Path
import tempfile
import unittest

from check_cards import parse


class CardPreflightTest(unittest.TestCase):
    def parse_text(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "card.input"
            path.write_text(text)
            return parse(path)

    def test_veto_configuration_preserves_paths_and_numeric_values(self):
        card = self.parse_text("""DataType 10
JetVetoMapFile jetVetoMaps/OO2025/Summer24Prompt25_RunCDEFG.root # selected map
JetVetoMapHistogram jetvetomap_all
JetVetoPtConvention stored
MixingFileList mixingFileList/pool.txt
JetVetoInputMinPt 15
WeightExponent 1 2
""")
        self.assertEqual(card["JetVetoMapFile"], ["jetVetoMaps/OO2025/Summer24Prompt25_RunCDEFG.root"])
        self.assertEqual(card["JetVetoMapHistogram"], ["jetvetomap_all"])
        self.assertEqual(card["JetVetoPtConvention"], ["stored"])
        self.assertEqual(card["MixingFileList"], ["mixingFileList/pool.txt"])
        self.assertEqual(card["JetVetoInputMinPt"], [15.0])
        self.assertEqual(card["WeightExponent"], [1.0, 2.0])

    def test_invalid_numeric_value_is_rejected(self):
        with self.assertRaises(ValueError):
            self.parse_text("JetVetoInputMinPt stored\n")

    def test_duplicate_string_key_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate key"):
            self.parse_text("JetVetoPtConvention raw\nJetVetoPtConvention stored\n")


if __name__ == "__main__":
    unittest.main()
