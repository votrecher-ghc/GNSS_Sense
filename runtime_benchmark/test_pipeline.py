"""Correctness checks for the isolated StarDial reference pipeline."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import stardial_pipeline as stardial


class StarDialPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.pipeline = stardial.StarDialPipeline()
        cls.request, cls.truth = stardial.generate_synthetic_input()

    def test_template_library_matches_repository_order(self) -> None:
        self.assertEqual(len(stardial.TEMPLATE_ORDER), 12)
        for label in stardial.TEMPLATE_ORDER:
            trace = stardial.template_trace(label, 160)
            self.assertEqual(trace.shape, (160, 2))
            self.assertTrue(np.all(np.isfinite(trace)))

    def test_batched_dtw_is_exact(self) -> None:
        rng = np.random.default_rng(42)
        first = rng.normal(size=(24, 2))
        templates = rng.normal(size=(4, 24, 2))
        expected = np.array(
            [stardial._dtw_distance(first, template) for template in templates]
        )
        actual = stardial._batched_dtw_distance(first, templates)
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-12)

    def test_genuine_accept_and_wrong_claim_reject(self) -> None:
        genuine = self.pipeline.authenticate(self.request)
        self.assertTrue(genuine.accepted)
        self.assertEqual(genuine.predicted_template, "Star")
        self.assertTrue(genuine.gesture_pass)
        self.assertTrue(genuine.physical_pass)
        self.assertEqual(genuine.trajectory_xy_m.shape, (100, 2))

        wrong_claim = stardial.AuthenticationInput(
            cn0_dbhz=self.request.cn0_dbhz,
            los_contact_xy_m=self.request.los_contact_xy_m,
            claimed_template="A",
        )
        rejected = self.pipeline.authenticate(wrong_claim)
        self.assertFalse(rejected.accepted)
        self.assertFalse(rejected.gesture_pass)
        self.assertTrue(rejected.physical_pass)

    def test_fixed_seed_is_deterministic(self) -> None:
        second, second_truth = stardial.generate_synthetic_input()
        np.testing.assert_array_equal(second.cn0_dbhz, self.request.cn0_dbhz)
        np.testing.assert_array_equal(second.los_contact_xy_m, self.request.los_contact_xy_m)
        np.testing.assert_array_equal(second_truth, self.truth)

    def test_satellite_geometry_mismatch_fails_physical_check(self) -> None:
        # Preserve the waveform shape but break its satellite-to-LoS mapping,
        # analogous to replaying evidence under a different constellation.
        permutation = np.random.default_rng(1).permutation(self.request.cn0_dbhz.shape[1])
        mismatched = stardial.AuthenticationInput(
            cn0_dbhz=self.request.cn0_dbhz[:, permutation],
            los_contact_xy_m=self.request.los_contact_xy_m,
            claimed_template="Star",
        )
        result = self.pipeline.authenticate(mismatched)
        self.assertEqual(result.predicted_template, "Star")
        self.assertTrue(result.gesture_pass)
        self.assertFalse(result.physical_pass)
        self.assertFalse(result.accepted)


if __name__ == "__main__":
    unittest.main()
