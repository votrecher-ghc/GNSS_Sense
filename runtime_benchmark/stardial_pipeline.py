"""Minimal, reproducible StarDial online authentication pipeline.

This module is intentionally isolated from the plotting and offline evaluation
code in the MATLAB repository.  It implements the online path described in
Section III of the paper:

    C/N0 cleaning -> diffraction evidence -> LoS geometric inversion
    -> gesture-template score -> physical consistency -> joint decision

The implementation consumes already parsed, time-aligned arrays.  RINEX I/O,
ephemeris parsing, plotting, threshold calibration, and result serialization
are deliberately outside the measured inference boundary.
"""

from __future__ import annotations

import os

# Fix numerical-library threading before NumPy is imported.  The benchmark is
# intended to describe one CPU process and must not change thread count between
# runs or machines without recording that fact.
for _thread_var in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

from dataclasses import dataclass
from time import perf_counter_ns
from typing import Any

import numpy as np


TEMPLATE_ORDER = (
    "A",
    "C",
    "M",
    "Star",
    "L",
    "X",
    "Z",
    "N",
    "V",
    "Rectangle",
    "LeftSwipe",
    "RightSwipe",
)


# Polyline definitions follow gesture_analysis/utils/gesture_template_library.m.
_TEMPLATE_VERTICES: dict[str, np.ndarray] = {
    "A": np.array(
        [[-0.40, -0.50], [0.00, 0.50], [0.40, -0.50], [-0.20, -0.10], [0.20, -0.10]],
        dtype=np.float64,
    ),
    "C": np.array(
        [[0.35, 0.40], [-0.35, 0.40], [-0.35, -0.40], [0.35, -0.40]],
        dtype=np.float64,
    ),
    "M": np.array(
        [[-0.40, -0.40], [-0.40, 0.40], [0.00, 0.00], [0.40, 0.40], [0.40, -0.40]],
        dtype=np.float64,
    ),
    "Star": np.array(
        [
            [-0.30, -0.45],
            [0.00, 0.55],
            [0.30, -0.45],
            [-0.48, 0.15],
            [0.48, 0.15],
            [-0.30, -0.45],
        ],
        dtype=np.float64,
    ),
    "L": np.array([[-0.30, 0.40], [-0.30, -0.40], [0.30, -0.40]], dtype=np.float64),
    "X": np.array(
        [[-0.30, 0.40], [0.30, -0.40], [0.30, 0.40], [-0.30, -0.40]], dtype=np.float64
    ),
    "Z": np.array(
        [[-0.30, 0.40], [0.30, 0.40], [-0.30, -0.40], [0.70, -0.40]], dtype=np.float64
    ),
    "N": np.array(
        [[-0.30, -0.40], [-0.30, 0.40], [0.30, -0.40], [0.30, 0.40]], dtype=np.float64
    ),
    "V": np.array([[-0.35, 0.40], [0.00, -0.40], [0.35, 0.40]], dtype=np.float64),
    "Rectangle": np.array(
        [[-0.35, 0.30], [0.35, 0.30], [0.35, -0.30], [-0.35, -0.30], [-0.35, 0.30]],
        dtype=np.float64,
    ),
    "LeftSwipe": np.array([[0.38, 0.00], [-0.38, 0.00]], dtype=np.float64),
    "RightSwipe": np.array([[-0.38, 0.00], [0.38, 0.00]], dtype=np.float64),
}


@dataclass(frozen=True)
class PipelineConfig:
    """Online parameters; offline-calibrated values are fixed for benchmarking."""

    sample_rate_hz: float = 25.0
    compare_points: int = 160
    grid_points_per_axis: int = 27
    grid_min_m: float = -0.52
    grid_max_m: float = 0.52
    body_anchor_x_m: float = 0.0
    body_anchor_y_m: float = -0.68
    fresnel_sigma_m: float = 0.050
    beam_width: int = 48
    max_jump_m: float = 0.18
    smoothness_weight: float = 16.0
    emission_weight: float = 18.0
    dtw_weight: float = 0.35
    rmse_weight: float = 0.45
    shape_weight: float = 0.20
    softmax_temperature: float = 0.16
    gesture_score_threshold: float = 0.12
    physical_score_threshold: float = 0.58
    event_threshold: float = 0.48


@dataclass(frozen=True)
class AuthenticationInput:
    """Time-aligned online input after RINEX and ephemeris parsing."""

    cn0_dbhz: np.ndarray
    los_contact_xy_m: np.ndarray
    claimed_template: str
    sample_rate_hz: float = 25.0

    def validate(self) -> None:
        if self.cn0_dbhz.ndim != 2:
            raise ValueError("cn0_dbhz must have shape [frames, satellites]")
        if self.los_contact_xy_m.shape != (*self.cn0_dbhz.shape, 2):
            raise ValueError("los_contact_xy_m must have shape [frames, satellites, 2]")
        if self.cn0_dbhz.shape[0] < 8 or self.cn0_dbhz.shape[1] < 4:
            raise ValueError("at least 8 frames and 4 visible satellites are required")
        if self.claimed_template not in TEMPLATE_ORDER:
            raise ValueError(f"unsupported claimed template: {self.claimed_template}")
        if not np.all(np.isfinite(self.los_contact_xy_m)):
            raise ValueError("LoS contact points must be finite")


@dataclass(frozen=True)
class AuthenticationResult:
    accepted: bool
    predicted_template: str
    claimed_template: str
    claimed_score: float
    best_score: float
    physical_score: float
    gesture_pass: bool
    physical_pass: bool
    trajectory_xy_m: np.ndarray
    template_scores: dict[str, float]
    diagnostics: dict[str, float]
    stage_ms: dict[str, float]


def resample_polyline(points: np.ndarray, count: int) -> np.ndarray:
    """Arc-length resampling equivalent to the MATLAB template utilities."""

    points = np.asarray(points, dtype=np.float64)
    finite = np.all(np.isfinite(points), axis=1)
    points = points[finite]
    if len(points) == 0:
        return np.full((count, 2), np.nan, dtype=np.float64)
    if len(points) == 1:
        return np.repeat(points, count, axis=0)

    segment = np.linalg.norm(np.diff(points, axis=0), axis=1)
    keep = np.r_[True, segment > 1e-12]
    points = points[keep]
    if len(points) == 1:
        return np.repeat(points, count, axis=0)
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    target = np.linspace(0.0, distance[-1], count)
    return np.column_stack(
        [np.interp(target, distance, points[:, dim]) for dim in range(2)]
    )


def template_trace(label: str, count: int) -> np.ndarray:
    if label not in _TEMPLATE_VERTICES:
        raise ValueError(f"unsupported template: {label}")
    return resample_polyline(_TEMPLATE_VERTICES[label], count)


def _normalize_xy(points: np.ndarray) -> np.ndarray:
    centered = points - np.mean(points, axis=0, keepdims=True)
    span = float(np.max(np.ptp(centered, axis=0)))
    return centered / max(span, 1e-9)


def _point_to_segments_distance(
    points: np.ndarray, segment_end: np.ndarray, segment_start: np.ndarray
) -> np.ndarray:
    """Distance from S points to M body-to-hand candidate segments -> [M, S]."""

    ab = segment_end - segment_start[None, :]
    ap = points - segment_start[None, :]
    length2 = np.sum(ab * ab, axis=1) + 1e-12
    fraction = (ab @ ap.T) / length2[:, None]
    fraction = np.clip(fraction, 0.0, 1.0)
    closest = segment_start[None, None, :] + fraction[:, :, None] * ab[:, None, :]
    return np.linalg.norm(points[None, :, :] - closest, axis=2)


def _median3(values: np.ndarray) -> np.ndarray:
    padded = np.pad(values, ((1, 1), (0, 0)), mode="edge")
    return np.median(np.stack((padded[:-2], padded[1:-1], padded[2:])), axis=0)


def _moving_mean3(values: np.ndarray) -> np.ndarray:
    padded = np.pad(values, ((1, 1), (0, 0)), mode="edge")
    return (padded[:-2] + padded[1:-1] + padded[2:]) / 3.0


def _trace_features(trace_xy: np.ndarray) -> dict[str, float | bool]:
    delta = np.diff(trace_xy, axis=0)
    segment_len = np.linalg.norm(delta, axis=1)
    span_xy = np.ptp(trace_xy, axis=0)
    scale = max(float(np.max(span_xy)), 1e-12)
    start, stop = trace_xy[0], trace_xy[-1]

    turn_count = 0
    for idx in range(1, len(trace_xy) - 1):
        v1 = trace_xy[idx] - trace_xy[idx - 1]
        v2 = trace_xy[idx + 1] - trace_xy[idx]
        n1, n2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
        if n1 < 0.01 or n2 < 0.01:
            continue
        cosine = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
        if np.degrees(np.arccos(cosine)) >= 48.0:
            turn_count += 1

    start_sign = 1 if start[0] > 1e-6 else (-1 if start[0] < -1e-6 else 0)
    stop_sign = 1 if stop[0] > 1e-6 else (-1 if stop[0] < -1e-6 else 0)
    return {
        "same_x_sign": start_sign * stop_sign > 0,
        "x_progress": float((stop[0] - start[0]) / scale),
        "y_progress": float((stop[1] - start[1]) / scale),
        "horizontal_ratio": float(span_xy[1] / max(span_xy[0], 1e-12)),
        "path_ratio": float(np.sum(segment_len) / scale),
        "corner_count": float(turn_count),
        "end_gap": float(np.linalg.norm(stop - start) / scale),
        "endpoint_x_min": float(min(start[0], stop[0])),
        "endpoint_x_max": float(max(start[0], stop[0])),
    }


def _shape_penalty(label: str, feature: dict[str, float | bool]) -> float:
    f = feature
    penalty = 0.0
    if label == "A":
        penalty += 0.35 if f["corner_count"] < 2 else 0.0
        penalty += 0.30 if f["x_progress"] < 0.08 else 0.0
        penalty += 0.25 if f["path_ratio"] < 1.7 else 0.0
    elif label == "RightSwipe":
        penalty += 0.90 if f["x_progress"] < 0.22 else 0.0
        penalty += 0.60 if f["horizontal_ratio"] > 0.28 else 0.0
        penalty += 0.25 if f["corner_count"] > 1 else 0.0
    elif label == "LeftSwipe":
        penalty += 0.90 if f["x_progress"] > -0.22 else 0.0
        penalty += 0.60 if f["horizontal_ratio"] > 0.28 else 0.0
    elif label == "C":
        penalty += 0.80 if not f["same_x_sign"] else 0.0
        penalty += 0.45 if f["endpoint_x_min"] < 0.01 else 0.0
        penalty += 0.25 if abs(f["x_progress"]) > 0.20 else 0.0
    elif label == "L":
        penalty += 0.45 if f["y_progress"] > -0.18 else 0.0
        penalty += 0.30 if f["x_progress"] < 0.12 else 0.0
    elif label == "V":
        penalty += 0.45 if f["y_progress"] < 0.18 else 0.0
        penalty += 0.20 if f["corner_count"] < 1 else 0.0
    elif label == "X":
        penalty += 0.55 if f["corner_count"] < 2 else 0.0
        penalty += 0.35 if f["end_gap"] < 0.20 else 0.0
    elif label == "Star":
        penalty += 0.55 if f["corner_count"] < 4 else 0.0
        penalty += 0.30 if f["path_ratio"] < 2.4 else 0.0
    elif label == "Rectangle":
        penalty += 0.55 if f["corner_count"] < 3 else 0.0
        penalty += 0.60 if f["end_gap"] > 0.24 else 0.0
    elif label == "M":
        penalty += 0.25 if f["corner_count"] < 3 else 0.0
    elif label in {"N", "Z"}:
        penalty += 0.20 if f["corner_count"] < 2 else 0.0
    return penalty


def _dtw_distance(first: np.ndarray, second: np.ndarray) -> float:
    """Exact dynamic-programming DTW, normalized as in the MATLAB module."""

    cost = np.linalg.norm(first[:, None, :] - second[None, :, :], axis=2)
    n_rows, n_cols = cost.shape
    previous = np.full(n_cols + 1, np.inf, dtype=np.float64)
    previous[0] = 0.0
    for row_idx in range(n_rows):
        current = np.full(n_cols + 1, np.inf, dtype=np.float64)
        row = cost[row_idx]
        for col_idx in range(1, n_cols + 1):
            current[col_idx] = row[col_idx - 1] + min(
                previous[col_idx], current[col_idx - 1], previous[col_idx - 1]
            )
        previous = current
    return float(previous[-1] / max(n_rows + n_cols, 1))


def _batched_dtw_distance(first: np.ndarray, templates: np.ndarray) -> np.ndarray:
    """Exact DTW for all templates, evaluated on independent anti-diagonals."""

    cost = np.linalg.norm(first[None, :, None, :] - templates[:, None, :, :], axis=3)
    template_count, n_rows, n_cols = cost.shape
    cumulative = np.full(
        (template_count, n_rows + 1, n_cols + 1), np.inf, dtype=np.float64
    )
    cumulative[:, 0, 0] = 0.0
    template_index = np.arange(template_count)[:, None]
    for diagonal in range(2, n_rows + n_cols + 1):
        row_min = max(1, diagonal - n_cols)
        row_max = min(n_rows, diagonal - 1)
        rows = np.arange(row_min, row_max + 1)
        cols = diagonal - rows
        predecessor = np.minimum(
            np.minimum(
                cumulative[template_index, rows[None, :] - 1, cols[None, :]],
                cumulative[template_index, rows[None, :], cols[None, :] - 1],
            ),
            cumulative[template_index, rows[None, :] - 1, cols[None, :] - 1],
        )
        cumulative[template_index, rows[None, :], cols[None, :]] = (
            cost[template_index, rows[None, :] - 1, cols[None, :] - 1] + predecessor
        )
    return cumulative[:, -1, -1] / max(n_rows + n_cols, 1)


class StarDialPipeline:
    """Reusable online authenticator with templates and search grid preloaded."""

    def __init__(self, config: PipelineConfig | None = None) -> None:
        self.config = config or PipelineConfig()
        axis = np.linspace(
            self.config.grid_min_m,
            self.config.grid_max_m,
            self.config.grid_points_per_axis,
            dtype=np.float64,
        )
        gx, gy = np.meshgrid(axis, axis, indexing="xy")
        self.grid_xy = np.column_stack((gx.ravel(), gy.ravel()))
        self.body_anchor = np.array(
            [self.config.body_anchor_x_m, self.config.body_anchor_y_m], dtype=np.float64
        )
        transition_delta = self.grid_xy[:, None, :] - self.grid_xy[None, :, :]
        transition_distance2 = np.sum(transition_delta * transition_delta, axis=2)
        self.transition_score = -self.config.smoothness_weight * transition_distance2
        self.transition_score[transition_distance2 > self.config.max_jump_m**2] = -np.inf
        self.templates = {
            label: _normalize_xy(template_trace(label, self.config.compare_points))
            for label in TEMPLATE_ORDER
        }

    def _preprocess_cn0(self, cn0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(cn0, dtype=np.float64)
        if np.any(~np.isfinite(values)):
            # Linear interpolation is kept inside the online stage because raw
            # receiver epochs can contain isolated missing C/N0 samples.
            values = values.copy()
            frame = np.arange(values.shape[0])
            for sat_idx in range(values.shape[1]):
                finite = np.isfinite(values[:, sat_idx])
                if np.count_nonzero(finite) < 2:
                    values[:, sat_idx] = np.nanmedian(values)
                else:
                    values[:, sat_idx] = np.interp(
                        frame, frame[finite], values[finite, sat_idx]
                    )
        cleaned = _median3(values)
        baseline = np.percentile(cleaned, 88.0, axis=0, keepdims=True)
        departure = np.clip(baseline - cleaned, 0.0, None)
        scale = np.percentile(departure, 95.0, axis=0, keepdims=True)
        scale = np.maximum(scale, 1.0)
        normalized_departure = np.clip(departure / scale, 0.0, 1.0)
        return cleaned, normalized_departure

    def _extract_diffraction(self, departure: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        second_difference = np.zeros_like(departure)
        second_difference[1:-1] = np.abs(
            departure[2:] - 2.0 * departure[1:-1] + departure[:-2]
        )
        oscillation_scale = np.percentile(second_difference, 95.0, axis=0, keepdims=True)
        oscillation = np.clip(second_difference / np.maximum(oscillation_scale, 0.12), 0.0, 1.0)
        continuity = _moving_mean3((departure >= 0.28).astype(np.float64))
        saliency = np.clip(0.60 * departure + 0.20 * oscillation + 0.20 * continuity, 0.0, 1.0)
        levels = np.digitize(saliency, bins=(0.28, 0.58)).astype(np.uint8)
        return saliency, levels

    def _recover_trajectory(self, evidence: np.ndarray, contact_xy: np.ndarray) -> np.ndarray:
        cfg = self.config
        n_frames = evidence.shape[0]
        n_states = len(self.grid_xy)
        beam = min(cfg.beam_width, n_states)

        # Evaluate all frame/state/satellite likelihoods in one array.  This is
        # algebraically identical to the frame loop in the MATLAB reference but
        # avoids Python dispatch overhead in the measured implementation.
        ab = self.grid_xy - self.body_anchor[None, :]
        ap = contact_xy - self.body_anchor[None, None, :]
        length2 = np.sum(ab * ab, axis=1) + 1e-12
        fraction = np.einsum("md,tsd->tms", ab, ap, optimize=True) / length2[None, :, None]
        fraction = np.clip(fraction, 0.0, 1.0)
        closest = (
            self.body_anchor[None, None, None, :]
            + fraction[:, :, :, None] * ab[None, :, None, :]
        )
        distance2 = np.sum((contact_xy[:, None, :, :] - closest) ** 2, axis=3)
        predicted = np.exp(-0.5 * distance2 / (cfg.fresnel_sigma_m**2))
        sat_weight = 0.20 + 0.80 * evidence
        sat_weight /= np.sum(sat_weight, axis=1, keepdims=True) + 1e-12
        mismatch = np.sum(
            sat_weight[:, None, :] * (predicted - evidence[:, None, :]) ** 2, axis=2
        )
        state_radius = np.linalg.norm(self.grid_xy - self.body_anchor[None, :], axis=1)
        radius_penalty = np.maximum(state_radius - 1.20, 0.0) ** 2
        emission = -cfg.emission_weight * mismatch - 0.05 * radius_penalty[None, :]

        first_states = np.argpartition(emission[0], -beam)[-beam:]
        order = np.argsort(emission[0, first_states])[::-1]
        state_history = np.empty((n_frames, beam), dtype=np.int32)
        parent_history = np.full((n_frames, beam), -1, dtype=np.int32)
        state_history[0] = first_states[order]
        previous_score = emission[0, state_history[0]].copy()

        for frame_idx in range(1, n_frames):
            transition = (
                previous_score[:, None]
                + self.transition_score[state_history[frame_idx - 1]]
            )
            best_parent = np.argmax(transition, axis=0)
            score = emission[frame_idx] + transition[best_parent, np.arange(n_states)]
            chosen = np.argpartition(score, -beam)[-beam:]
            chosen = chosen[np.argsort(score[chosen])[::-1]]
            state_history[frame_idx] = chosen
            parent_history[frame_idx] = best_parent[chosen]
            previous_score = score[chosen]

        path_index = np.empty(n_frames, dtype=np.int32)
        beam_position = 0
        for frame_idx in range(n_frames - 1, -1, -1):
            path_index[frame_idx] = state_history[frame_idx, beam_position]
            if frame_idx > 0:
                beam_position = parent_history[frame_idx, beam_position]
        trajectory = self.grid_xy[path_index]

        # A short edge-preserving moving average matches the final smoothing
        # used by the MATLAB trajectory implementation.
        padded = np.pad(trajectory, ((2, 2), (0, 0)), mode="edge")
        smoothed = sum(padded[offset : offset + n_frames] for offset in range(5)) / 5.0
        return smoothed

    def _match_templates(
        self, trajectory: np.ndarray
    ) -> tuple[str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        cfg = self.config
        trace = _normalize_xy(resample_polyline(trajectory, cfg.compare_points))
        feature = _trace_features(trace)
        template_array = np.stack([self.templates[label] for label in TEMPLATE_ORDER])
        error = np.linalg.norm(trace[None, :, :] - template_array, axis=2)
        rmse = np.sqrt(np.mean(error * error, axis=1))
        dtw = _batched_dtw_distance(trace, template_array)
        shape = np.empty_like(rmse)
        for idx, label in enumerate(TEMPLATE_ORDER):
            shape[idx] = _shape_penalty(label, feature)
        distance = cfg.dtw_weight * dtw + cfg.rmse_weight * rmse + cfg.shape_weight * shape
        logits = -distance / max(cfg.softmax_temperature, 1e-9)
        weights = np.exp(logits - np.max(logits))
        scores = weights / np.sum(weights)
        predicted = TEMPLATE_ORDER[int(np.argmax(scores))]
        return predicted, scores, distance, rmse, dtw

    def _physical_consistency(
        self, trajectory: np.ndarray, evidence: np.ndarray, contact_xy: np.ndarray
    ) -> tuple[float, dict[str, float]]:
        cfg = self.config
        ab = trajectory - self.body_anchor[None, :]
        ap = contact_xy - self.body_anchor[None, None, :]
        length2 = np.sum(ab * ab, axis=1) + 1e-12
        fraction = np.einsum("td,tsd->ts", ab, ap, optimize=True) / length2[:, None]
        fraction = np.clip(fraction, 0.0, 1.0)
        closest = self.body_anchor[None, None, :] + fraction[:, :, None] * ab[:, None, :]
        distance2 = np.sum((contact_xy - closest) ** 2, axis=2)
        predicted = np.exp(-0.5 * distance2 / (cfg.fresnel_sigma_m**2))

        active = evidence >= cfg.event_threshold
        silent = evidence <= 0.18
        event_score = float(np.mean(predicted[active])) if np.any(active) else 0.0
        exclusion_score = float(np.mean(1.0 - predicted[silent])) if np.any(silent) else 0.0
        agreement_score = float(np.exp(-np.mean((predicted - evidence) ** 2) / 0.16))

        timing_scores: list[float] = []
        duration = max(evidence.shape[0] - 1, 1)
        for sat_idx in range(evidence.shape[1]):
            if np.max(evidence[:, sat_idx]) < cfg.event_threshold:
                continue
            observed_peak = int(np.argmax(evidence[:, sat_idx]))
            predicted_peak = int(np.argmax(predicted[:, sat_idx]))
            timing_scores.append(np.exp(-abs(observed_peak - predicted_peak) / (0.10 * duration)))
        temporal_score = float(np.mean(timing_scores)) if timing_scores else 0.0

        physical_score = (
            0.35 * event_score
            + 0.25 * exclusion_score
            + 0.25 * agreement_score
            + 0.15 * temporal_score
        )
        detail = {
            "event_consistency": event_score,
            "exclusion_consistency": exclusion_score,
            "waveform_agreement": agreement_score,
            "temporal_consistency": temporal_score,
            "active_event_fraction": float(np.mean(active)),
        }
        return float(physical_score), detail

    def authenticate(self, request: AuthenticationInput) -> AuthenticationResult:
        """Run one complete online authentication inference."""

        request.validate()
        stage_ms: dict[str, float] = {}
        total_start = perf_counter_ns()

        start = perf_counter_ns()
        _, departure = self._preprocess_cn0(request.cn0_dbhz)
        stage_ms["preprocessing"] = (perf_counter_ns() - start) / 1e6

        start = perf_counter_ns()
        evidence, levels = self._extract_diffraction(departure)
        stage_ms["diffraction_features"] = (perf_counter_ns() - start) / 1e6

        start = perf_counter_ns()
        trajectory = self._recover_trajectory(evidence, request.los_contact_xy_m)
        stage_ms["geometric_inversion"] = (perf_counter_ns() - start) / 1e6

        start = perf_counter_ns()
        predicted, scores, distance, rmse, dtw = self._match_templates(trajectory)
        stage_ms["template_matching"] = (perf_counter_ns() - start) / 1e6

        start = perf_counter_ns()
        physical_score, physical_detail = self._physical_consistency(
            trajectory, evidence, request.los_contact_xy_m
        )
        stage_ms["physical_validation"] = (perf_counter_ns() - start) / 1e6

        start = perf_counter_ns()
        claim_idx = TEMPLATE_ORDER.index(request.claimed_template)
        claimed_score = float(scores[claim_idx])
        best_score = float(np.max(scores))
        gesture_pass = predicted == request.claimed_template and claimed_score >= self.config.gesture_score_threshold
        physical_pass = physical_score >= self.config.physical_score_threshold
        accepted = bool(gesture_pass and physical_pass)
        stage_ms["joint_decision"] = (perf_counter_ns() - start) / 1e6
        stage_ms["total_internal"] = (perf_counter_ns() - total_start) / 1e6

        diagnostics: dict[str, float] = {
            **physical_detail,
            "salient_segment_fraction": float(np.mean(levels > 0)),
            "best_distance": float(np.min(distance)),
            "best_rmse": float(rmse[np.argmax(scores)]),
            "best_dtw": float(dtw[np.argmax(scores)]),
        }
        return AuthenticationResult(
            accepted=accepted,
            predicted_template=predicted,
            claimed_template=request.claimed_template,
            claimed_score=claimed_score,
            best_score=best_score,
            physical_score=physical_score,
            gesture_pass=bool(gesture_pass),
            physical_pass=bool(physical_pass),
            trajectory_xy_m=trajectory,
            template_scores={label: float(scores[idx]) for idx, label in enumerate(TEMPLATE_ORDER)},
            diagnostics=diagnostics,
            stage_ms=stage_ms,
        )


def generate_synthetic_input(
    *,
    label: str = "Star",
    claimed_template: str | None = None,
    frames: int = 100,
    satellites: int = 16,
    sample_rate_hz: float = 25.0,
    seed: int = 20260811,
    config: PipelineConfig | None = None,
) -> tuple[AuthenticationInput, np.ndarray]:
    """Create deterministic, physically coupled C/N0 and LoS benchmark input.

    Satellite contact points remain nearly fixed during the short window, as
    real GNSS geometry does.  Their placement is selected from different rays
    of the gesture fan; diffraction strength follows the distance from each LoS
    contact point to the body-to-hand segment.  Ground truth is returned only
    for smoke tests and is never passed to ``authenticate``.
    """

    if label not in TEMPLATE_ORDER:
        raise ValueError(f"unsupported label: {label}")
    if satellites < 4:
        raise ValueError("satellites must be at least 4")
    cfg = config or PipelineConfig()
    rng = np.random.default_rng(seed)
    body = np.array([cfg.body_anchor_x_m, cfg.body_anchor_y_m], dtype=np.float64)
    truth = template_trace(label, frames)
    truth -= np.mean(truth, axis=0, keepdims=True)
    truth *= 0.92

    # Select contact points along multiple body-to-gesture rays.  This produces
    # an instantaneous LoS configuration that constrains both direction and
    # radial extent without moving satellites with the hand.
    phase_index = np.linspace(0, frames - 1, satellites, endpoint=False)
    phase_index = np.rint(phase_index).astype(int)
    phase_index = np.clip(phase_index + rng.integers(-2, 3, satellites), 0, frames - 1)
    radial_fraction = rng.uniform(0.30, 0.92, satellites)
    base_contact = body + radial_fraction[:, None] * (truth[phase_index] - body)
    base_contact += rng.normal(0.0, 0.008, size=base_contact.shape)

    # Satellite motion over four seconds is small but nonzero.  The drift keeps
    # the physical check explicitly tied to the current time-indexed geometry.
    drift_velocity = rng.normal(0.0, 0.00035, size=(satellites, 2))
    centered_time = (np.arange(frames) - (frames - 1) / 2.0) / sample_rate_hz
    contacts = base_contact[None, :, :] + centered_time[:, None, None] * drift_velocity[None, :, :]

    occultation = np.empty((frames, satellites), dtype=np.float64)
    for frame_idx in range(frames):
        distance = _point_to_segments_distance(
            contacts[frame_idx], truth[frame_idx][None, :], body
        )[0]
        occultation[frame_idx] = np.exp(-0.5 * (distance / cfg.fresnel_sigma_m) ** 2)

    baseline = rng.uniform(42.0, 47.0, size=(1, satellites))
    phase = rng.uniform(0.0, 2.0 * np.pi, size=(1, satellites))
    time_axis = np.arange(frames, dtype=np.float64)[:, None] / sample_rate_hz
    diffraction = 0.55 * occultation * np.sin(2.0 * np.pi * 3.2 * time_axis + phase)
    thermal_noise = rng.normal(0.0, 0.16, size=(frames, satellites))
    cn0 = baseline - 8.0 * occultation + diffraction + thermal_noise
    cn0 = np.rint(cn0)  # u-blox-style quantized C/N0 reports

    request = AuthenticationInput(
        cn0_dbhz=cn0,
        los_contact_xy_m=contacts,
        claimed_template=claimed_template or label,
        sample_rate_hz=sample_rate_hz,
    )
    return request, truth


def result_to_dict(result: AuthenticationResult) -> dict[str, Any]:
    """JSON-friendly compact result (trajectory intentionally omitted)."""

    return {
        "accepted": result.accepted,
        "predicted_template": result.predicted_template,
        "claimed_template": result.claimed_template,
        "claimed_score": result.claimed_score,
        "best_score": result.best_score,
        "physical_score": result.physical_score,
        "gesture_pass": result.gesture_pass,
        "physical_pass": result.physical_pass,
        "diagnostics": result.diagnostics,
        "stage_ms": result.stage_ms,
    }
