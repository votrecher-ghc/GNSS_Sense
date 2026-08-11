"""Measure latency, CPU, memory, throughput, cold start, and input scaling.

Run from the repository root:

    python runtime_benchmark/benchmark_runtime.py

Only the online ``StarDialPipeline.authenticate`` call is included in steady
latency.  Input generation, template/grid construction, serialization,
plotting, and offline calibration are excluded and reported separately.
"""

from __future__ import annotations

import os

for _thread_var in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_var] = "1"

import argparse
import csv
import ctypes
import gc
import json
import platform
import subprocess
import sys
import threading
import time
import tracemalloc
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from stardial_pipeline import (
    AuthenticationInput,
    PipelineConfig,
    StarDialPipeline,
    generate_synthetic_input,
    result_to_dict,
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
MIB = 1024.0 * 1024.0


def _percentile_summary(values: list[float] | np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(array)),
        "std": float(np.std(array, ddof=1)) if len(array) > 1 else 0.0,
        "min": float(np.min(array)),
        "median": float(np.median(array)),
        "p95": float(np.percentile(array, 95)),
        "p99": float(np.percentile(array, 99)),
        "max": float(np.max(array)),
    }


def _windows_memory() -> dict[str, int]:
    class ProcessMemoryCountersEx(ctypes.Structure):
        _fields_ = [
            ("cb", ctypes.c_ulong),
            ("PageFaultCount", ctypes.c_ulong),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
            ("PrivateUsage", ctypes.c_size_t),
        ]

    counters = ProcessMemoryCountersEx()
    counters.cb = ctypes.sizeof(counters)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    psapi.GetProcessMemoryInfo.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ProcessMemoryCountersEx),
        ctypes.c_ulong,
    ]
    psapi.GetProcessMemoryInfo.restype = ctypes.c_bool
    handle = kernel32.GetCurrentProcess()
    ok = psapi.GetProcessMemoryInfo(
        handle, ctypes.byref(counters), counters.cb
    )
    if not ok:
        raise ctypes.WinError()
    return {
        "rss_bytes": int(counters.WorkingSetSize),
        "peak_rss_bytes": int(counters.PeakWorkingSetSize),
        "private_bytes": int(counters.PrivateUsage),
    }


def _process_memory() -> dict[str, int]:
    if sys.platform == "win32":
        return _windows_memory()
    try:
        import resource

        peak = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        if sys.platform != "darwin":
            peak *= 1024
        return {"rss_bytes": peak, "peak_rss_bytes": peak, "private_bytes": 0}
    except (ImportError, OSError):
        return {"rss_bytes": 0, "peak_rss_bytes": 0, "private_bytes": 0}


def _total_memory_bytes() -> int:
    if sys.platform != "win32":
        return 0

    class MemoryStatusEx(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    status = MemoryStatusEx()
    status.dwLength = ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        return 0
    return int(status.ullTotalPhys)


def _parse_cpu_affinity(specification: str) -> list[int]:
    if not specification.strip():
        return []
    selected: set[int] = set()
    for item in specification.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            start_text, end_text = item.split("-", 1)
            start, end = int(start_text), int(end_text)
            if end < start:
                raise ValueError(f"invalid CPU affinity range: {item}")
            selected.update(range(start, end + 1))
        else:
            selected.add(int(item))
    logical_cpus = os.cpu_count() or 1
    if not selected or min(selected) < 0 or max(selected) >= logical_cpus:
        raise ValueError(f"CPU affinity must select processors in [0, {logical_cpus - 1}]")
    return sorted(selected)


def _set_cpu_affinity(processors: list[int]) -> None:
    if not processors:
        return
    if sys.platform == "win32":
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        kernel32.SetProcessAffinityMask.restype = ctypes.c_bool
        mask = sum(1 << processor for processor in processors)
        if not kernel32.SetProcessAffinityMask(kernel32.GetCurrentProcess(), mask):
            raise ctypes.WinError(ctypes.get_last_error())
        return
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(processors))


def _cpu_model() -> str:
    if sys.platform == "win32":
        try:
            import winreg

            key_path = r"HARDWARE\DESCRIPTION\System\CentralProcessor\0"
            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, key_path) as key:
                return str(winreg.QueryValueEx(key, "ProcessorNameString")[0]).strip()
        except OSError:
            pass
    return platform.processor() or "unknown"


def _git_commit() -> str:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_DIR,
            text=True,
            capture_output=True,
            timeout=5,
            check=True,
        )
        return completed.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _environment_info(cpu_affinity: list[int]) -> dict[str, Any]:
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "os": platform.platform(),
        "cpu_model": _cpu_model(),
        "logical_cpu_count": os.cpu_count() or 1,
        "cpu_affinity_logical_processors": cpu_affinity or "operating-system default",
        "total_memory_gib": _total_memory_bytes() / (1024.0**3),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "git_commit": _git_commit(),
        "timer": "time.perf_counter_ns",
        "thread_environment": {
            key: os.environ.get(key, "unset")
            for key in (
                "OPENBLAS_NUM_THREADS",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
    }


def _measure_iterations(
    pipeline: StarDialPipeline,
    request: AuthenticationInput,
    iterations: int,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    logical_cpus = os.cpu_count() or 1
    rows: list[dict[str, Any]] = []
    batch_wall_start = time.perf_counter_ns()
    batch_cpu_start = time.process_time_ns()
    for iteration in range(iterations):
        wall_start = time.perf_counter_ns()
        cpu_start = time.process_time_ns()
        result = pipeline.authenticate(request)
        cpu_ms = (time.process_time_ns() - cpu_start) / 1e6
        wall_ms = (time.perf_counter_ns() - wall_start) / 1e6
        row: dict[str, Any] = {
            "iteration": iteration,
            "wall_ms": wall_ms,
            "cpu_ms": cpu_ms,
            "cpu_core_equivalent_pct": 100.0 * cpu_ms / max(wall_ms, 1e-12),
            "cpu_system_normalized_pct": 100.0 * cpu_ms / max(wall_ms * logical_cpus, 1e-12),
            "accepted": int(result.accepted),
            "predicted_template": result.predicted_template,
            "claimed_score": result.claimed_score,
            "physical_score": result.physical_score,
        }
        for stage, value in result.stage_ms.items():
            row[f"stage_{stage}_ms"] = value
        rows.append(row)
    batch_cpu_ms = (time.process_time_ns() - batch_cpu_start) / 1e6
    batch_wall_ms = (time.perf_counter_ns() - batch_wall_start) / 1e6
    batch = {
        "wall_ms": batch_wall_ms,
        "cpu_ms": batch_cpu_ms,
        "cpu_time_per_inference_ms": batch_cpu_ms / iterations,
        "cpu_core_equivalent_pct": 100.0 * batch_cpu_ms / batch_wall_ms,
        "cpu_system_normalized_pct": 100.0 * batch_cpu_ms / (batch_wall_ms * logical_cpus),
        "throughput_inferences_per_second": 1000.0 * iterations / batch_wall_ms,
    }
    return rows, batch


def _measure_memory(
    pipeline: StarDialPipeline, request: AuthenticationInput, iterations: int = 8
) -> dict[str, float]:
    gc.collect()
    before = _process_memory()
    samples = [before["rss_bytes"]]
    stop_event = threading.Event()

    def sample_working_set() -> None:
        while not stop_event.is_set():
            samples.append(_process_memory()["rss_bytes"])
            time.sleep(0.0005)

    sampler = threading.Thread(target=sample_working_set, daemon=True)
    sampler.start()
    try:
        for _ in range(iterations):
            pipeline.authenticate(request)
    finally:
        stop_event.set()
        sampler.join(timeout=2.0)

    # Track Python/NumPy allocations in a separate pass so tracemalloc's own
    # bookkeeping does not inflate the sampled OS working-set peak above.
    gc.collect()
    tracemalloc.start()
    pipeline.authenticate(request)
    _, traced_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    after = _process_memory()
    sampled_peak = max(samples + [after["rss_bytes"]])
    return {
        "baseline_rss_mib": before["rss_bytes"] / MIB,
        "post_measurement_rss_mib": after["rss_bytes"] / MIB,
        "sampled_peak_rss_mib": sampled_peak / MIB,
        "sampled_incremental_rss_mib": max(0, sampled_peak - before["rss_bytes"]) / MIB,
        "process_lifetime_peak_rss_mib": after["peak_rss_bytes"] / MIB,
        "private_memory_mib": after["private_bytes"] / MIB,
        "traced_peak_allocation_mib": traced_peak / MIB,
        "measurement_iterations": float(iterations),
    }


def _cold_worker(args: argparse.Namespace) -> int:
    input_start = time.perf_counter_ns()
    request, _ = generate_synthetic_input(
        label=args.label,
        frames=args.frames,
        satellites=args.satellites,
        seed=args.seed,
    )
    input_ms = (time.perf_counter_ns() - input_start) / 1e6
    init_start = time.perf_counter_ns()
    pipeline = StarDialPipeline()
    init_ms = (time.perf_counter_ns() - init_start) / 1e6
    inference_start = time.perf_counter_ns()
    result = pipeline.authenticate(request)
    inference_ms = (time.perf_counter_ns() - inference_start) / 1e6
    print(
        json.dumps(
            {
                "input_generation_ms": input_ms,
                "pipeline_init_ms": init_ms,
                "first_inference_ms": inference_ms,
                "pipeline_to_decision_ms": init_ms + inference_ms,
                "accepted": result.accepted,
                "predicted_template": result.predicted_template,
            },
            separators=(",", ":"),
        )
    )
    return 0


def _measure_cold_start(args: argparse.Namespace) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--cold-worker",
        "--label",
        args.label,
        "--frames",
        str(args.frames),
        "--satellites",
        str(args.satellites),
        "--seed",
        str(args.seed),
        "--cpu-affinity",
        args.cpu_affinity,
    ]
    child_env = os.environ.copy()
    for key in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        child_env[key] = "1"
    for run_index in range(args.cold_runs):
        start = time.perf_counter_ns()
        completed = subprocess.run(
            command,
            cwd=REPO_DIR,
            env=child_env,
            text=True,
            capture_output=True,
            timeout=120,
            check=True,
        )
        process_ms = (time.perf_counter_ns() - start) / 1e6
        output_lines = [line for line in completed.stdout.splitlines() if line.strip()]
        record = json.loads(output_lines[-1])
        record["run"] = run_index
        record["process_start_to_exit_ms"] = process_ms
        records.append(record)
    summaries = {
        key: _percentile_summary([float(record[key]) for record in records])
        for key in (
            "input_generation_ms",
            "pipeline_init_ms",
            "first_inference_ms",
            "pipeline_to_decision_ms",
            "process_start_to_exit_ms",
        )
    }
    return {"runs": records, "summary": summaries}


def _measure_scaling(
    pipeline: StarDialPipeline, args: argparse.Namespace
) -> list[dict[str, Any]]:
    cases = [
        ("satellites", 8, args.frames, 8),
        ("satellites", 16, args.frames, 16),
        ("satellites", 24, args.frames, 24),
        ("window_frames", 50, 50, args.satellites),
        ("window_frames", 100, 100, args.satellites),
        ("window_frames", 150, 150, args.satellites),
    ]
    output: list[dict[str, Any]] = []
    for dimension, value, frames, satellites in cases:
        request, _ = generate_synthetic_input(
            label=args.label,
            frames=frames,
            satellites=satellites,
            seed=args.seed,
        )
        for _ in range(2):
            pipeline.authenticate(request)
        rows, batch = _measure_iterations(pipeline, request, args.scaling_iterations)
        latency = _percentile_summary([float(row["wall_ms"]) for row in rows])
        output.append(
            {
                "dimension": dimension,
                "value": value,
                "frames": frames,
                "satellites": satellites,
                "window_seconds": frames / 25.0,
                "iterations": args.scaling_iterations,
                "mean_ms": latency["mean"],
                "median_ms": latency["median"],
                "p95_ms": latency["p95"],
                "p99_ms": latency["p99"],
                "throughput_inferences_per_second": batch["throughput_inferences_per_second"],
                "accepted_fraction": float(np.mean([row["accepted"] for row in rows])),
                "predicted_template": str(rows[-1]["predicted_template"]),
            }
        )
    return output


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _build_metric_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    latency = summary["steady_state"]["latency_ms"]
    batch = summary["steady_state"]["batch_cpu_and_throughput"]
    memory = summary["memory"]
    rows = [
        {"category": "latency", "metric": key, "value": value, "unit": "ms"}
        for key, value in latency.items()
    ]
    rows.extend(
        [
            {
                "category": "cpu",
                "metric": "cpu_time_per_inference",
                "value": batch["cpu_time_per_inference_ms"],
                "unit": "ms",
            },
            {
                "category": "cpu",
                "metric": "one_core_equivalent",
                "value": batch["cpu_core_equivalent_pct"],
                "unit": "%",
            },
            {
                "category": "cpu",
                "metric": "whole_machine_normalized",
                "value": batch["cpu_system_normalized_pct"],
                "unit": "%",
            },
            {
                "category": "throughput",
                "metric": "sequential_inferences",
                "value": batch["throughput_inferences_per_second"],
                "unit": "inferences/s",
            },
        ]
    )
    rows.extend(
        {"category": "memory", "metric": key, "value": value, "unit": "MiB"}
        for key, value in memory.items()
        if key.endswith("_mib")
    )
    if summary.get("cold_start"):
        cold = summary["cold_start"]["summary"]
        rows.extend(
            [
                {
                    "category": "cold_start",
                    "metric": "first_inference_mean",
                    "value": cold["first_inference_ms"]["mean"],
                    "unit": "ms",
                },
                {
                    "category": "cold_start",
                    "metric": "pipeline_to_decision_mean",
                    "value": cold["pipeline_to_decision_ms"]["mean"],
                    "unit": "ms",
                },
            ]
        )
    return rows


def _paper_markdown(summary: dict[str, Any]) -> str:
    env = summary["environment"]
    inp = summary["input"]
    steady = summary["steady_state"]
    latency = steady["latency_ms"]
    cpu = steady["batch_cpu_and_throughput"]
    memory = summary["memory"]
    stage = steady["stage_latency_ms"]
    affinity = env["cpu_affinity_logical_processors"]
    if isinstance(affinity, list):
        affinity_zh = f"为避免混合架构的核心迁移，进程固定在逻辑处理器 {','.join(map(str, affinity))}"
        affinity_en = (
            f" Process affinity was restricted to logical processors {','.join(map(str, affinity))} "
            "to avoid P-core/E-core migration."
        )
    else:
        affinity_zh = ""
        affinity_en = ""
    cold_text = "not measured"
    if summary.get("cold_start"):
        cold_text = f"{summary['cold_start']['summary']['first_inference_ms']['mean']:.2f} ms"

    lines = [
        "# StarDial Runtime Results",
        "",
        "## Paper-ready table",
        "",
        "| Implementation / input | Mean (ms) | Median (ms) | P95 (ms) | P99 (ms) | CPU time (ms) | CPU (% of one core) | CPU (% of machine) | Peak RSS (MiB) | Throughput (auth/s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        (
            f"| Python/NumPy reference; {inp['frames']} frames, {inp['satellites']} satellites, "
            f"{inp['template_count']} templates | {latency['mean']:.2f} | {latency['median']:.2f} | "
            f"{latency['p95']:.2f} | {latency['p99']:.2f} | "
            f"{cpu['cpu_time_per_inference_ms']:.2f} | {cpu['cpu_core_equivalent_pct']:.2f} | "
            f"{cpu['cpu_system_normalized_pct']:.2f} | {memory['sampled_peak_rss_mib']:.2f} | "
            f"{cpu['throughput_inferences_per_second']:.2f} |"
        ),
        "",
        "## Stage breakdown",
        "",
        "| Online stage | Mean (ms) | Median (ms) | P95 (ms) | Share of internal latency (%) |",
        "|---|---:|---:|---:|---:|",
    ]
    total_mean = stage["total_internal"]["mean"]
    stage_labels = {
        "preprocessing": "C/N0 preprocessing",
        "diffraction_features": "Diffraction feature extraction",
        "geometric_inversion": "LoS geometric inversion",
        "template_matching": "Template matching (12 templates)",
        "physical_validation": "Physical consistency validation",
        "joint_decision": "Joint decision",
    }
    for key, label in stage_labels.items():
        value = stage[key]
        share = 100.0 * value["mean"] / total_mean
        lines.append(
            f"| {label} | {value['mean']:.3f} | {value['median']:.3f} | {value['p95']:.3f} | {share:.2f} |"
        )
    if summary.get("cold_start"):
        cold = summary["cold_start"]["summary"]
        lines.extend(
            [
                "",
                "## Cold-start breakdown",
                "",
                "| Cold-start metric | Mean (ms) | P95 (ms) |",
                "|---|---:|---:|",
                f"| First inference | {cold['first_inference_ms']['mean']:.2f} | {cold['first_inference_ms']['p95']:.2f} |",
                f"| Pipeline initialization + first inference | {cold['pipeline_to_decision_ms']['mean']:.2f} | {cold['pipeline_to_decision_ms']['p95']:.2f} |",
                f"| New process start to exit | {cold['process_start_to_exit_ms']['mean']:.2f} | {cold['process_start_to_exit_ms']['p95']:.2f} |",
            ]
        )
    if summary.get("scaling"):
        lines.extend(
            [
                "",
                "## Input scaling",
                "",
                "| Dimension | Frames | Satellites | Mean (ms) | P95 (ms) | Throughput (auth/s) | Accepted fraction |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in summary["scaling"]:
            lines.append(
                f"| {row['dimension']}={row['value']} | {row['frames']} | {row['satellites']} | "
                f"{row['mean_ms']:.2f} | {row['p95_ms']:.2f} | "
                f"{row['throughput_inferences_per_second']:.2f} | {row['accepted_fraction']:.2f} |"
            )
    lines.extend(
        [
            "",
            "## 中文论文表述（可直接修改后使用）",
            "",
            (
                f"我们在配备 {env['cpu_model']} 和 {env['total_memory_gib']:.1f} GiB 内存的主机上，"
                f"采用 Python {env['python_version']} / NumPy {env['numpy_version']} 对 StarDial 的精简在线认证流水线进行开销测试。"
                f"{affinity_zh}。"
                f"输入为 {inp['sample_rate_hz']:.0f} Hz、{inp['window_seconds']:.1f} s 的观测窗口，"
                f"包含 {inp['satellites']} 颗可见卫星，并与 {inp['template_count']} 个、每个 {inp['compare_points']} 点的手势模板进行匹配。"
                f"在预热 {steady['warmup_iterations']} 次后连续执行 {steady['measurement_iterations']} 次，"
                f"端到端认证时延的均值、中位数和 P95 分别为 {latency['mean']:.2f} ms、"
                f"{latency['median']:.2f} ms 和 {latency['p95']:.2f} ms；顺序吞吐量为 "
                f"{cpu['throughput_inferences_per_second']:.2f} 次认证/秒。每次认证平均消耗 "
                f"{cpu['cpu_time_per_inference_ms']:.2f} ms CPU 时间，相当于单个逻辑核的 "
                f"{cpu['cpu_core_equivalent_pct']:.2f}%（折合整机 {cpu['cpu_system_normalized_pct']:.2f}%）；"
                f"实测进程峰值驻留内存为 {memory['sampled_peak_rss_mib']:.2f} MiB。"
                f"首次推理的平均时延为 {cold_text}。"
            ),
            "",
            "## English paper wording",
            "",
            (
                f"We measured the computational overhead of the streamlined StarDial online authentication pipeline on a host equipped with "
                f"{env['cpu_model']} and {env['total_memory_gib']:.1f} GiB RAM, using Python {env['python_version']} and NumPy {env['numpy_version']}. "
                f"{affinity_en} "
                f"Each input contains a {inp['window_seconds']:.1f}-s window sampled at {inp['sample_rate_hz']:.0f} Hz from {inp['satellites']} visible satellites, "
                f"and is compared against {inp['template_count']} gesture templates resampled to {inp['compare_points']} points. "
                f"After {steady['warmup_iterations']} warm-up runs, {steady['measurement_iterations']} sequential trials achieved mean, median, and 95th-percentile "
                f"end-to-end latencies of {latency['mean']:.2f}, {latency['median']:.2f}, and {latency['p95']:.2f} ms, respectively, corresponding to "
                f"{cpu['throughput_inferences_per_second']:.2f} authentications/s. Each inference consumed {cpu['cpu_time_per_inference_ms']:.2f} ms of CPU time "
                f"({cpu['cpu_core_equivalent_pct']:.2f}% of one logical core, or {cpu['cpu_system_normalized_pct']:.2f}% of the full machine), and the sampled peak resident set size was "
                f"{memory['sampled_peak_rss_mib']:.2f} MiB."
            ),
            "",
            "## Measurement boundary and caveat",
            "",
            "Steady-state latency includes C/N0 cleaning, diffraction feature extraction, event-conditioned LoS geometric inversion, 12-template RMSE/DTW/shape scoring, physical consistency validation, and the final joint decision. RINEX parsing, ephemeris download, plots, file I/O, synthetic-input construction, and offline threshold calibration are excluded.",
            "",
            "These numbers describe the isolated Python/NumPy reference implementation in `runtime_benchmark`, not the original MATLAB research scripts. The benchmark uses deterministic synthetic observations because raw RINEX obs/nav files are not stored in this repository. Results should be rerun on the final deployment receiver before publication claims are frozen.",
            "",
        ]
    )
    return "\n".join(lines)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=200)
    parser.add_argument("--warmup", type=int, default=15)
    parser.add_argument("--scaling-iterations", type=int, default=12)
    parser.add_argument("--cold-runs", type=int, default=5)
    parser.add_argument("--frames", type=int, default=100)
    parser.add_argument("--satellites", type=int, default=16)
    parser.add_argument("--label", default="Star")
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument(
        "--cpu-affinity",
        default="",
        help="optional logical CPUs, e.g. 0-11 or 2,4,6 (recorded in results)",
    )
    parser.add_argument(
        "--max-latency-drift-percent",
        type=float,
        default=15.0,
        help="fail if last-quarter mean drifts this far from first-quarter mean",
    )
    parser.add_argument("--output", type=Path, default=SCRIPT_DIR / "results" / "latest")
    parser.add_argument("--skip-scaling", action="store_true")
    parser.add_argument("--skip-cold", action="store_true")
    parser.add_argument("--cold-worker", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    cpu_affinity = _parse_cpu_affinity(args.cpu_affinity)
    _set_cpu_affinity(cpu_affinity)
    if args.cold_worker:
        return _cold_worker(args)
    if args.iterations < 2 or args.warmup < 0 or args.scaling_iterations < 2:
        raise ValueError("iterations must be >=2, warmup >=0, and scaling iterations >=2")

    config = PipelineConfig()
    pipeline = StarDialPipeline(config)
    request, _ = generate_synthetic_input(
        label=args.label,
        frames=args.frames,
        satellites=args.satellites,
        seed=args.seed,
        config=config,
    )

    genuine = pipeline.authenticate(request)
    wrong_claim = "A" if args.label != "A" else "Star"
    impostor_request = AuthenticationInput(
        cn0_dbhz=request.cn0_dbhz,
        los_contact_xy_m=request.los_contact_xy_m,
        claimed_template=wrong_claim,
        sample_rate_hz=request.sample_rate_hz,
    )
    impostor = pipeline.authenticate(impostor_request)
    if not genuine.accepted or impostor.accepted:
        raise RuntimeError(
            "correctness smoke test failed: expected genuine accept and wrong-claim reject"
        )

    for _ in range(args.warmup):
        pipeline.authenticate(request)
    raw_rows, batch = _measure_iterations(pipeline, request, args.iterations)
    latency = _percentile_summary([float(row["wall_ms"]) for row in raw_rows])
    quarter = max(1, args.iterations // 4)
    first_quarter_mean = float(
        np.mean([float(row["wall_ms"]) for row in raw_rows[:quarter]])
    )
    last_quarter_mean = float(
        np.mean([float(row["wall_ms"]) for row in raw_rows[-quarter:]])
    )
    latency_drift_percent = 100.0 * (last_quarter_mean - first_quarter_mean) / first_quarter_mean
    stationarity = {
        "first_quarter_mean_ms": first_quarter_mean,
        "last_quarter_mean_ms": last_quarter_mean,
        "last_vs_first_drift_percent": latency_drift_percent,
        "maximum_allowed_absolute_drift_percent": args.max_latency_drift_percent,
        "passed": abs(latency_drift_percent) <= args.max_latency_drift_percent,
    }
    if not stationarity["passed"]:
        raise RuntimeError(
            f"non-stationary latency: first quarter {first_quarter_mean:.2f} ms, "
            f"last quarter {last_quarter_mean:.2f} ms ({latency_drift_percent:+.1f}%)"
        )
    stage_names = [
        key.removeprefix("stage_").removesuffix("_ms")
        for key in raw_rows[0]
        if key.startswith("stage_")
    ]
    stage_latency = {
        stage: _percentile_summary(
            [float(row[f"stage_{stage}_ms"]) for row in raw_rows]
        )
        for stage in stage_names
    }
    memory = _measure_memory(pipeline, request)
    cold = None if args.skip_cold else _measure_cold_start(args)
    scaling = [] if args.skip_scaling else _measure_scaling(pipeline, args)

    summary: dict[str, Any] = {
        "scope": {
            "implementation": "isolated Python/NumPy StarDial reference pipeline",
            "included": [
                "C/N0 preprocessing",
                "diffraction feature extraction",
                "event-conditioned LoS geometric inversion",
                "RMSE + DTW + structured template scoring",
                "physical consistency validation",
                "joint authentication decision",
            ],
            "excluded": [
                "RINEX parsing",
                "ephemeris acquisition",
                "synthetic input generation",
                "plotting and file I/O",
                "offline threshold calibration",
            ],
        },
        "environment": _environment_info(cpu_affinity),
        "input": {
            "label": args.label,
            "frames": args.frames,
            "satellites": args.satellites,
            "sample_rate_hz": request.sample_rate_hz,
            "window_seconds": args.frames / request.sample_rate_hz,
            "template_count": 12,
            "compare_points": config.compare_points,
            "seed": args.seed,
            "source": "deterministic synthetic C/N0 coupled to time-indexed LoS geometry",
        },
        "configuration": asdict(config),
        "correctness_smoke_test": {
            "genuine": result_to_dict(genuine),
            "wrong_claim": result_to_dict(impostor),
        },
        "steady_state": {
            "warmup_iterations": args.warmup,
            "measurement_iterations": args.iterations,
            "latency_ms": latency,
            "stage_latency_ms": stage_latency,
            "stationarity_check": stationarity,
            "batch_cpu_and_throughput": batch,
            "accepted_fraction": float(np.mean([row["accepted"] for row in raw_rows])),
        },
        "memory": memory,
        "cold_start": cold,
        "scaling": scaling,
    }

    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "raw_iterations.csv", raw_rows)
    _write_csv(output_dir / "summary_metrics.csv", _build_metric_rows(summary))
    _write_csv(output_dir / "scaling.csv", scaling)
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    (output_dir / "PAPER_RUNTIME_TEXT.md").write_text(
        _paper_markdown(summary), encoding="utf-8"
    )

    print(json.dumps({"output": str(output_dir), "latency_ms": latency, "batch": batch}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
