# StarDial Runtime Results

## Paper-ready table

| Implementation / input | Mean (ms) | Median (ms) | P95 (ms) | P99 (ms) | CPU time (ms) | CPU (% of one core) | CPU (% of machine) | Peak RSS (MiB) | Throughput (auth/s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Python/NumPy reference; 100 frames, 16 satellites, 12 templates | 100.26 | 100.31 | 104.39 | 106.75 | 97.03 | 96.78 | 4.84 | 109.52 | 9.97 |

## Stage breakdown

| Online stage | Mean (ms) | Median (ms) | P95 (ms) | Share of internal latency (%) |
|---|---:|---:|---:|---:|
| C/N0 preprocessing | 0.330 | 0.306 | 0.440 | 0.33 |
| Diffraction feature extraction | 0.158 | 0.143 | 0.242 | 0.16 |
| LoS geometric inversion | 79.125 | 79.203 | 82.507 | 78.96 |
| Template matching (12 templates) | 20.232 | 20.107 | 21.578 | 20.19 |
| Physical consistency validation | 0.362 | 0.343 | 0.470 | 0.36 |
| Joint decision | 0.004 | 0.003 | 0.005 | 0.00 |

## Cold-start breakdown

| Cold-start metric | Mean (ms) | P95 (ms) |
|---|---:|---:|
| First inference | 103.41 | 104.36 |
| Pipeline initialization + first inference | 116.66 | 117.46 |
| New process start to exit | 262.85 | 267.09 |

## Input scaling

| Dimension | Frames | Satellites | Mean (ms) | P95 (ms) | Throughput (auth/s) | Accepted fraction |
|---|---:|---:|---:|---:|---:|---:|
| satellites=8 | 100 | 8 | 63.45 | 64.35 | 15.76 | 1.00 |
| satellites=16 | 100 | 16 | 100.03 | 101.67 | 10.00 | 1.00 |
| satellites=24 | 100 | 24 | 136.25 | 138.52 | 7.34 | 1.00 |
| window_frames=50 | 50 | 16 | 59.76 | 60.79 | 16.73 | 1.00 |
| window_frames=100 | 100 | 16 | 99.32 | 101.75 | 10.07 | 1.00 |
| window_frames=150 | 150 | 16 | 139.19 | 141.69 | 7.18 | 1.00 |

## 中文论文表述（可直接修改后使用）

我们在配备 Intel(R) Core(TM) i5-14500 和 31.7 GiB 内存的主机上，采用 Python 3.12.10 / NumPy 2.5.1 对 StarDial 的精简在线认证流水线进行开销测试。为避免混合架构的核心迁移，进程固定在逻辑处理器 0,1,2,3,4,5,6,7,8,9,10,11。输入为 25 Hz、4.0 s 的观测窗口，包含 16 颗可见卫星，并与 12 个、每个 160 点的手势模板进行匹配。在预热 15 次后连续执行 200 次，端到端认证时延的均值、中位数和 P95 分别为 100.26 ms、100.31 ms 和 104.39 ms；顺序吞吐量为 9.97 次认证/秒。每次认证平均消耗 97.03 ms CPU 时间，相当于单个逻辑核的 96.78%（折合整机 4.84%）；实测进程峰值驻留内存为 109.52 MiB。首次推理的平均时延为 103.41 ms。

## English paper wording

We measured the computational overhead of the streamlined StarDial online authentication pipeline on a host equipped with Intel(R) Core(TM) i5-14500 and 31.7 GiB RAM, using Python 3.12.10 and NumPy 2.5.1.  Process affinity was restricted to logical processors 0,1,2,3,4,5,6,7,8,9,10,11 to avoid P-core/E-core migration. Each input contains a 4.0-s window sampled at 25 Hz from 16 visible satellites, and is compared against 12 gesture templates resampled to 160 points. After 15 warm-up runs, 200 sequential trials achieved mean, median, and 95th-percentile end-to-end latencies of 100.26, 100.31, and 104.39 ms, respectively, corresponding to 9.97 authentications/s. Each inference consumed 97.03 ms of CPU time (96.78% of one logical core, or 4.84% of the full machine), and the sampled peak resident set size was 109.52 MiB.

## Measurement boundary and caveat

Steady-state latency includes C/N0 cleaning, diffraction feature extraction, event-conditioned LoS geometric inversion, 12-template RMSE/DTW/shape scoring, physical consistency validation, and the final joint decision. RINEX parsing, ephemeris download, plots, file I/O, synthetic-input construction, and offline threshold calibration are excluded.

These numbers describe the isolated Python/NumPy reference implementation in `runtime_benchmark`, not the original MATLAB research scripts. The benchmark uses deterministic synthetic observations because raw RINEX obs/nav files are not stored in this repository. Results should be rerun on the final deployment receiver before publication claims are frozen.
