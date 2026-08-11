# StarDial online authentication runtime benchmark

This directory contains a small, dependency-light reference implementation of
the online authentication path described in Section III of the StarDial paper.
It is isolated from the repository's MATLAB plotting, gallery generation, and
offline evaluation logic so those activities cannot contaminate runtime data.

## What is measured

One steady-state inference contains:

1. baseline-based C/N0 cleaning;
2. amplitude/oscillation/continuity diffraction evidence extraction;
3. event-conditioned LoS geometric inversion with exclusion evidence;
4. comparison against the repository's 12 gesture templates using 160-point
   RMSE, exact DTW, and the structured feature penalty;
5. event, exclusion, waveform, and temporal physical-consistency validation;
6. the paper's joint gesture-and-physical acceptance decision.

RINEX parsing, ephemeris acquisition, synthetic input creation, plots, file
output, and offline threshold selection are outside the timed inference call.
Template and search-grid construction are also excluded from steady state, but
are included in the separately reported cold pipeline-to-decision result.

## Run

From the repository root:

```powershell
python runtime_benchmark\benchmark_runtime.py
```

The default experiment uses a deterministic 4 s / 25 Hz window, 16 satellites,
12 templates, 15 warm-up runs, 200 measured runs, five clean-process cold runs,
and two scaling sweeps. It requires only Python and NumPy; Windows memory data
are read directly from the operating system.

On a hybrid P-core/E-core CPU, pin the process to the performance-core logical
processors to avoid scheduler migration changing the result halfway through a
run. For the i5-14500 test host used here, logical processors 0-11 are the six
P-cores (including their SMT siblings):

```powershell
python runtime_benchmark\benchmark_runtime.py --cpu-affinity 0-11
```

The selected logical processors are stored in `summary.json`. The benchmark
also rejects a run automatically if the last-quarter mean latency differs from
the first-quarter mean by more than 15%.

For a quick functional check:

```powershell
python -m unittest runtime_benchmark\test_pipeline.py
python runtime_benchmark\benchmark_runtime.py --iterations 10 --warmup 2 --skip-cold --skip-scaling
```

## Outputs

The default output directory is `runtime_benchmark/results/latest`:

- `raw_iterations.csv`: every latency, CPU, decision, and stage sample;
- `summary_metrics.csv`: compact headline metrics;
- `scaling.csv`: satellite-count and observation-window scaling;
- `summary.json`: complete configuration, environment, correctness, and result
  record;
- `PAPER_RUNTIME_TEXT.md`: paper-ready table plus Chinese and English wording.

CPU time is measured using the process CPU clock. “CPU (% of one core)” is
`process CPU time / wall time`; 100% means one logical CPU was fully occupied.
“CPU (% of machine)” divides that value by the number of logical CPUs and is
therefore comparable to whole-machine utilization shown by common system
monitors. Memory and latency are measured in separate passes.

## Scientific scope

The repository does not contain raw RINEX obs/nav recordings, and MATLAB or
Octave is not installed on this machine. The benchmark therefore uses
fixed-seed synthetic C/N0 observations whose diffraction strength is coupled to
time-indexed LoS contact geometry. The ground-truth trajectory is used only to
generate the observations; it is never supplied to the authenticator.

Consequently, reported numbers are the cost of this Python/NumPy reference
pipeline, not the cost of the original MATLAB research scripts or of RINEX
parsing. Rerun the same benchmark with real parsed observation arrays, and on
the final receiver hardware, before freezing a publication claim.
