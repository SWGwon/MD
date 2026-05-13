# MD Analysis

ROOT-based analysis tools and experiment-specific macros for PMT/liquid-scintillator measurements.

## Repository Layout

```text
analysis_code/
  view_waveforms.C      # General waveform viewer
  npe/                  # Cs137 NPE/background-subtraction analysis tool
  reflection_legacy/    # Older reflection-test helper macros kept for reference
  spe/                  # SPE calibration macros

SPE/
  ...                   # Local SPE data/plots only. Ignored by git except already tracked history.

docs/
  ...                   # Notes, run logs, and analysis documentation

data/                   # Local data only. Ignored by git.
results/                # Local plots/outputs only. Ignored by git.
```

## Data Policy

Large or generated files are intentionally not tracked:

- ROOT files (`*.root`)
- plots (`*.png`, `*.pdf`)
- ROOT ACLiC/build products (`*_C.so`, `*_C.d`, `*.pcm`)
- local outputs under `data/` and `results/`

Keep reproducible analysis code in `analysis_code/`, and keep raw/processed run files under `data/` or another external data directory.

## Cs137 NPE Analysis

The current user-facing analysis entry point is:

```bash
analysis_code/npe/run_npe_gui.sh
```

It automatically uses the Tkinter desktop GUI when available, otherwise it falls back to a browser-based GUI that only needs standard Python.

Command-line usage is also available:

```bash
analysis_code/npe/run_npe_analysis.sh \
  -d /path/to/data \
  -s source.root \
  -b background.root \
  -O /path/to/results \
  -o run_prefix
```

See [analysis_code/npe/README.md](analysis_code/npe/README.md) for details.
