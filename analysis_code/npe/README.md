# Cs137 NPE Analysis Tool

This directory contains the Cs137 liquid-scintillator NPE comparison workflow.

## Files

```text
plot_npe_subtracted.C   ROOT macro that builds source/BG overlays and subtraction plots
run_npe_analysis.sh     Command-line wrapper around the ROOT macro
run_npe_gui.sh          Recommended GUI launcher
run_npe_gui.py          Tkinter desktop GUI
run_npe_web.py          Browser GUI fallback, no external Python packages
setup_env.sh            Optional environment check/venv helper
```

## Quick Start

From the repository root:

```bash
analysis_code/npe/run_npe_gui.sh
```

If Tkinter is available, a desktop GUI opens. If not, the launcher starts a local browser GUI and prints a URL like:

```text
http://127.0.0.1:8765
```

## Command-Line Usage

```bash
analysis_code/npe/run_npe_analysis.sh \
  -d /path/to/data \
  -s Cs137_nocollimator.root \
  -b background_1hr_prod.root \
  -O /path/to/results \
  -o nocollimator
```

Useful options:

```text
-x MIN -X MAX   Manual NPE x-axis range
-n BINS         Histogram bin count
-g GAIN         PMT gain, default 1.0e7
-q QUANTILE     X-axis quantile, 1.0 means full range
```

Example zoomed plot:

```bash
analysis_code/npe/run_npe_analysis.sh \
  -d /path/to/data \
  -s Cs137_nocollimator.root \
  -b background_1hr_prod.root \
  -O /path/to/results \
  -o nocollimator_zoom \
  -x 0 -X 200000
```

## Main Outputs

```text
<prefix>_overlay_total_rate_log.png   Live-time and bin-width normalized source/BG overlay
<prefix>_overlay_total_log.png        Raw-count source/BG overlay
<prefix>_subtracted_total.png         Source - scaled background
<prefix>_overlay_channels_log.png     Channel-by-channel source/BG overlay
<prefix>_subtracted_channels.png      Channel-by-channel subtraction
```

Use `*_overlay_total_rate_log.png` first when comparing runs with different trigger rates or live times.
