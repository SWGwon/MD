# NPE Analysis Tool

This directory contains the liquid-scintillator NPE comparison workflow.

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

The GUI always runs and previews two total-subtraction plots:

```text
<prefix>_full_subtracted_total.png        Full-range result
<prefix>_range_<xmin>_<xmax>_subtracted_total.png
                                           User-selected NPE range
```

Use the full-range preview first, then adjust `Range X Min` and `Range X Max` to inspect the region of interest for the source being analyzed.

By default the selected range is also full range (`X Max = -1`). Set a finite `X Max` after inspecting the full-range preview.

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

## Charge and NPE Conversion

For `phys_tree` files produced by `/home/sgwon/CPNR_dt5730s/src/production_dt5730.cpp`, `Charge_CH*` is not pC. It is a baseline-subtracted ADC integral:

```text
Charge_CH = sum(baseline - ADC) over the pulse window
unit      = ADC count * sample
```

The macro converts this to charge and NPE using the DT5730 assumptions:

```text
Q[pC] = Charge_CH * (2.0 V / (2^14 - 1)) * (2 ns) / (50 ohm) * 1e12
NPE   = Q[pC] / (gain * 1.60217663e-7 pC)
```

With the default `gain = 1e7`, this is:

```text
Q[pC] = Charge_CH * 0.00488311
NPE   = Charge_CH * 0.00304780
```

If the digitizer range, sampling period, impedance, or ADC bit depth changes, update the optional constants in `plot_npe_subtracted.C`.

Example zoomed plot:

```bash
analysis_code/npe/run_npe_analysis.sh \
  -d /path/to/data \
  -s Cs137_nocollimator.root \
  -b background_1hr_prod.root \
  -O /path/to/results \
  -o nocollimator_zoom \
  -x 0 -X 6000
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
