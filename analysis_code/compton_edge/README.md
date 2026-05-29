# Compton Edge Analysis Tool

This directory contains the Cs137 source/background Compton-edge workflow built from liquid-scintillator NPE spectra.

## Files

```text
plot_compton_edge.C   ROOT macro that builds source/BG overlays and subtraction plots
run_compton_edge_analysis.sh     Command-line wrapper around the ROOT macro
run_compton_edge_fit.sh          Compton-edge erfc fit wrapper for generated histograms
run_compton_edge_scan.sh         Fit range stability scan wrapper
run_compton_edge_model_compare.sh
                                  Compares erfc_linear and erfc_gaussian fit results
run_compton_edge_gui.sh          Recommended GUI launcher
run_compton_edge_gui.py          Tkinter desktop GUI
run_compton_edge_web.py          Browser GUI fallback, no external Python packages
check_env.sh            Checks Python, Tkinter, ROOT, and ROOT macro loading
setup_env.sh            Optional environment check/venv helper
```

## Environment

Recommended runtime:

```text
ROOT    6.36 stable series for reproducible analysis
Python  3.10 or newer
Bash    4 or newer
Tkinter optional; without it, run_compton_edge_gui.sh uses the browser GUI fallback
```

Check the local environment before analyzing data:

```bash
analysis_code/compton_edge/check_env.sh
```

If Tkinter is missing, install the OS package such as `python3-tk` on Ubuntu/Debian. ROOT must be available as the `root` command in `PATH`.

## Quick Start

From the repository root:

```bash
analysis_code/compton_edge/run_compton_edge_gui.sh
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

`Source Label` controls the source name shown in plot legends. If left empty in command-line use, the label is derived from the source ROOT file name.

## Command-Line Usage

```bash
analysis_code/compton_edge/run_compton_edge_analysis.sh \
  -d /path/to/data \
  -s source.root \
  -b background_1hr_prod.root \
  -O /path/to/results \
  -o source_run
```

Useful options:

```text
-x MIN -X MAX   Manual NPE x-axis range
-n BINS         Histogram bin count
-g GAIN         PMT gain, default 1.0e7
-q QUANTILE     X-axis quantile, 1.0 means full range
-L LABEL        Source name shown in plot legends
-C LIST         Channels to analyze, e.g. 0,1
-T LIST         Event NPE thresholds, e.g. 0:5,1:5
```

## Compton Edge Fit

After creating `<prefix>_histograms.root`, fit the total background-subtracted histogram with a linear background plus smeared edge:

```text
f(x) = offset + slope*x + 0.5*amplitude*erfc((x - edge)/(sqrt(2)*sigma))
```

Command-line example:

```bash
analysis_code/compton_edge/run_compton_edge_fit.sh \
  -i /path/to/results/source_zoom_histograms.root \
  -x 3000 -X 6000 \
  -o /path/to/results/source_zoom \
  -M erfc_linear
```

Available models:

```text
erfc_linear     linear background + smeared falling edge
erfc_gaussian   erfc_linear plus one Gaussian component for local bump/peak structure
```

The fit plot includes a lower pull panel showing `(data - fit) / error` across the selected fit range. This makes it easier to see where the model misses the spectrum even when the overlay looks acceptable.

For `erfc_gaussian`, the fit plot also draws the `linear + erfc` component, the Gaussian contribution, and the Gaussian mean marker separately. Use these component curves to check whether the Gaussian term is modeling a physically intended peak or merely absorbing continuum-shape residuals.

The fit writes:

```text
<prefix>_compton_edge_fit.png    Histogram with fit curve and edge marker
<prefix>_compton_edge_fit.pdf    Vector version for publication workflows
<prefix>_compton_edge_fit.root   Fit input histogram, TF1, and canvas
<prefix>_compton_edge_fit.txt    Fit parameters and errors
```

To estimate the fit-range dependence, run the range scan. It repeats the same fit after shifting the lower and upper fit bounds by a fraction of the selected fit width:

```bash
analysis_code/compton_edge/run_compton_edge_scan.sh \
  -i /path/to/results/source_zoom_histograms.root \
  -x 3000 -X 6000 \
  -o /path/to/results/source_zoom \
  -M erfc_linear \
  -F 0.10
```

The scan writes:

```text
<prefix>_compton_edge_scan.png   Edge value for each shifted fit range
<prefix>_compton_edge_scan.pdf   Vector version for publication workflows
<prefix>_compton_edge_scan.root  Scan input histogram, graph, and canvas
<prefix>_compton_edge_scan.txt   Range table and edge mean/RMS summary
```

Use `edge_rms` from the scan summary as a first estimate of the systematic uncertainty from fit-range choice. It is not a replacement for model comparison, but it quickly shows whether the selected range is stable enough to trust.

To check model dependence, compare the two available fit models on the same histogram and range:

```bash
analysis_code/compton_edge/run_compton_edge_model_compare.sh \
  -i /path/to/results/source_zoom_histograms.root \
  -x 3000 -X 6000 \
  -o /path/to/results/source_zoom
```

This runs both `erfc_linear` and `erfc_gaussian`, saves each fit plot with a model-specific prefix, and writes:

```text
<prefix>_compton_edge_model_compare.txt
```

Use `edge_diff_abs` as a first estimate of model-choice sensitivity. If it is comparable to or larger than the statistical edge error, inspect both fit plots and the pull panels before quoting the edge.

The fit plot also marks the maximum falling-slope position. This derivative-based edge estimate is used as the initial edge seed and is written as `slope_edge` in the text summary, so it can be compared against the final fitted `edge`.

The fit summary includes warning tags when the result needs inspection:

```text
fit_status_N                         ROOT fit did not return status 0
edge_near_fit_boundary               fitted edge is too close to the selected fit range boundary
edge_error_large_or_invalid          edge uncertainty is zero, negative, or very large
sigma_at_limit_or_too_large          edge-smearing width is at a fit limit or too broad
large_chi2_ndf                       chi2/ndf is greater than 5
large_pull_rms                       pull RMS is greater than 3
gaussian_sigma_at_limit_or_too_large Gaussian component width is at a limit or too broad
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

If the digitizer range, sampling period, impedance, or ADC bit depth changes, update the optional constants in `plot_compton_edge.C`.

Example zoomed plot:

```bash
analysis_code/compton_edge/run_compton_edge_analysis.sh \
  -d /path/to/data \
  -s source.root \
  -b background_1hr_prod.root \
  -O /path/to/results \
  -o source_zoom \
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
