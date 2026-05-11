# SPE Calibration Macros

ROOT macros for single-photoelectron calibration studies.

```text
calculate_charge.C
view_waveforms.C
fitting/
  fit_all_channels.C
  fit_all_spe.C
  fit_channels_dark_led.C
  fit_led_only.C
  fit_spe.C
```

Raw SPE ROOT files and generated fit plots should stay outside `analysis_code/`, for example under `SPE/` locally or an external data directory.
