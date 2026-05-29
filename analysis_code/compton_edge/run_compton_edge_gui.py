#!/usr/bin/env python3
import subprocess
import threading
import math
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


APP_DIR = Path(__file__).resolve().parent
REPO_DIR = APP_DIR.parents[1]
DEFAULT_DATA_DIR = REPO_DIR / "data"
DEFAULT_OUT_DIR = REPO_DIR / "results"
WRAPPER = APP_DIR / "run_compton_edge_analysis.sh"
FIT_WRAPPER = APP_DIR / "run_compton_edge_fit.sh"
SCAN_WRAPPER = APP_DIR / "run_compton_edge_scan.sh"
COMPARE_WRAPPER = APP_DIR / "run_compton_edge_model_compare.sh"


class ComptonEdgeAnalysisGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Compton Edge Analysis")
        self.geometry("1180x820")
        self.minsize(1040, 740)
        self.process = None
        self.preview_images = []
        self.result_paths = []
        self.fit_result_path = None
        self.scan_result_path = None
        self.model_compare_paths = []

        self.source_var = tk.StringVar()
        self.bg_var = tk.StringVar()
        self.out_dir_var = tk.StringVar(value=str(DEFAULT_OUT_DIR))
        self.prefix_var = tk.StringVar(value="compton_edge")
        self.source_label_var = tk.StringVar()
        self.gain_var = tk.StringVar(value="1.0e7")
        self.quantile_var = tk.StringVar(value="1.0")
        self.bins_var = tk.StringVar(value="400")
        self.xmin_var = tk.StringVar(value="0")
        self.xmax_var = tk.StringVar(value="-1")
        self.clock_var = tk.StringVar(value="125.0e6")
        self.fit_xmin_var = tk.StringVar()
        self.fit_xmax_var = tk.StringVar()
        self.fit_model_var = tk.StringVar(value="erfc_linear")
        self.fit_result_var = tk.StringVar(value="No fit has been run.")
        self.channel_vars = [tk.BooleanVar(value=True) for _ in range(8)]
        self.threshold_vars = [tk.StringVar() for _ in range(8)]

        self._build()

    def _build(self):
        root = ttk.Frame(self, padding=12)
        root.pack(fill=tk.BOTH, expand=True)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(9, weight=1)
        root.rowconfigure(11, weight=1)

        self._file_row(root, 0, "Source ROOT", self.source_var, self.pick_source)
        self._file_row(root, 1, "Background ROOT (optional)", self.bg_var, self.pick_background)
        self._file_row(root, 2, "Output Dir", self.out_dir_var, self.pick_output_dir)

        ttk.Label(root, text="Output Prefix").grid(row=3, column=0, sticky="w", pady=4)
        ttk.Entry(root, textvariable=self.prefix_var).grid(row=3, column=1, columnspan=2, sticky="ew", pady=4)

        ttk.Label(root, text="Source Label").grid(row=4, column=0, sticky="w", pady=4)
        ttk.Entry(root, textvariable=self.source_label_var).grid(row=4, column=1, columnspan=2, sticky="ew", pady=4)

        options = ttk.LabelFrame(root, text="Options", padding=8)
        options.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(8, 4))
        for col in range(6):
            options.columnconfigure(col, weight=1)

        self._option(options, 0, 0, "Gain", self.gain_var)
        self._option(options, 0, 2, "X Quantile", self.quantile_var)
        self._option(options, 0, 4, "Bins", self.bins_var)
        self._option(options, 1, 0, "Range X Min (NPE)", self.xmin_var)
        self._option(options, 1, 2, "Range X Max (NPE)", self.xmax_var)
        self._option(options, 1, 4, "TTT Clock Hz", self.clock_var)

        channel_frame = ttk.Frame(options)
        channel_frame.grid(row=2, column=0, columnspan=6, sticky="w", pady=(6, 0))
        ttk.Label(channel_frame, text="Channels").grid(row=0, column=0, sticky="w", padx=(0, 8))
        for ch, var in enumerate(self.channel_vars):
            ttk.Checkbutton(channel_frame, text=f"CH{ch}", variable=var).grid(row=0, column=ch + 1, sticky="w", padx=(0, 8))

        threshold_frame = ttk.Frame(options)
        threshold_frame.grid(row=3, column=0, columnspan=6, sticky="w", pady=(6, 0))
        ttk.Label(threshold_frame, text="Min NPE").grid(row=0, column=0, sticky="w", padx=(0, 8))
        for ch, var in enumerate(self.threshold_vars):
            cell = ttk.Frame(threshold_frame)
            cell.grid(row=0, column=ch + 1, sticky="w", padx=(0, 8))
            ttk.Label(cell, text=f"CH{ch}").grid(row=0, column=0, sticky="w")
            ttk.Entry(cell, textvariable=var, width=6).grid(row=1, column=0, sticky="w")

        ttk.Label(options, text="Draw Histograms creates full-range and selected-range plots. Inspect full range first, then adjust the selected range.").grid(row=4, column=0, columnspan=6, sticky="w", pady=(6, 0))

        buttons = ttk.Frame(root)
        buttons.grid(row=6, column=0, columnspan=3, sticky="ew", pady=(8, 8))
        buttons.columnconfigure(5, weight=1)

        self.run_button = ttk.Button(buttons, text="Draw Histograms", command=self.run_analysis)
        self.run_button.grid(row=0, column=0, padx=(0, 8))
        self.stop_button = ttk.Button(buttons, text="Stop", command=self.stop_analysis, state=tk.DISABLED)
        self.stop_button.grid(row=0, column=1)
        self.fit_button = ttk.Button(buttons, text="Fit Edge", command=self.run_fit, state=tk.DISABLED)
        self.fit_button.grid(row=0, column=2, padx=(8, 0))
        self.scan_button = ttk.Button(buttons, text="Scan Range", command=self.run_scan, state=tk.DISABLED)
        self.scan_button.grid(row=0, column=3, padx=(8, 0))
        self.compare_button = ttk.Button(buttons, text="Compare Models", command=self.run_model_compare, state=tk.DISABLED)
        self.compare_button.grid(row=0, column=4, padx=(8, 0))
        ttk.Button(buttons, text="Clear Log", command=self.clear_log).grid(row=0, column=6)

        fit_options = ttk.LabelFrame(root, text="Compton Edge Fit", padding=8)
        fit_options.grid(row=7, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        fit_options.columnconfigure(7, weight=1)
        self._option(fit_options, 0, 0, "Fit X Min (NPE)", self.fit_xmin_var)
        self._option(fit_options, 0, 2, "Fit X Max (NPE)", self.fit_xmax_var)
        ttk.Label(fit_options, text="Model").grid(row=0, column=4, sticky="w", padx=(0, 4), pady=4)
        ttk.Combobox(
            fit_options,
            textvariable=self.fit_model_var,
            values=["erfc_linear", "erfc_gaussian"],
            width=14,
            state="readonly",
        ).grid(row=0, column=5, sticky="w", padx=(0, 12), pady=4)
        ttk.Label(fit_options, textvariable=self.fit_result_var).grid(row=0, column=6, columnspan=2, sticky="w", padx=(8, 0))

        ttk.Label(root, text="Result Preview").grid(row=8, column=0, sticky="w")
        self.preview_frame = ttk.Frame(root)
        self.preview_frame.grid(row=9, column=0, columnspan=3, sticky="nsew")
        self.preview_frame.columnconfigure(0, weight=1)
        self.preview_frame.columnconfigure(1, weight=1)
        self.preview_frame.rowconfigure(1, weight=1)

        ttk.Label(root, text="Log").grid(row=10, column=0, sticky="w", pady=(8, 0))
        log_frame = ttk.Frame(root)
        log_frame.grid(row=11, column=0, columnspan=3, sticky="nsew")
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

        self.log = tk.Text(log_frame, wrap="word", height=16)
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log.configure(yscrollcommand=scroll.set)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(root, textvariable=self.status_var).grid(row=12, column=0, columnspan=3, sticky="ew", pady=(8, 0))

    def _file_row(self, parent, row, label, variable, command):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4, padx=(8, 8))
        ttk.Button(parent, text="Browse", command=command).grid(row=row, column=2, sticky="ew", pady=4)

    def _option(self, parent, row, col, label, variable):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=(0, 4), pady=4)
        ttk.Entry(parent, textvariable=variable, width=14).grid(row=row, column=col + 1, sticky="ew", padx=(0, 12), pady=4)

    def pick_source(self):
        path = filedialog.askopenfilename(
            title="Select source ROOT file",
            initialdir=self.default_file_dialog_dir(),
            filetypes=[("ROOT files", "*.root"), ("All files", "*")],
        )
        if path:
            self.source_var.set(path)
            if self.prefix_var.get() in {"", "npe", "compton_edge"}:
                self.prefix_var.set(Path(path).stem)
            if not self.source_label_var.get().strip():
                self.source_label_var.set(Path(path).stem)

    def pick_background(self):
        path = filedialog.askopenfilename(
            title="Select background ROOT file",
            initialdir=self.default_file_dialog_dir(),
            filetypes=[("ROOT files", "*.root"), ("All files", "*")],
        )
        if path:
            self.bg_var.set(path)

    def pick_output_dir(self):
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            self.out_dir_var.set(path)

    def default_file_dialog_dir(self):
        if DEFAULT_DATA_DIR.is_dir():
            return str(DEFAULT_DATA_DIR)
        return str(REPO_DIR)

    def validate_inputs(self):
        if not WRAPPER.exists():
            messagebox.showerror("Missing wrapper", f"Cannot find:\n{WRAPPER}")
            return False
        if not self.source_var.get() or not Path(self.source_var.get()).is_file():
            messagebox.showerror("Missing source", "Select a valid source ROOT file.")
            return False
        if self.bg_var.get().strip() and not Path(self.bg_var.get()).is_file():
            messagebox.showerror("Missing background", "Select a valid background ROOT file.")
            return False
        if not self.prefix_var.get().strip():
            messagebox.showerror("Missing prefix", "Enter an output prefix.")
            return False
        try:
            self.validated_options()
        except ValueError as exc:
            messagebox.showerror("Invalid option", str(exc))
            return False
        Path(self.out_dir_var.get()).mkdir(parents=True, exist_ok=True)
        return True

    def validated_options(self):
        gain = parse_float(self.gain_var.get(), "Gain")
        quantile = parse_float(self.quantile_var.get(), "X Quantile")
        bins = parse_int(self.bins_var.get(), "Bins")
        xmin = parse_float(self.xmin_var.get(), "Range X Min")
        xmax = parse_float(self.xmax_var.get(), "Range X Max")
        clock = parse_float(self.clock_var.get(), "TTT Clock Hz")
        channels = [str(ch) for ch, var in enumerate(self.channel_vars) if var.get()]
        thresholds = []
        for ch, var in enumerate(self.threshold_vars):
            value = var.get().strip()
            if not value:
                continue
            threshold = parse_float(value, f"CH{ch} Min NPE")
            thresholds.append(f"{ch}:{format_number(threshold)}")

        if gain <= 0:
            raise ValueError("Gain must be greater than 0.")
        if not (0 < quantile <= 1):
            raise ValueError("X Quantile must be greater than 0 and less than or equal to 1.")
        if bins <= 0:
            raise ValueError("Bins must be a positive integer.")
        if xmax != -1 and xmax <= xmin:
            raise ValueError("Range X Max must be -1 for full range or greater than Range X Min.")
        if clock <= 0:
            raise ValueError("TTT Clock Hz must be greater than 0.")
        if not channels:
            raise ValueError("Select at least one channel.")

        return {
            "gain": format_number(gain),
            "quantile": format_number(quantile),
            "bins": str(bins),
            "xmin": format_number(xmin),
            "xmax": format_number(xmax),
            "clock": format_number(clock),
            "channels": ",".join(channels),
            "thresholds": ",".join(thresholds),
        }

    def run_analysis(self):
        if self.process is not None:
            return
        if not self.validate_inputs():
            return

        base_prefix = self.prefix_var.get().strip()
        options = self.validated_options()
        range_prefix = self.range_prefix(base_prefix)
        full_prefix = f"{base_prefix}_full"

        full_cmd = [
            str(WRAPPER),
            "-s", self.source_var.get(),
            "-O", self.out_dir_var.get(),
            "-o", full_prefix,
            "-L", self.source_label_var.get().strip(),
            "-g", options["gain"],
            "-q", options["quantile"],
            "-n", "400",
            "-x", "0",
            "-X", "-1",
            "-C", options["channels"],
            "-T", options["thresholds"],
            "-c", options["clock"],
        ]

        range_cmd = [
            str(WRAPPER),
            "-s", self.source_var.get(),
            "-O", self.out_dir_var.get(),
            "-o", range_prefix,
            "-L", self.source_label_var.get().strip(),
            "-g", options["gain"],
            "-q", "1.0",
            "-n", options["bins"],
            "-x", options["xmin"],
            "-X", options["xmax"],
            "-C", options["channels"],
            "-T", options["thresholds"],
            "-c", options["clock"],
        ]
        bg_path = self.bg_var.get().strip()
        if bg_path:
            full_cmd[3:3] = ["-b", bg_path]
            range_cmd[3:3] = ["-b", bg_path]

        commands = [("full range", full_cmd), ("selected range", range_cmd)]
        self.result_paths = [
            Path(self.out_dir_var.get()) / f"{full_prefix}_subtracted_total.png",
            Path(self.out_dir_var.get()) / f"{range_prefix}_subtracted_total.png",
        ]
        self.full_hist_path = Path(self.out_dir_var.get()) / f"{full_prefix}_histograms.root"
        self.range_hist_path = Path(self.out_dir_var.get()) / f"{range_prefix}_histograms.root"
        self.fit_output_prefix = Path(self.out_dir_var.get()) / f"{range_prefix}"
        self.fit_result_path = None
        self.scan_result_path = None
        self.model_compare_paths = []
        self.fit_button.configure(state=tk.DISABLED)
        self.scan_button.configure(state=tk.DISABLED)
        self.compare_button.configure(state=tk.DISABLED)
        self.fit_result_var.set("Draw histograms before fitting.")

        self.append_log("$ " + " ".join(full_cmd) + "\n")
        self.append_log("$ " + " ".join(range_cmd) + "\n\n")
        self.status_var.set("Drawing histograms...")
        self.run_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)

        thread = threading.Thread(target=self._run_subprocess, args=(commands,), daemon=True)
        thread.start()

    def validate_fit_inputs(self):
        if not FIT_WRAPPER.exists():
            messagebox.showerror("Missing fit wrapper", f"Cannot find:\n{FIT_WRAPPER}")
            return False
        return self.validate_common_fit_inputs()

    def validate_scan_inputs(self):
        if not SCAN_WRAPPER.exists():
            messagebox.showerror("Missing scan wrapper", f"Cannot find:\n{SCAN_WRAPPER}")
            return False
        return self.validate_common_fit_inputs()

    def validate_compare_inputs(self):
        if not COMPARE_WRAPPER.exists():
            messagebox.showerror("Missing model-compare wrapper", f"Cannot find:\n{COMPARE_WRAPPER}")
            return False
        return self.validate_common_fit_inputs()

    def validate_common_fit_inputs(self):
        hist_path = self.current_fit_hist_path()
        if hist_path is None or not hist_path.is_file():
            messagebox.showerror("Missing histogram", "Draw histograms first so the histogram ROOT file is available.")
            return False
        try:
            xmin = parse_float(self.fit_xmin_var.get(), "Fit X Min")
            xmax = parse_float(self.fit_xmax_var.get(), "Fit X Max")
        except ValueError as exc:
            messagebox.showerror("Invalid fit option", str(exc))
            return False
        if xmax <= xmin:
            messagebox.showerror("Invalid fit option", "Fit X Max must be greater than Fit X Min.")
            return False
        return True

    def current_fit_hist_path(self):
        range_hist = getattr(self, "range_hist_path", None)
        if range_hist and range_hist.is_file():
            return range_hist
        full_hist = getattr(self, "full_hist_path", None)
        if full_hist and full_hist.is_file():
            return full_hist
        return None

    def run_fit(self):
        if self.process is not None:
            return
        if not self.validate_fit_inputs():
            return

        hist_path = self.current_fit_hist_path()
        fit_prefix = getattr(self, "fit_output_prefix", None)
        if fit_prefix is None:
            fit_prefix = Path(self.out_dir_var.get()) / self.prefix_var.get().strip()
        cmd = [
            str(FIT_WRAPPER),
            "-i", str(hist_path),
            "-x", format_number(parse_float(self.fit_xmin_var.get(), "Fit X Min")),
            "-X", format_number(parse_float(self.fit_xmax_var.get(), "Fit X Max")),
            "-o", str(fit_prefix),
            "-M", self.fit_model_var.get(),
        ]

        self.fit_result_path = Path(f"{fit_prefix}_compton_edge_fit.png")
        self.append_log("$ " + " ".join(cmd) + "\n\n")
        self.status_var.set("Fitting...")
        self.run_button.configure(state=tk.DISABLED)
        self.fit_button.configure(state=tk.DISABLED)
        self.scan_button.configure(state=tk.DISABLED)
        self.compare_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)

        thread = threading.Thread(target=self._run_fit_subprocess, args=(cmd,), daemon=True)
        thread.start()

    def run_scan(self):
        if self.process is not None:
            return
        if not self.validate_scan_inputs():
            return

        hist_path = self.current_fit_hist_path()
        fit_prefix = getattr(self, "fit_output_prefix", None)
        if fit_prefix is None:
            fit_prefix = Path(self.out_dir_var.get()) / self.prefix_var.get().strip()
        cmd = [
            str(SCAN_WRAPPER),
            "-i", str(hist_path),
            "-x", format_number(parse_float(self.fit_xmin_var.get(), "Fit X Min")),
            "-X", format_number(parse_float(self.fit_xmax_var.get(), "Fit X Max")),
            "-o", str(fit_prefix),
            "-M", self.fit_model_var.get(),
        ]

        self.scan_result_path = Path(f"{fit_prefix}_compton_edge_scan.png")
        self.append_log("$ " + " ".join(cmd) + "\n\n")
        self.status_var.set("Scanning fit range...")
        self.run_button.configure(state=tk.DISABLED)
        self.fit_button.configure(state=tk.DISABLED)
        self.scan_button.configure(state=tk.DISABLED)
        self.compare_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)

        thread = threading.Thread(target=self._run_scan_subprocess, args=(cmd,), daemon=True)
        thread.start()

    def run_model_compare(self):
        if self.process is not None:
            return
        if not self.validate_compare_inputs():
            return

        hist_path = self.current_fit_hist_path()
        fit_prefix = getattr(self, "fit_output_prefix", None)
        if fit_prefix is None:
            fit_prefix = Path(self.out_dir_var.get()) / self.prefix_var.get().strip()
        cmd = [
            str(COMPARE_WRAPPER),
            "-i", str(hist_path),
            "-x", format_number(parse_float(self.fit_xmin_var.get(), "Fit X Min")),
            "-X", format_number(parse_float(self.fit_xmax_var.get(), "Fit X Max")),
            "-o", str(fit_prefix),
        ]

        self.model_compare_paths = [
            Path(f"{fit_prefix}_erfc_linear_compton_edge_fit.png"),
            Path(f"{fit_prefix}_erfc_gaussian_compton_edge_fit.png"),
        ]
        self.append_log("$ " + " ".join(cmd) + "\n\n")
        self.status_var.set("Comparing fit models...")
        self.run_button.configure(state=tk.DISABLED)
        self.fit_button.configure(state=tk.DISABLED)
        self.scan_button.configure(state=tk.DISABLED)
        self.compare_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)

        thread = threading.Thread(target=self._run_model_compare_subprocess, args=(cmd,), daemon=True)
        thread.start()

    def range_prefix(self, base_prefix):
        xmin = self.xmin_var.get().strip().replace(".", "p").replace("-", "m")
        xmax = self.xmax_var.get().strip().replace(".", "p").replace("-", "m")
        return f"{base_prefix}_range_{xmin}_{xmax}"

    def _run_subprocess(self, commands):
        final_code = 0
        try:
            for label, cmd in commands:
                self.after(0, self.append_log, f"\n--- Running {label} ---\n")
                self.process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert self.process.stdout is not None
                for line in self.process.stdout:
                    self.after(0, self.append_log, line)
                final_code = self.process.wait()
                if final_code != 0:
                    break
            self.after(0, self.finish_run, final_code)
        except Exception as exc:
            self.after(0, self.append_log, f"\nError: {exc}\n")
            self.after(0, self.finish_run, 1)

    def _run_fit_subprocess(self, cmd):
        final_code = 0
        output = []
        try:
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert self.process.stdout is not None
            for line in self.process.stdout:
                output.append(line)
                self.after(0, self.append_log, line)
            final_code = self.process.wait()
            self.after(0, self.finish_fit, final_code, "".join(output))
        except Exception as exc:
            self.after(0, self.append_log, f"\nError: {exc}\n")
            self.after(0, self.finish_fit, 1, "")

    def _run_scan_subprocess(self, cmd):
        final_code = 0
        output = []
        try:
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert self.process.stdout is not None
            for line in self.process.stdout:
                output.append(line)
                self.after(0, self.append_log, line)
            final_code = self.process.wait()
            self.after(0, self.finish_scan, final_code, "".join(output))
        except Exception as exc:
            self.after(0, self.append_log, f"\nError: {exc}\n")
            self.after(0, self.finish_scan, 1, "")

    def _run_model_compare_subprocess(self, cmd):
        final_code = 0
        output = []
        try:
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert self.process.stdout is not None
            for line in self.process.stdout:
                output.append(line)
                self.after(0, self.append_log, line)
            final_code = self.process.wait()
            self.after(0, self.finish_model_compare, final_code, "".join(output))
        except Exception as exc:
            self.after(0, self.append_log, f"\nError: {exc}\n")
            self.after(0, self.finish_model_compare, 1, "")

    def stop_analysis(self):
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            self.append_log("\nTerminated by user.\n")

    def finish_run(self, code):
        self.process = None
        self.run_button.configure(state=tk.NORMAL)
        self.stop_button.configure(state=tk.DISABLED)
        if code == 0:
            self.status_var.set("Histograms done")
            self.append_log("\nHistogram drawing done.\n")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.show_result_previews()
        else:
            self.status_var.set(f"Failed with exit code {code}")
            self.append_log(f"\nFailed with exit code {code}.\n")

    def finish_fit(self, code, output):
        self.process = None
        self.run_button.configure(state=tk.NORMAL)
        self.stop_button.configure(state=tk.DISABLED)
        if code == 0:
            self.status_var.set("Fit done")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.fit_result_var.set(self.format_fit_result(output))
            self.append_log("\nFit done.\n")
            self.show_result_previews()
        else:
            self.status_var.set(f"Fit failed with exit code {code}")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.fit_result_var.set(f"Fit failed with exit code {code}.")
            self.append_log(f"\nFit failed with exit code {code}.\n")

    def finish_scan(self, code, output):
        self.process = None
        self.run_button.configure(state=tk.NORMAL)
        self.stop_button.configure(state=tk.DISABLED)
        if code == 0:
            self.status_var.set("Range scan done")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.fit_result_var.set(self.format_scan_result(output))
            self.append_log("\nRange scan done.\n")
            self.show_result_previews()
        else:
            self.status_var.set(f"Range scan failed with exit code {code}")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.fit_result_var.set(f"Range scan failed with exit code {code}.")
            self.append_log(f"\nRange scan failed with exit code {code}.\n")

    def finish_model_compare(self, code, output):
        self.process = None
        self.run_button.configure(state=tk.NORMAL)
        self.stop_button.configure(state=tk.DISABLED)
        if code == 0:
            self.status_var.set("Model comparison done")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.fit_result_var.set(self.format_model_compare_result(output))
            self.append_log("\nModel comparison done.\n")
            self.show_result_previews()
        else:
            self.status_var.set(f"Model comparison failed with exit code {code}")
            self.fit_button.configure(state=tk.NORMAL)
            self.scan_button.configure(state=tk.NORMAL)
            self.compare_button.configure(state=tk.NORMAL)
            self.fit_result_var.set(f"Model comparison failed with exit code {code}.")
            self.append_log(f"\nModel comparison failed with exit code {code}.\n")

    def format_fit_result(self, output):
        values = {}
        for line in output.splitlines():
            if not line.startswith("FIT_RESULT "):
                continue
            for item in line.split()[1:]:
                if "=" not in item:
                    continue
                key, value = item.split("=", 1)
                values[key] = value
        if not values:
            return "Fit completed. See log for details."
        warning_text = values.get("warnings", "none")
        warning_suffix = "" if warning_text == "none" else f"; warnings = {warning_text}"
        return (
            f"{values.get('model', 'fit')}: "
            f"Edge = {values.get('edge', '?')} +/- {values.get('edge_error', '?')} NPE; "
            f"slope seed = {values.get('slope_edge', '?')} NPE; "
            f"sigma = {values.get('sigma', '?')} +/- {values.get('sigma_error', '?')} NPE; "
            f"chi2/ndf = {values.get('chi2_ndf', '?')}; "
            f"pull RMS = {values.get('pull_rms', '?')}"
            f"{warning_suffix}"
        )

    def format_scan_result(self, output):
        values = {}
        for line in output.splitlines():
            if not line.startswith("SCAN_RESULT "):
                continue
            for item in line.split()[1:]:
                if "=" not in item:
                    continue
                key, value = item.split("=", 1)
                values[key] = value
        if not values:
            return "Range scan completed. See log for details."
        warning_text = values.get("warnings", "none")
        warning_suffix = "" if warning_text == "none" else f"; warnings = {warning_text}"
        return (
            f"Range scan {values.get('model', 'fit')}: "
            f"edge mean = {values.get('edge_mean', '?')} NPE; "
            f"range RMS = {values.get('edge_rms', '?')} NPE; "
            f"min/max = {values.get('edge_min', '?')} / {values.get('edge_max', '?')} NPE; "
            f"good fits = {values.get('n_good', '?')} / {values.get('n_total', '?')}"
            f"{warning_suffix}"
        )

    def format_model_compare_result(self, output):
        values = {}
        for line in output.splitlines():
            if not line.startswith("MODEL_COMPARE_RESULT "):
                continue
            for item in line.split()[1:]:
                if "=" not in item:
                    continue
                key, value = item.split("=", 1)
                values[key] = value
        if not values:
            return "Model comparison completed. See log for details."
        warning_parts = []
        if values.get("linear_warnings", "none") != "none":
            warning_parts.append(f"linear={values.get('linear_warnings')}")
        if values.get("gaussian_warnings", "none") != "none":
            warning_parts.append(f"gaussian={values.get('gaussian_warnings')}")
        warning_suffix = "" if not warning_parts else "; warnings = " + ", ".join(warning_parts)
        return (
            f"Model compare: "
            f"linear edge = {values.get('linear_edge', '?')} +/- {values.get('linear_error', '?')} NPE; "
            f"gaussian edge = {values.get('gaussian_edge', '?')} +/- {values.get('gaussian_error', '?')} NPE; "
            f"|diff| = {values.get('edge_diff_abs', '?')} NPE; "
            f"diff pull = {values.get('edge_diff_pull', '?')}"
            f"{warning_suffix}"
        )

    def show_result_previews(self):
        for child in self.preview_frame.winfo_children():
            child.destroy()
        self.preview_images = []

        labels = ["Full range", "Selected range"]
        preview_paths = [(label, path) for label, path in zip(labels, getattr(self, "result_paths", []))]
        if self.fit_result_path and self.fit_result_path.exists():
            preview_paths.append(("Compton edge fit", self.fit_result_path))
        if self.scan_result_path and self.scan_result_path.exists():
            preview_paths.append(("Fit range scan", self.scan_result_path))
        for path in getattr(self, "model_compare_paths", []):
            if path.exists():
                label = "Model compare linear" if "erfc_linear" in path.name else "Model compare gaussian"
                preview_paths.append((label, path))
        for col, (label, path) in enumerate(preview_paths):
            self.preview_frame.columnconfigure(col, weight=1)
            pane = ttk.Frame(self.preview_frame)
            pane.grid(row=0, column=col, sticky="nsew", padx=(0, 12))
            pane.columnconfigure(0, weight=1)
            pane.rowconfigure(1, weight=1)

            header = ttk.Frame(pane)
            header.grid(row=0, column=0, sticky="ew")
            header.columnconfigure(0, weight=1)
            ttk.Label(header, text=f"{label}: {path.name}").grid(row=0, column=0, sticky="w")
            if not path.exists():
                ttk.Label(pane, text="Image not found").grid(row=1, column=0, sticky="nsew")
                continue
            image = tk.PhotoImage(file=str(path))
            factor = max(1, math.ceil(max(image.width() / 520, image.height() / 380)))
            if factor > 1:
                image = image.subsample(factor, factor)
            self.preview_images.append(image)
            ttk.Button(header, text="Open Full Image", command=lambda p=path: self.open_image_window(p)).grid(row=0, column=1, sticky="e", padx=(8, 0))
            ttk.Label(pane, image=image).grid(row=1, column=0, sticky="nsew")

    def open_image_window(self, path):
        if not path.exists():
            messagebox.showerror("Missing image", f"Cannot find:\n{path}")
            return

        top = tk.Toplevel(self)
        top.title(path.name)
        top.geometry("980x720")
        top.minsize(640, 480)

        frame = ttk.Frame(top, padding=8)
        frame.pack(fill=tk.BOTH, expand=True)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        canvas = tk.Canvas(frame, background="white")
        canvas.grid(row=0, column=0, sticky="nsew")
        yscroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=canvas.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=canvas.xview)
        xscroll.grid(row=1, column=0, sticky="ew")
        canvas.configure(xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)

        image = tk.PhotoImage(file=str(path))
        self.preview_images.append(image)
        canvas.create_image(0, 0, anchor="nw", image=image)
        canvas.configure(scrollregion=(0, 0, image.width(), image.height()))

    def append_log(self, text):
        self.log.insert(tk.END, text)
        self.log.see(tk.END)

    def clear_log(self):
        self.log.delete("1.0", tk.END)


def parse_float(value, label):
    try:
        return float(str(value).strip())
    except ValueError as exc:
        raise ValueError(f"{label} must be a number.") from exc


def parse_int(value, label):
    try:
        text = str(value).strip()
        number = int(text)
    except ValueError as exc:
        raise ValueError(f"{label} must be an integer.") from exc
    if str(number) != text:
        raise ValueError(f"{label} must be an integer.")
    return number


def format_number(value):
    return f"{value:.12g}"


if __name__ == "__main__":
    app = ComptonEdgeAnalysisGui()
    app.mainloop()
