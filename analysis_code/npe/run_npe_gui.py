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
WRAPPER = APP_DIR / "run_npe_analysis.sh"


class NpeAnalysisGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("NPE Analysis")
        self.geometry("1180x820")
        self.minsize(1040, 740)
        self.process = None
        self.preview_images = []

        self.source_var = tk.StringVar()
        self.bg_var = tk.StringVar()
        self.out_dir_var = tk.StringVar(value=str(DEFAULT_OUT_DIR))
        self.prefix_var = tk.StringVar(value="npe")
        self.source_label_var = tk.StringVar()
        self.gain_var = tk.StringVar(value="1.0e7")
        self.quantile_var = tk.StringVar(value="1.0")
        self.bins_var = tk.StringVar(value="400")
        self.xmin_var = tk.StringVar(value="0")
        self.xmax_var = tk.StringVar(value="-1")
        self.clock_var = tk.StringVar(value="125.0e6")
        self.channel_vars = [tk.BooleanVar(value=True) for _ in range(8)]

        self._build()

    def _build(self):
        root = ttk.Frame(self, padding=12)
        root.pack(fill=tk.BOTH, expand=True)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(8, weight=1)
        root.rowconfigure(10, weight=1)

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

        ttk.Label(options, text="Run Analysis creates full-range and selected-range plots. Inspect full range first, then adjust the selected range.").grid(row=3, column=0, columnspan=6, sticky="w", pady=(6, 0))

        buttons = ttk.Frame(root)
        buttons.grid(row=6, column=0, columnspan=3, sticky="ew", pady=(8, 8))
        buttons.columnconfigure(2, weight=1)

        self.run_button = ttk.Button(buttons, text="Run Analysis", command=self.run_analysis)
        self.run_button.grid(row=0, column=0, padx=(0, 8))
        self.stop_button = ttk.Button(buttons, text="Stop", command=self.stop_analysis, state=tk.DISABLED)
        self.stop_button.grid(row=0, column=1)
        ttk.Button(buttons, text="Clear Log", command=self.clear_log).grid(row=0, column=3)

        ttk.Label(root, text="Result Preview").grid(row=7, column=0, sticky="w")
        self.preview_frame = ttk.Frame(root)
        self.preview_frame.grid(row=8, column=0, columnspan=3, sticky="nsew")
        self.preview_frame.columnconfigure(0, weight=1)
        self.preview_frame.columnconfigure(1, weight=1)
        self.preview_frame.rowconfigure(1, weight=1)

        ttk.Label(root, text="Log").grid(row=9, column=0, sticky="w", pady=(8, 0))
        log_frame = ttk.Frame(root)
        log_frame.grid(row=10, column=0, columnspan=3, sticky="nsew")
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

        self.log = tk.Text(log_frame, wrap="word", height=16)
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log.configure(yscrollcommand=scroll.set)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(root, textvariable=self.status_var).grid(row=11, column=0, columnspan=3, sticky="ew", pady=(8, 0))

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
            if self.prefix_var.get() in {"", "npe"}:
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

        self.append_log("$ " + " ".join(full_cmd) + "\n")
        self.append_log("$ " + " ".join(range_cmd) + "\n\n")
        self.status_var.set("Running...")
        self.run_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)

        thread = threading.Thread(target=self._run_subprocess, args=(commands,), daemon=True)
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

    def stop_analysis(self):
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            self.append_log("\nTerminated by user.\n")

    def finish_run(self, code):
        self.process = None
        self.run_button.configure(state=tk.NORMAL)
        self.stop_button.configure(state=tk.DISABLED)
        if code == 0:
            self.status_var.set("Done")
            self.append_log("\nDone.\n")
            self.show_result_previews()
        else:
            self.status_var.set(f"Failed with exit code {code}")
            self.append_log(f"\nFailed with exit code {code}.\n")

    def show_result_previews(self):
        for child in self.preview_frame.winfo_children():
            child.destroy()
        self.preview_images = []

        labels = ["Full range", "Selected range"]
        for col, path in enumerate(getattr(self, "result_paths", [])):
            pane = ttk.Frame(self.preview_frame)
            pane.grid(row=0, column=col, sticky="nsew", padx=(0, 12))
            pane.columnconfigure(0, weight=1)
            pane.rowconfigure(1, weight=1)

            header = ttk.Frame(pane)
            header.grid(row=0, column=0, sticky="ew")
            header.columnconfigure(0, weight=1)
            ttk.Label(header, text=f"{labels[col]}: {path.name}").grid(row=0, column=0, sticky="w")
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
    app = NpeAnalysisGui()
    app.mainloop()
