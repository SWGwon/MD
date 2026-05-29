#!/usr/bin/env python3
import os
import queue
import shlex
import subprocess
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


HERE = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = Path("/home/sgwon/MD/data/SPE/5inch_set2")
DEFAULT_OUT_DIR = HERE / "set2_spe_results"
MACRO = HERE / "fit_gain_curve_set2.C"


class SpeGainGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("SPE Gain Curve")
        self.geometry("1260x760")
        self.minsize(1020, 640)

        self.log_queue = queue.Queue()
        self.worker = None
        self.preview_image = None

        self.data_dir = tk.StringVar(value=str(DEFAULT_DATA_DIR))
        self.out_dir = tk.StringVar(value=str(DEFAULT_OUT_DIR))
        self.preview_file = tk.StringVar()
        self.channels = tk.StringVar(value="0,1")
        self.v_start = tk.IntVar(value=1600)
        self.v_stop = tk.IntVar(value=2200)
        self.v_step = tk.IntVar(value=100)
        self.bins = tk.IntVar(value=600)
        self.x_min = tk.DoubleVar(value=-1.0)
        self.x_max = tk.StringVar(value="auto")
        self.dynamic_range = tk.DoubleVar(value=2.0)
        self.sampling_ns = tk.DoubleVar(value=2.0)
        self.resistance = tk.DoubleVar(value=50.0)
        self.adc_bits = tk.IntVar(value=14)

        self._build_ui()
        self.after(100, self._drain_log_queue)

    def _build_ui(self):
        root = ttk.Frame(self, padding=12)
        root.pack(fill="both", expand=True)
        root.columnconfigure(0, weight=0)
        root.columnconfigure(1, weight=1)
        root.columnconfigure(2, weight=1)
        root.rowconfigure(2, weight=1)

        paths = ttk.LabelFrame(root, text="Input / Output", padding=10)
        paths.grid(row=0, column=0, columnspan=3, sticky="ew")
        paths.columnconfigure(1, weight=1)

        ttk.Label(paths, text="Data directory").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(paths, textvariable=self.data_dir).grid(row=0, column=1, sticky="ew", pady=3)
        ttk.Button(paths, text="Browse", command=self._browse_data).grid(row=0, column=2, padx=(8, 0), pady=3)

        ttk.Label(paths, text="Output directory").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(paths, textvariable=self.out_dir).grid(row=1, column=1, sticky="ew", pady=3)
        ttk.Button(paths, text="Browse", command=self._browse_out).grid(row=1, column=2, padx=(8, 0), pady=3)

        options = ttk.LabelFrame(root, text="Fit Options", padding=10)
        options.grid(row=1, column=0, sticky="nsew", pady=(10, 10), padx=(0, 10))

        self._labeled_entry(options, "Channels", self.channels, 0)
        self._labeled_spin(options, "Start V", self.v_start, 1, 0, 5000, 100)
        self._labeled_spin(options, "Stop V", self.v_stop, 2, 0, 5000, 100)
        self._labeled_spin(options, "Step V", self.v_step, 3, 1, 1000, 10)
        self._labeled_spin(options, "Bins", self.bins, 4, 50, 5000, 50)
        self._labeled_entry(options, "X min pC", self.x_min, 5)
        self._labeled_entry(options, "X max pC", self.x_max, 6)

        digitizer = ttk.LabelFrame(root, text="Charge Conversion", padding=10)
        digitizer.grid(row=1, column=1, sticky="nsew", pady=(10, 10))
        self._labeled_entry(digitizer, "Dynamic range V", self.dynamic_range, 0)
        self._labeled_entry(digitizer, "Sampling ns", self.sampling_ns, 1)
        self._labeled_entry(digitizer, "Resistance ohm", self.resistance, 2)
        self._labeled_spin(digitizer, "ADC bits", self.adc_bits, 3, 8, 24, 1)

        actions = ttk.Frame(digitizer)
        actions.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(16, 0))
        actions.columnconfigure(0, weight=1)
        actions.columnconfigure(1, weight=1)
        actions.columnconfigure(2, weight=1)
        actions.columnconfigure(3, weight=1)
        self.scan_button = ttk.Button(actions, text="Scan Files", command=self.scan_files)
        self.scan_button.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        self.run_button = ttk.Button(actions, text="Run Fit", command=self.run_fit)
        self.run_button.grid(row=0, column=1, sticky="ew", padx=6)
        self.debug_button = ttk.Button(actions, text="Debug Plots", command=self.make_debug_plots)
        self.debug_button.grid(row=0, column=2, sticky="ew", padx=6)
        ttk.Button(actions, text="Open Output", command=self.open_output).grid(row=0, column=3, sticky="ew", padx=(6, 0))

        preview = ttk.LabelFrame(root, text="Preview", padding=8)
        preview.grid(row=1, column=2, rowspan=2, sticky="nsew", pady=(10, 0))
        preview.rowconfigure(1, weight=1)
        preview.columnconfigure(0, weight=1)

        preview_top = ttk.Frame(preview)
        preview_top.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        preview_top.columnconfigure(0, weight=1)
        self.preview_combo = ttk.Combobox(preview_top, textvariable=self.preview_file, state="readonly")
        self.preview_combo.grid(row=0, column=0, sticky="ew", padx=(0, 6))
        self.preview_combo.bind("<<ComboboxSelected>>", lambda _event: self.load_preview())
        ttk.Button(preview_top, text="Refresh", command=self.refresh_preview_files).grid(row=0, column=1)

        self.preview_label = ttk.Label(preview, anchor="center")
        self.preview_label.grid(row=1, column=0, sticky="nsew")

        log_frame = ttk.LabelFrame(root, text="Log", padding=8)
        log_frame.grid(row=2, column=0, columnspan=2, sticky="nsew")
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

        self.log = tk.Text(log_frame, wrap="word", height=18)
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(log_frame, orient="vertical", command=self.log.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log.configure(yscrollcommand=scroll.set)
        self.refresh_preview_files()

    def _labeled_entry(self, parent, label, variable, row):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(parent, textvariable=variable, width=18).grid(row=row, column=1, sticky="ew", pady=4)
        parent.columnconfigure(1, weight=1)

    def _labeled_spin(self, parent, label, variable, row, from_, to, increment):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Spinbox(parent, textvariable=variable, from_=from_, to=to, increment=increment, width=16).grid(
            row=row, column=1, sticky="ew", pady=4
        )
        parent.columnconfigure(1, weight=1)

    def _browse_data(self):
        selected = filedialog.askdirectory(initialdir=self.data_dir.get())
        if selected:
            self.data_dir.set(selected)

    def _browse_out(self):
        selected = filedialog.askdirectory(initialdir=self.out_dir.get())
        if selected:
            self.out_dir.set(selected)

    def append_log(self, text):
        self.log.insert("end", text)
        self.log.see("end")

    def _drain_log_queue(self):
        try:
            while True:
                self.append_log(self.log_queue.get_nowait())
        except queue.Empty:
            pass
        self.after(100, self._drain_log_queue)

    def voltages(self):
        start = int(self.v_start.get())
        stop = int(self.v_stop.get())
        step = int(self.v_step.get())
        if step <= 0:
            raise ValueError("Voltage step must be positive.")
        if stop < start:
            raise ValueError("Stop voltage must be greater than or equal to start voltage.")
        return list(range(start, stop + 1, step))

    def scan_files(self):
        data_dir = Path(self.data_dir.get()).expanduser()
        self.append_log(f"\nScanning {data_dir}\n")
        missing = 0
        empty = 0
        for voltage in self.voltages():
            for kind in ("LED", "DARK"):
                path = data_dir / f"{kind}_RUN_{voltage}V_prod.root"
                if not path.exists():
                    self.append_log(f"  MISSING {path.name}\n")
                    missing += 1
                elif path.stat().st_size <= 0:
                    self.append_log(f"  EMPTY   {path.name}\n")
                    empty += 1
                else:
                    self.append_log(f"  OK      {path.name} ({path.stat().st_size / 1024 / 1024:.1f} MB)\n")
        self.append_log(f"Scan done: missing={missing}, empty={empty}\n")

    def root_command(self):
        x_max_text = self.x_max.get().strip().lower()
        x_max = -1.0 if x_max_text in ("", "auto") else float(x_max_text)
        args = [
            str(Path(self.data_dir.get()).expanduser()),
            str(Path(self.out_dir.get()).expanduser()),
            self.channels.get(),
            int(self.v_start.get()),
            int(self.v_stop.get()),
            int(self.v_step.get()),
            int(self.bins.get()),
            float(self.x_min.get()),
            x_max,
            float(self.dynamic_range.get()),
            float(self.sampling_ns.get()),
            float(self.resistance.get()),
            int(self.adc_bits.get()),
        ]
        macro_args = ",".join(self._root_literal(arg) for arg in args)
        return ["root", "-l", "-b", "-q", f"{MACRO}({macro_args})"]

    def debug_root_command(self):
        x_max_text = self.x_max.get().strip().lower()
        x_max = 45.0 if x_max_text in ("", "auto") else float(x_max_text)
        args = [
            str(Path(self.data_dir.get()).expanduser()),
            str(Path(self.out_dir.get()).expanduser()),
            self.channels.get(),
            int(self.v_start.get()),
            int(self.v_stop.get()),
            int(self.v_step.get()),
            int(self.bins.get()),
            float(self.x_min.get()),
            x_max,
            float(self.dynamic_range.get()),
            float(self.sampling_ns.get()),
            float(self.resistance.get()),
            int(self.adc_bits.get()),
        ]
        macro_args = ",".join(self._root_literal(arg) for arg in args)
        return [
            "root",
            "-l",
            "-b",
            "-q",
            "-e",
            f".L {MACRO}",
            "-e",
            f"make_set2_debug_charge_plots({macro_args})",
        ]

    @staticmethod
    def _root_literal(value):
        if isinstance(value, str):
            return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
        return str(value)

    def run_fit(self):
        if self.worker and self.worker.is_alive():
            messagebox.showinfo("Running", "A fit is already running.")
            return
        if not MACRO.exists():
            messagebox.showerror("Missing macro", f"Cannot find {MACRO}")
            return

        try:
            cmd = self.root_command()
        except Exception as exc:
            messagebox.showerror("Invalid option", str(exc))
            return

        Path(self.out_dir.get()).expanduser().mkdir(parents=True, exist_ok=True)
        self.append_log("\nRunning:\n  " + " ".join(shlex.quote(part) for part in cmd) + "\n")
        self.run_button.configure(state="disabled")
        self.worker = threading.Thread(
            target=self._run_subprocess,
            args=(cmd, [self.run_button], self.refresh_preview_files),
            daemon=True,
        )
        self.worker.start()

    def make_debug_plots(self):
        if self.worker and self.worker.is_alive():
            messagebox.showinfo("Running", "A ROOT job is already running.")
            return
        if not MACRO.exists():
            messagebox.showerror("Missing macro", f"Cannot find {MACRO}")
            return

        try:
            cmd = self.debug_root_command()
        except Exception as exc:
            messagebox.showerror("Invalid option", str(exc))
            return

        Path(self.out_dir.get()).expanduser().mkdir(parents=True, exist_ok=True)
        self.append_log("\nMaking debug plots:\n  " + " ".join(shlex.quote(part) for part in cmd) + "\n")
        self.debug_button.configure(state="disabled")
        self.worker = threading.Thread(
            target=self._run_subprocess,
            args=(cmd, [self.debug_button], self.refresh_preview_files),
            daemon=True,
        )
        self.worker.start()

    def _run_subprocess(self, cmd, buttons, after_callback=None):
        try:
            proc = subprocess.Popen(
                cmd,
                cwd=str(HERE),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            assert proc.stdout is not None
            for line in proc.stdout:
                self.log_queue.put(line)
            code = proc.wait()
            self.log_queue.put(f"\nROOT exited with code {code}\n")
        except FileNotFoundError:
            self.log_queue.put("\nCannot find root command. Source your ROOT environment first.\n")
        except Exception as exc:
            self.log_queue.put(f"\nRun failed: {exc}\n")
        finally:
            def finish():
                for button in buttons:
                    button.configure(state="normal")
                if after_callback:
                    after_callback()
            self.after(0, finish)

    def refresh_preview_files(self):
        out_dir = Path(self.out_dir.get()).expanduser()
        names = []
        if out_dir.exists():
            patterns = [
                "debug_raw_led_ch*.png",
                "debug_raw_dark_ch*.png",
                "debug_led_minus_dark_ch*.png",
                "gain_curve_set2.png",
                "fit_ch*_montage.png",
            ]
            for pattern in patterns:
                names.extend(path.name for path in sorted(out_dir.glob(pattern)))
        self.preview_combo.configure(values=names)
        if names and self.preview_file.get() not in names:
            self.preview_file.set(names[0])
            self.load_preview()

    def load_preview(self):
        name = self.preview_file.get()
        if not name:
            return
        path = Path(self.out_dir.get()).expanduser() / name
        if not path.exists():
            return
        try:
            image = tk.PhotoImage(file=str(path))
            max_w = max(self.preview_label.winfo_width(), 420)
            max_h = max(self.preview_label.winfo_height(), 360)
            factor = max(
                1,
                int(max((image.width() + max_w - 1) // max_w,
                        (image.height() + max_h - 1) // max_h)),
            )
            if factor > 1:
                image = image.subsample(factor, factor)
            self.preview_image = image
            self.preview_label.configure(image=self.preview_image, text="")
        except Exception as exc:
            self.preview_image = None
            self.preview_label.configure(image="", text=f"Cannot preview {name}\n{exc}")

    def open_output(self):
        out_dir = Path(self.out_dir.get()).expanduser()
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            subprocess.Popen(["xdg-open", str(out_dir)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception as exc:
            messagebox.showerror("Open failed", str(exc))


if __name__ == "__main__":
    app = SpeGainGui()
    app.mainloop()
