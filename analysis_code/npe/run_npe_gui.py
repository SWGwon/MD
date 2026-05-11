#!/usr/bin/env python3
import subprocess
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


APP_DIR = Path(__file__).resolve().parent
WRAPPER = APP_DIR / "run_npe_analysis.sh"


class NpeAnalysisGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("NPE Analysis")
        self.geometry("860x620")
        self.minsize(760, 540)
        self.process = None

        self.source_var = tk.StringVar()
        self.bg_var = tk.StringVar()
        self.out_dir_var = tk.StringVar(value=str(Path.cwd()))
        self.prefix_var = tk.StringVar(value="nocollimator")
        self.gain_var = tk.StringVar(value="1.0e7")
        self.quantile_var = tk.StringVar(value="1.0")
        self.bins_var = tk.StringVar(value="400")
        self.xmin_var = tk.StringVar(value="0")
        self.xmax_var = tk.StringVar(value="-1")
        self.clock_var = tk.StringVar(value="125.0e6")

        self._build()

    def _build(self):
        root = ttk.Frame(self, padding=12)
        root.pack(fill=tk.BOTH, expand=True)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(7, weight=1)

        self._file_row(root, 0, "Source ROOT", self.source_var, self.pick_source)
        self._file_row(root, 1, "Background ROOT", self.bg_var, self.pick_background)
        self._file_row(root, 2, "Output Dir", self.out_dir_var, self.pick_output_dir)

        ttk.Label(root, text="Output Prefix").grid(row=3, column=0, sticky="w", pady=4)
        ttk.Entry(root, textvariable=self.prefix_var).grid(row=3, column=1, columnspan=2, sticky="ew", pady=4)

        options = ttk.LabelFrame(root, text="Options", padding=8)
        options.grid(row=4, column=0, columnspan=3, sticky="ew", pady=(8, 4))
        for col in range(6):
            options.columnconfigure(col, weight=1)

        self._option(options, 0, 0, "Gain", self.gain_var)
        self._option(options, 0, 2, "X Quantile", self.quantile_var)
        self._option(options, 0, 4, "Bins", self.bins_var)
        self._option(options, 1, 0, "X Min", self.xmin_var)
        self._option(options, 1, 2, "X Max", self.xmax_var)
        self._option(options, 1, 4, "TTT Clock Hz", self.clock_var)

        buttons = ttk.Frame(root)
        buttons.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(8, 8))
        buttons.columnconfigure(2, weight=1)

        self.run_button = ttk.Button(buttons, text="Run Analysis", command=self.run_analysis)
        self.run_button.grid(row=0, column=0, padx=(0, 8))
        self.stop_button = ttk.Button(buttons, text="Stop", command=self.stop_analysis, state=tk.DISABLED)
        self.stop_button.grid(row=0, column=1)
        ttk.Button(buttons, text="Clear Log", command=self.clear_log).grid(row=0, column=3)

        ttk.Label(root, text="Log").grid(row=6, column=0, sticky="w")
        log_frame = ttk.Frame(root)
        log_frame.grid(row=7, column=0, columnspan=3, sticky="nsew")
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

        self.log = tk.Text(log_frame, wrap="word", height=16)
        self.log.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log.configure(yscrollcommand=scroll.set)

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(root, textvariable=self.status_var).grid(row=8, column=0, columnspan=3, sticky="ew", pady=(8, 0))

    def _file_row(self, parent, row, label, variable, command):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4, padx=(8, 8))
        ttk.Button(parent, text="Browse", command=command).grid(row=row, column=2, sticky="ew", pady=4)

    def _option(self, parent, row, col, label, variable):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=(0, 4), pady=4)
        ttk.Entry(parent, textvariable=variable, width=14).grid(row=row, column=col + 1, sticky="ew", padx=(0, 12), pady=4)

    def pick_source(self):
        path = filedialog.askopenfilename(title="Select source ROOT file", filetypes=[("ROOT files", "*.root"), ("All files", "*")])
        if path:
            self.source_var.set(path)
            if self.prefix_var.get() in {"", "nocollimator"}:
                self.prefix_var.set(Path(path).stem)

    def pick_background(self):
        path = filedialog.askopenfilename(title="Select background ROOT file", filetypes=[("ROOT files", "*.root"), ("All files", "*")])
        if path:
            self.bg_var.set(path)

    def pick_output_dir(self):
        path = filedialog.askdirectory(title="Select output directory")
        if path:
            self.out_dir_var.set(path)

    def validate_inputs(self):
        if not WRAPPER.exists():
            messagebox.showerror("Missing wrapper", f"Cannot find:\n{WRAPPER}")
            return False
        if not self.source_var.get() or not Path(self.source_var.get()).is_file():
            messagebox.showerror("Missing source", "Select a valid source ROOT file.")
            return False
        if not self.bg_var.get() or not Path(self.bg_var.get()).is_file():
            messagebox.showerror("Missing background", "Select a valid background ROOT file.")
            return False
        if not self.prefix_var.get().strip():
            messagebox.showerror("Missing prefix", "Enter an output prefix.")
            return False
        Path(self.out_dir_var.get()).mkdir(parents=True, exist_ok=True)
        return True

    def run_analysis(self):
        if self.process is not None:
            return
        if not self.validate_inputs():
            return

        cmd = [
            str(WRAPPER),
            "-s", self.source_var.get(),
            "-b", self.bg_var.get(),
            "-O", self.out_dir_var.get(),
            "-o", self.prefix_var.get().strip(),
            "-g", self.gain_var.get().strip(),
            "-q", self.quantile_var.get().strip(),
            "-n", self.bins_var.get().strip(),
            "-x", self.xmin_var.get().strip(),
            "-X", self.xmax_var.get().strip(),
            "-c", self.clock_var.get().strip(),
        ]

        self.append_log("$ " + " ".join(cmd) + "\n\n")
        self.status_var.set("Running...")
        self.run_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)

        thread = threading.Thread(target=self._run_subprocess, args=(cmd,), daemon=True)
        thread.start()

    def _run_subprocess(self, cmd):
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
                self.after(0, self.append_log, line)
            code = self.process.wait()
            self.after(0, self.finish_run, code)
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
        else:
            self.status_var.set(f"Failed with exit code {code}")
            self.append_log(f"\nFailed with exit code {code}.\n")

    def append_log(self, text):
        self.log.insert(tk.END, text)
        self.log.see(tk.END)

    def clear_log(self):
        self.log.delete("1.0", tk.END)


if __name__ == "__main__":
    app = NpeAnalysisGui()
    app.mainloop()
