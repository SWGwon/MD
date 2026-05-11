#!/usr/bin/env python3
from html import escape
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs
import subprocess
import sys


APP_DIR = Path(__file__).resolve().parent
WRAPPER = APP_DIR / "run_npe_analysis.sh"
DEFAULT_DATA_DIR = Path.cwd()
DEFAULT_OUT_DIR = Path.cwd()


def root_files(data_dir):
    try:
        return sorted(p.name for p in Path(data_dir).expanduser().glob("*.root") if p.is_file())
    except OSError:
        return []


def page(data=None, log=""):
    data = data or {}
    data_dir = data.get("data_dir", str(DEFAULT_DATA_DIR))
    out_dir = data.get("out_dir", str(DEFAULT_OUT_DIR))
    source = data.get("source", "")
    bg = data.get("background", "")
    files = root_files(data_dir)

    def val(name, default=""):
        return escape(str(data.get(name, default)), quote=True)

    options = "\n".join(
        f'<option value="{escape(name, quote=True)}" {"selected" if name == source else ""}>{escape(name)}</option>'
        for name in files
    )
    bg_options = "\n".join(
        f'<option value="{escape(name, quote=True)}" {"selected" if name == bg else ""}>{escape(name)}</option>'
        for name in files
    )

    return f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>NPE Analysis</title>
  <style>
    body {{ font-family: sans-serif; margin: 24px; max-width: 980px; }}
    label {{ display: block; font-weight: 600; margin: 12px 0 4px; }}
    input, select {{ width: 100%; padding: 7px; box-sizing: border-box; }}
    .grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; }}
    .actions {{ margin-top: 18px; display: flex; gap: 10px; }}
    button {{ padding: 9px 14px; cursor: pointer; }}
    pre {{ background: #111; color: #eee; padding: 14px; overflow: auto; min-height: 180px; }}
    .hint {{ color: #555; font-size: 0.92em; }}
  </style>
</head>
<body>
  <h1>NPE Analysis</h1>
  <form method="post" action="/run">
    <label>Data directory</label>
    <input name="data_dir" value="{escape(data_dir, quote=True)}">
    <div class="hint">Enter the directory containing source/background ROOT files, then press Refresh file list.</div>

    <div class="actions">
      <button formaction="/" formmethod="post">Refresh file list</button>
    </div>

    <label>Source ROOT</label>
    <select name="source">{options}</select>
    <label>Background ROOT</label>
    <select name="background">{bg_options}</select>

    <label>Output directory</label>
    <input name="out_dir" value="{escape(out_dir, quote=True)}">

    <label>Output prefix</label>
    <input name="prefix" value="{val("prefix", "nocollimator")}">

    <div class="grid">
      <div>
        <label>Gain</label>
        <input name="gain" value="{val("gain", "1.0e7")}">
      </div>
      <div>
        <label>X quantile</label>
        <input name="quantile" value="{val("quantile", "1.0")}">
      </div>
      <div>
        <label>Bins</label>
        <input name="bins" value="{val("bins", "400")}">
      </div>
      <div>
        <label>X min</label>
        <input name="xmin" value="{val("xmin", "0")}">
      </div>
      <div>
        <label>X max</label>
        <input name="xmax" value="{val("xmax", "-1")}">
      </div>
      <div>
        <label>TTT clock Hz</label>
        <input name="clock" value="{val("clock", "125.0e6")}">
      </div>
    </div>

    <div class="actions">
      <button type="submit">Run Analysis</button>
    </div>
  </form>

  <h2>Log</h2>
  <pre>{escape(log)}</pre>
</body>
</html>"""


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.respond(page())

    def do_POST(self):
        length = int(self.headers.get("content-length", "0"))
        fields = parse_qs(self.rfile.read(length).decode("utf-8"))
        data = {key: values[0] for key, values in fields.items() if values}

        if self.path == "/run":
            log = self.run_analysis(data)
            self.respond(page(data, log))
        else:
            self.respond(page(data))

    def run_analysis(self, data):
        data_dir = Path(data.get("data_dir", ".")).expanduser().resolve()
        out_dir = Path(data.get("out_dir", ".")).expanduser().resolve()
        source = data_dir / data.get("source", "")
        bg = data_dir / data.get("background", "")
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(WRAPPER),
            "-d", str(data_dir),
            "-s", str(source),
            "-b", str(bg),
            "-O", str(out_dir),
            "-o", data.get("prefix", "npe"),
            "-g", data.get("gain", "1.0e7"),
            "-q", data.get("quantile", "1.0"),
            "-n", data.get("bins", "400"),
            "-x", data.get("xmin", "0"),
            "-X", data.get("xmax", "-1"),
            "-c", data.get("clock", "125.0e6"),
        ]

        try:
            result = subprocess.run(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            return "$ " + " ".join(cmd) + "\n\n" + result.stdout
        except Exception as exc:
            return "$ " + " ".join(cmd) + f"\n\nError: {exc}\n"

    def respond(self, body):
        encoded = body.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)

    def log_message(self, fmt, *args):
        sys.stderr.write(fmt % args + "\n")


def main():
    host = "127.0.0.1"
    port = 8765
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"Open this URL in your browser: http://{host}:{port}", flush=True)
    print("Press Ctrl+C to stop.", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
