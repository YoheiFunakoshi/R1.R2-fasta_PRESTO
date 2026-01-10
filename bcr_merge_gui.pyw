#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
import queue
import shutil
import struct
import subprocess
import threading
import zlib
import sys

import tkinter as tk
from tkinter import filedialog, messagebox, ttk


DEFAULT_ADAPTER_5P_R1 = "TGAGTTCCACGACACCGTCA"
DEFAULT_ADAPTER_5P_R2 = "CTAATACGACTCCGAATTCC"
DEFAULT_ADAPTER_3P_R1 = "GGAATTCGGAGTCGTATTAG"
DEFAULT_ADAPTER_3P_R2 = "TGACGGTGTCGTGGAACTCA"

MODE_LABELS = {
    "merge-only": "withoutTrim",
    "trim-merge": "trimThenMerge",
    "merge-trim": "mergeThenTrim",
}


def infer_prefix(r1_path: Path) -> str:
    name = r1_path.name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for marker in ["_R1_001", "_R1", "-R1", ".R1"]:
        if name.endswith(marker):
            return name[: -len(marker)]
    if "_R1_" in name:
        return name.split("_R1_")[0]
    return name


def mode_label(mode: str) -> str:
    return MODE_LABELS.get(mode, mode)


def png_chunk(tag: bytes, data: bytes) -> bytes:
    return struct.pack("!I", len(data)) + tag + data + struct.pack(
        "!I", zlib.crc32(tag + data) & 0xFFFFFFFF
    )


def write_icon_png(path: Path, size: int = 64) -> None:
    width = height = size
    bg = (18, 40, 60)
    border = (235, 240, 245)
    accent = (0, 180, 160)
    accent2 = (60, 130, 210)

    rows = []
    for y in range(height):
        row = bytearray()
        row.append(0)
        for x in range(width):
            color = bg
            if x < 3 or x >= width - 3 or y < 3 or y >= height - 3:
                color = border
            elif abs(x - y) < 2:
                color = accent
            elif abs((width - 1 - x) - y) < 2:
                color = accent2
            row.extend(color)
        rows.append(bytes(row))

    raw = b"".join(rows)
    compressed = zlib.compress(raw, level=9)

    signature = b"\x89PNG\r\n\x1a\n"
    ihdr = struct.pack("!IIBBBBB", width, height, 8, 2, 0, 0, 0)
    png_bytes = signature + png_chunk(b"IHDR", ihdr) + png_chunk(
        b"IDAT", compressed
    ) + png_chunk(b"IEND", b"")
    path.write_bytes(png_bytes)


def write_icon_ico(ico_path: Path, png_bytes: bytes, size: int = 64) -> None:
    width = size if size < 256 else 0
    height = size if size < 256 else 0
    header = struct.pack("<HHH", 0, 1, 1)
    entry = struct.pack(
        "<BBBBHHII",
        width,
        height,
        0,
        0,
        1,
        32,
        len(png_bytes),
        6 + 16,
    )
    ico_path.write_bytes(header + entry + png_bytes)


def ensure_icon_files(base_dir: Path) -> tuple[Path, Path]:
    icon_png = base_dir / "bcr_merge_icon.png"
    icon_ico = base_dir / "bcr_merge_icon.ico"
    if not icon_png.exists():
        write_icon_png(icon_png)
    if not icon_ico.exists():
        png_bytes = icon_png.read_bytes()
        write_icon_ico(icon_ico, png_bytes, size=64)
    return icon_png, icon_ico


def guess_conda_exe() -> str | None:
    home = Path.home()
    candidates = [
        Path(r"C:\miniforge3\Scripts\conda.exe"),
        Path(r"C:\ProgramData\miniforge3\Scripts\conda.exe"),
        home / "miniforge3" / "Scripts" / "conda.exe",
        home / "miniconda3" / "Scripts" / "conda.exe",
    ]
    for path in candidates:
        if path.exists():
            return str(path)
    found = shutil.which("conda")
    if found:
        return found
    return None


def guess_conda_env() -> str:
    candidates = [
        Path(r"C:\miniforge3\envs\presto_env"),
        Path(r"C:\ProgramData\miniforge3\envs\presto_env"),
        Path.home() / "miniforge3" / "envs" / "presto_env",
    ]
    for path in candidates:
        if path.exists():
            return str(path)
    return "presto_env"


class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("BCR Merge GUI")
        base_dir = Path(__file__).resolve().parent
        icon_png, _ = ensure_icon_files(base_dir)
        try:
            icon_img = tk.PhotoImage(file=str(icon_png))
            self.iconphoto(True, icon_img)
            self._icon_img = icon_img
        except tk.TclError:
            pass

        self.queue: "queue.Queue[tuple[str, str | int]]" = queue.Queue()
        self.worker_thread: threading.Thread | None = None
        self.running = False

        self.var_r1 = tk.StringVar()
        self.var_r2 = tk.StringVar()
        self.var_outdir = tk.StringVar()
        self.var_prefix = tk.StringVar()
        self.var_assemble_prefix = tk.StringVar()
        self.var_mode = tk.StringVar(value="merge-only")
        self.var_threads = tk.StringVar(value="8")
        self.var_error_rate = tk.StringVar(value="0.15")
        self.var_overlap = tk.StringVar(value="10")
        self.var_min_len = tk.StringVar(value="50")
        self.var_adapt_5p_r1 = tk.StringVar(value=DEFAULT_ADAPTER_5P_R1)
        self.var_adapt_5p_r2 = tk.StringVar(value=DEFAULT_ADAPTER_5P_R2)
        self.var_adapt_3p_r1 = tk.StringVar(value=DEFAULT_ADAPTER_3P_R1)
        self.var_adapt_3p_r2 = tk.StringVar(value=DEFAULT_ADAPTER_3P_R2)
        self.var_failed = tk.BooleanVar(value=False)
        self.var_dry_run = tk.BooleanVar(value=False)

        default_use_conda = "presto_env" not in os.environ.get("CONDA_PREFIX", "")
        self.var_use_conda = tk.BooleanVar(value=default_use_conda)
        self.var_conda_exe = tk.StringVar(value=guess_conda_exe() or "")
        self.var_conda_env = tk.StringVar(value=guess_conda_env())

        self._build_ui()
        self.after(100, self._poll_queue)

    def _build_ui(self) -> None:
        frame = ttk.Frame(self, padding=10)
        frame.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        row = 0
        ttk.Label(frame, text="R1 FASTQ").grid(row=row, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.var_r1, width=60).grid(
            row=row, column=1, sticky="ew", padx=5
        )
        ttk.Button(frame, text="Browse", command=self._browse_r1).grid(
            row=row, column=2, sticky="ew"
        )
        row += 1

        ttk.Label(frame, text="R2 FASTQ").grid(row=row, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.var_r2, width=60).grid(
            row=row, column=1, sticky="ew", padx=5
        )
        ttk.Button(frame, text="Browse", command=self._browse_r2).grid(
            row=row, column=2, sticky="ew"
        )
        row += 1

        ttk.Label(frame, text="Output base folder").grid(row=row, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.var_outdir, width=60).grid(
            row=row, column=1, sticky="ew", padx=5
        )
        ttk.Button(frame, text="Browse", command=self._browse_outdir).grid(
            row=row, column=2, sticky="ew"
        )
        row += 1

        ttk.Label(frame, text="Prefix (optional)").grid(row=row, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.var_prefix, width=30).grid(
            row=row, column=1, sticky="w", padx=5
        )
        row += 1

        ttk.Label(frame, text="Assemble prefix (optional)").grid(
            row=row, column=0, sticky="w"
        )
        ttk.Entry(frame, textvariable=self.var_assemble_prefix, width=30).grid(
            row=row, column=1, sticky="w", padx=5
        )
        row += 1

        mode_frame = ttk.LabelFrame(frame, text="Mode")
        mode_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=5)
        ttk.Radiobutton(
            mode_frame, text="merge-only", variable=self.var_mode, value="merge-only"
        ).grid(row=0, column=0, sticky="w", padx=5)
        ttk.Radiobutton(
            mode_frame, text="trim-merge", variable=self.var_mode, value="trim-merge"
        ).grid(row=0, column=1, sticky="w", padx=5)
        ttk.Radiobutton(
            mode_frame, text="merge-trim", variable=self.var_mode, value="merge-trim"
        ).grid(row=0, column=2, sticky="w", padx=5)
        row += 1

        adapter_frame = ttk.LabelFrame(frame, text="Adapters")
        adapter_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=5)
        ttk.Label(adapter_frame, text="R1 5'").grid(row=0, column=0, sticky="w")
        ttk.Entry(adapter_frame, textvariable=self.var_adapt_5p_r1, width=30).grid(
            row=0, column=1, sticky="w", padx=5
        )
        ttk.Label(adapter_frame, text="R2 5'").grid(row=0, column=2, sticky="w")
        ttk.Entry(adapter_frame, textvariable=self.var_adapt_5p_r2, width=30).grid(
            row=0, column=3, sticky="w", padx=5
        )
        ttk.Label(adapter_frame, text="R1 3'").grid(row=1, column=0, sticky="w")
        ttk.Entry(adapter_frame, textvariable=self.var_adapt_3p_r1, width=30).grid(
            row=1, column=1, sticky="w", padx=5
        )
        ttk.Label(adapter_frame, text="R2 3'").grid(row=1, column=2, sticky="w")
        ttk.Entry(adapter_frame, textvariable=self.var_adapt_3p_r2, width=30).grid(
            row=1, column=3, sticky="w", padx=5
        )
        row += 1

        opt_frame = ttk.LabelFrame(frame, text="cutadapt options")
        opt_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=5)
        ttk.Label(opt_frame, text="Threads").grid(row=0, column=0, sticky="w")
        ttk.Entry(opt_frame, textvariable=self.var_threads, width=8).grid(
            row=0, column=1, sticky="w", padx=5
        )
        ttk.Label(opt_frame, text="Error rate (-e)").grid(
            row=0, column=2, sticky="w"
        )
        ttk.Entry(opt_frame, textvariable=self.var_error_rate, width=8).grid(
            row=0, column=3, sticky="w", padx=5
        )
        ttk.Label(opt_frame, text="Overlap (-O)").grid(row=0, column=4, sticky="w")
        ttk.Entry(opt_frame, textvariable=self.var_overlap, width=8).grid(
            row=0, column=5, sticky="w", padx=5
        )
        ttk.Label(opt_frame, text="Min length").grid(row=0, column=6, sticky="w")
        ttk.Entry(opt_frame, textvariable=self.var_min_len, width=8).grid(
            row=0, column=7, sticky="w", padx=5
        )
        row += 1

        exec_frame = ttk.LabelFrame(frame, text="Execution")
        exec_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=5)
        exec_frame.columnconfigure(1, weight=1)
        ttk.Checkbutton(
            exec_frame, text="Run inside conda env", variable=self.var_use_conda
        ).grid(row=0, column=0, sticky="w", padx=5)
        ttk.Label(exec_frame, text="Conda exe").grid(row=1, column=0, sticky="w")
        ttk.Entry(exec_frame, textvariable=self.var_conda_exe, width=70).grid(
            row=1, column=1, sticky="ew", padx=5
        )
        ttk.Button(exec_frame, text="Browse", command=self._browse_conda).grid(
            row=1, column=2, sticky="w"
        )
        ttk.Label(exec_frame, text="Env name or path").grid(row=2, column=0, sticky="w")
        ttk.Entry(exec_frame, textvariable=self.var_conda_env, width=70).grid(
            row=2, column=1, sticky="ew", padx=5
        )
        row += 1

        flags_frame = ttk.Frame(frame)
        flags_frame.grid(row=row, column=0, columnspan=3, sticky="w")
        ttk.Checkbutton(
            flags_frame, text="Emit assemble-fail files", variable=self.var_failed
        ).grid(row=0, column=0, sticky="w", padx=5)
        ttk.Checkbutton(flags_frame, text="Dry run", variable=self.var_dry_run).grid(
            row=0, column=1, sticky="w", padx=5
        )
        row += 1

        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=5)
        self.run_btn = ttk.Button(btn_frame, text="Run", command=self._start_run)
        self.run_btn.grid(row=0, column=0, sticky="w", padx=5)
        ttk.Button(btn_frame, text="Open output folder", command=self._open_outdir).grid(
            row=0, column=1, sticky="w", padx=5
        )
        ttk.Button(btn_frame, text="Quit", command=self.destroy).grid(
            row=0, column=2, sticky="w", padx=5
        )
        row += 1

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(frame, textvariable=self.status_var).grid(
            row=row, column=0, columnspan=3, sticky="w"
        )
        row += 1

        log_frame = ttk.LabelFrame(frame, text="Log")
        log_frame.grid(row=row, column=0, columnspan=3, sticky="nsew", pady=5)
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)
        self.log_text = tk.Text(log_frame, height=12, wrap="none")
        self.log_text.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=scroll.set)

        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(row, weight=1)

    def _browse_r1(self) -> None:
        path = filedialog.askopenfilename(
            title="Select R1 FASTQ",
            filetypes=[
                ("FASTQ", "*.fastq *.fq *.fastq.gz *.fq.gz"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.var_r1.set(path)
            self._auto_fill_from_r1(Path(path))

    def _browse_r2(self) -> None:
        path = filedialog.askopenfilename(
            title="Select R2 FASTQ",
            filetypes=[
                ("FASTQ", "*.fastq *.fq *.fastq.gz *.fq.gz"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.var_r2.set(path)

    def _browse_outdir(self) -> None:
        path = filedialog.askdirectory(title="Select output folder")
        if path:
            self.var_outdir.set(path)

    def _browse_conda(self) -> None:
        path = filedialog.askopenfilename(
            title="Select conda.exe",
            filetypes=[("conda.exe", "conda.exe"), ("All files", "*.*")],
        )
        if path:
            self.var_conda_exe.set(path)

    def _auto_fill_from_r1(self, r1_path: Path) -> None:
        if not self.var_outdir.get():
            self.var_outdir.set(str(r1_path.parent))
        if not self.var_prefix.get():
            self.var_prefix.set(infer_prefix(r1_path))

    def _current_prefix_base(self) -> str:
        prefix = self.var_prefix.get().strip()
        if prefix:
            return prefix
        r1 = self.var_r1.get().strip()
        if r1:
            return infer_prefix(Path(r1))
        return "sample"

    def _current_run_dir(self) -> Path:
        r1 = self.var_r1.get().strip()
        outdir = self.var_outdir.get().strip()
        if outdir:
            base = Path(outdir)
        elif r1:
            base = Path(r1).parent
        else:
            base = Path.cwd()
        run_name = f"{self._current_prefix_base()}_{mode_label(self.var_mode.get())}"
        return base / run_name

    def _append_log(self, text: str) -> None:
        self.log_text.insert("end", text)
        self.log_text.see("end")

    def _build_command(self) -> list[str]:
        script_dir = Path(__file__).resolve().parent
        cli_path = script_dir / "bcr_merge_cli.py"
        if not cli_path.exists():
            raise FileNotFoundError(f"CLI not found: {cli_path}")

        r1 = self.var_r1.get().strip()
        r2 = self.var_r2.get().strip()
        outdir = self.var_outdir.get().strip()
        if not outdir:
            outdir = str(Path(r1).parent)

        cmd = [
            "--mode",
            self.var_mode.get(),
            "--r1",
            r1,
            "--r2",
            r2,
            "--outdir",
            outdir,
            "--threads",
            self.var_threads.get().strip(),
            "--cutadapt-error-rate",
            self.var_error_rate.get().strip(),
            "--cutadapt-overlap",
            self.var_overlap.get().strip(),
            "--cutadapt-min-length",
            self.var_min_len.get().strip(),
            "--adapter-5p-r1",
            self.var_adapt_5p_r1.get().strip(),
            "--adapter-5p-r2",
            self.var_adapt_5p_r2.get().strip(),
            "--adapter-3p-r1",
            self.var_adapt_3p_r1.get().strip(),
            "--adapter-3p-r2",
            self.var_adapt_3p_r2.get().strip(),
        ]

        prefix = self.var_prefix.get().strip()
        if prefix:
            cmd += ["--prefix", prefix]

        assemble_prefix = self.var_assemble_prefix.get().strip()
        if assemble_prefix:
            cmd += ["--assemble-prefix", assemble_prefix]

        if self.var_failed.get():
            cmd.append("--failed")
        if self.var_dry_run.get():
            cmd.append("--dry-run")

        if self.var_use_conda.get():
            conda_exe = self.var_conda_exe.get().strip()
            if not conda_exe:
                raise FileNotFoundError("conda.exe path is empty.")
            env_value = self.var_conda_env.get().strip() or "presto_env"
            env_path = Path(env_value)
            if env_path.exists():
                return [
                    conda_exe,
                    "run",
                    "-p",
                    str(env_path),
                    "python",
                    "-u",
                    str(cli_path),
                ] + cmd
            return [
                conda_exe,
                "run",
                "-n",
                env_value,
                "python",
                "-u",
                str(cli_path),
            ] + cmd

        return [sys.executable, "-u", str(cli_path)] + cmd

    def _start_run(self) -> None:
        if self.running:
            return
        r1 = self.var_r1.get().strip()
        r2 = self.var_r2.get().strip()
        if not r1 or not Path(r1).exists():
            messagebox.showerror("Error", "R1 FASTQ not found.")
            return
        if not r2 or not Path(r2).exists():
            messagebox.showerror("Error", "R2 FASTQ not found.")
            return

        try:
            cmd = self._build_command()
        except Exception as exc:
            messagebox.showerror("Error", str(exc))
            return

        self.running = True
        self.run_btn.configure(state="disabled")
        self.status_var.set("Running...")
        self._append_log("\n>> " + subprocess.list2cmdline(cmd) + "\n")
        self._append_log(f">> output folder: {self._current_run_dir()}\n")

        def worker() -> None:
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert proc.stdout is not None
                for line in proc.stdout:
                    self.queue.put(("line", line))
                rc = proc.wait()
                self.queue.put(("done", str(rc)))
            except Exception as exc:
                self.queue.put(("line", f"[error] {exc}\n"))
                self.queue.put(("done", "error"))

        self.worker_thread = threading.Thread(target=worker, daemon=True)
        self.worker_thread.start()

    def _poll_queue(self) -> None:
        while True:
            try:
                kind, payload = self.queue.get_nowait()
            except queue.Empty:
                break
            if kind == "line":
                self._append_log(str(payload))
            elif kind == "done":
                self.running = False
                self.run_btn.configure(state="normal")
                if payload == "0":
                    self.status_var.set("Completed.")
                else:
                    self.status_var.set(f"Finished with exit code: {payload}")
        self.after(100, self._poll_queue)

    def _open_outdir(self) -> None:
        run_dir = self._current_run_dir()
        if run_dir.exists():
            os.startfile(run_dir)
            return
        base = self.var_outdir.get().strip()
        if not base:
            messagebox.showinfo("Info", "Output folder is empty.")
            return
        path = Path(base)
        if not path.exists():
            messagebox.showerror("Error", f"Output folder not found: {path}")
            return
        os.startfile(path)


def main() -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.parse_args()
    app = App()
    app.minsize(900, 650)
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
