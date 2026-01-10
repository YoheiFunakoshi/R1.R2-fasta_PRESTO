#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys

DEFAULT_ADAPTER_5P_R1 = "TGAGTTCCACGACACCGTCA"
DEFAULT_ADAPTER_5P_R2 = "CTAATACGACTCCGAATTCC"
DEFAULT_ADAPTER_3P_R1 = "GGAATTCGGAGTCGTATTAG"
DEFAULT_ADAPTER_3P_R2 = "TGACGGTGTCGTGGAACTCA"

MODE_LABELS = {
    "merge-only": "withoutTrim",
    "trim-merge": "trimThenMerge",
    "merge-trim": "mergeThenTrim",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge paired-end BCR reads with PRESTO AssemblePairs and optional cutadapt trimming."
        )
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["merge-only", "trim-merge", "merge-trim"],
        help="Workflow to run: merge-only, trim-merge, or merge-trim.",
    )
    parser.add_argument("--r1", required=True, help="R1 FASTQ path.")
    parser.add_argument("--r2", required=True, help="R2 FASTQ path.")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to the R1 directory.",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Sample prefix used for outputs. Defaults to an inferred name from R1.",
    )
    parser.add_argument(
        "--assemble-prefix",
        default=None,
        help="Prefix passed to AssemblePairs --outname. Defaults to <prefix>_presto.",
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Optional run identifier used for output folder naming (default: timestamp).",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="AssemblePairs log filename. Defaults to AP_align.log in outdir.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for cutadapt (-j).",
    )
    parser.add_argument(
        "--cutadapt-error-rate",
        type=float,
        default=0.15,
        help="cutadapt -e error rate.",
    )
    parser.add_argument(
        "--cutadapt-overlap",
        type=int,
        default=10,
        help="cutadapt -O overlap length.",
    )
    parser.add_argument(
        "--cutadapt-min-length",
        type=int,
        default=50,
        help="cutadapt --minimum-length.",
    )
    parser.add_argument(
        "--adapter-5p-r1",
        default=DEFAULT_ADAPTER_5P_R1,
        help="R1 5' adapter sequence (without '^').",
    )
    parser.add_argument(
        "--adapter-5p-r2",
        default=DEFAULT_ADAPTER_5P_R2,
        help="R2 5' adapter sequence (without '^').",
    )
    parser.add_argument(
        "--adapter-3p-r1",
        default=DEFAULT_ADAPTER_3P_R1,
        help="R1 3' adapter sequence (without '$').",
    )
    parser.add_argument(
        "--adapter-3p-r2",
        default=DEFAULT_ADAPTER_3P_R2,
        help="R2 3' adapter sequence (without '$').",
    )
    parser.add_argument(
        "--coord",
        default="illumina",
        help="AssemblePairs --coord (default: illumina).",
    )
    parser.add_argument(
        "--rc",
        default="tail",
        help="AssemblePairs --rc (default: tail).",
    )
    parser.add_argument(
        "--failed",
        action="store_true",
        help="Emit assemble-fail FASTQ files.",
    )
    parser.add_argument(
        "--swap-r1r2",
        action="store_true",
        help="Swap R1/R2 for AssemblePairs (-1 R1, -2 R2).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    return parser.parse_args()


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


def validate_run_id(run_id: str) -> None:
    invalid = '<>:"/\\|?*'
    if any(ch in invalid for ch in run_id):
        raise ValueError(f"--run-id contains invalid characters: {invalid}")


def sanitize_adapter(seq: str) -> str:
    seq = seq.strip()
    return seq.replace("^", "").replace("$", "")


def ensure_prefix(seq: str, prefix: str) -> str:
    seq = seq.strip()
    if seq.startswith(prefix):
        return seq
    return f"{prefix}{seq}"


def ensure_suffix(seq: str, suffix: str) -> str:
    seq = seq.strip()
    if seq.endswith(suffix):
        return seq
    return f"{seq}{suffix}"


def resolve_assemblepairs_path(explicit: str | None = None) -> Path:
    if explicit:
        path = Path(explicit)
        if path.exists():
            return path
        raise FileNotFoundError(f"AssemblePairs.py not found at: {explicit}")

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidate = Path(conda_prefix) / "Scripts" / "AssemblePairs.py"
        if candidate.exists():
            return candidate

    found = shutil.which("AssemblePairs.py")
    if found:
        return Path(found)

    raise FileNotFoundError(
        "AssemblePairs.py not found. Activate presto_env or install presto."
    )


def resolve_cutadapt_cmd() -> list[str]:
    found = shutil.which("cutadapt")
    if found:
        return [found]
    return [sys.executable, "-m", "cutadapt"]


def run_command(cmd: list[str], dry_run: bool, cwd: Path | None = None) -> None:
    cmd_str = subprocess.list2cmdline(cmd)
    print(f">> {cmd_str}")
    if dry_run:
        return
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None)


def fastq_to_fasta(in_path: Path, out_path: Path) -> None:
    opener = gzip.open if in_path.suffix == ".gz" else open
    with opener(in_path, "rt", encoding="utf-8", errors="replace") as fin:
        with open(out_path, "w", encoding="utf-8") as fout:
            while True:
                header = fin.readline()
                if not header:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                if not qual:
                    raise ValueError("Truncated FASTQ record detected.")
                header = header.strip()
                if header.startswith("@"):
                    header = header[1:]
                header = header.split()[0]
                fout.write(f">{header}\n")
                fout.write(seq.strip() + "\n")


def run_cutadapt_paired(
    r1: Path,
    r2: Path,
    out_r1: Path,
    out_r2: Path,
    adapters: dict[str, str],
    threads: int,
    error_rate: float,
    overlap: int,
    min_length: int,
    dry_run: bool,
) -> None:
    cmd = resolve_cutadapt_cmd()
    cmd += [
        "-j",
        str(threads),
        "-g",
        ensure_prefix(adapters["r1_5p"], "^"),
        "-G",
        ensure_prefix(adapters["r2_5p"], "^"),
        "-a",
        ensure_suffix(adapters["r1_3p"], "$"),
        "-A",
        ensure_suffix(adapters["r2_3p"], "$"),
        "-e",
        str(error_rate),
        "-O",
        str(overlap),
        "--minimum-length",
        str(min_length),
        "-o",
        str(out_r1),
        "-p",
        str(out_r2),
        str(r1),
        str(r2),
    ]
    run_command(cmd, dry_run=dry_run)


def run_cutadapt_single(
    in_fastq: Path,
    out_fastq: Path,
    adapters: dict[str, str],
    threads: int,
    error_rate: float,
    overlap: int,
    min_length: int,
    dry_run: bool,
) -> None:
    cmd = resolve_cutadapt_cmd()
    cmd += [
        "-j",
        str(threads),
        "-g",
        ensure_prefix(adapters["r1_5p"], "^"),
        "-g",
        ensure_prefix(adapters["r2_5p"], "^"),
        "-a",
        ensure_suffix(adapters["r1_3p"], "$"),
        "-a",
        ensure_suffix(adapters["r2_3p"], "$"),
        "-e",
        str(error_rate),
        "-O",
        str(overlap),
        "--minimum-length",
        str(min_length),
        "-o",
        str(out_fastq),
        str(in_fastq),
    ]
    run_command(cmd, dry_run=dry_run)


def run_assemblepairs(
    r1: Path,
    r2: Path,
    outname: str,
    coord: str,
    rc: str,
    log_path: Path,
    failed: bool,
    swap_r1r2: bool,
    dry_run: bool,
    outdir: Path,
) -> None:
    assemble_path = resolve_assemblepairs_path()
    first = r1 if swap_r1r2 else r2
    second = r2 if swap_r1r2 else r1
    cmd = [
        sys.executable,
        str(assemble_path),
        "align",
        "-1",
        str(first),
        "-2",
        str(second),
        "--coord",
        coord,
        "--rc",
        rc,
        "--outname",
        outname,
        "--log",
        str(log_path),
    ]
    if failed:
        cmd.append("--failed")
    run_command(cmd, dry_run=dry_run, cwd=outdir)


def main() -> int:
    args = parse_args()
    r1 = Path(args.r1)
    r2 = Path(args.r2)
    if not r1.exists():
        print(f"R1 not found: {r1}", file=sys.stderr)
        return 2
    if not r2.exists():
        print(f"R2 not found: {r2}", file=sys.stderr)
        return 2

    base_outdir = Path(args.outdir) if args.outdir else r1.parent
    base_outdir.mkdir(parents=True, exist_ok=True)

    prefix_base = args.prefix or infer_prefix(r1)
    run_id = args.run_id.strip() if args.run_id else datetime.now().strftime("%Y%m%d_%H%M%S")
    validate_run_id(run_id)
    run_name = f"{prefix_base}_{mode_label(args.mode)}_{run_id}"
    run_dir = base_outdir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    assemble_prefix_arg = args.assemble_prefix or f"{run_name}_presto"
    assemble_prefix_path = Path(assemble_prefix_arg)
    if assemble_prefix_path.is_absolute():
        assemble_outname = assemble_prefix_path
    else:
        assemble_outname = run_dir / assemble_prefix_arg
    assemble_outname.parent.mkdir(parents=True, exist_ok=True)

    if args.log:
        log_path = Path(args.log)
        if not log_path.is_absolute():
            log_path = run_dir / log_path
    else:
        log_path = run_dir / f"{run_name}_AP_align.log"

    adapters = {
        "r1_5p": sanitize_adapter(args.adapter_5p_r1),
        "r2_5p": sanitize_adapter(args.adapter_5p_r2),
        "r1_3p": sanitize_adapter(args.adapter_3p_r1),
        "r2_3p": sanitize_adapter(args.adapter_3p_r2),
    }

    trimmed_r1 = run_dir / f"{run_name}_trim_R1.fastq"
    trimmed_r2 = run_dir / f"{run_name}_trim_R2.fastq"
    merge_then_trim_fastq = run_dir / f"{run_name}.fastq"

    if args.mode == "trim-merge":
        run_cutadapt_paired(
            r1=r1,
            r2=r2,
            out_r1=trimmed_r1,
            out_r2=trimmed_r2,
            adapters=adapters,
            threads=args.threads,
            error_rate=args.cutadapt_error_rate,
            overlap=args.cutadapt_overlap,
            min_length=args.cutadapt_min_length,
            dry_run=args.dry_run,
        )
        assemble_r1 = trimmed_r1
        assemble_r2 = trimmed_r2
    else:
        assemble_r1 = r1
        assemble_r2 = r2

    run_assemblepairs(
        r1=assemble_r1,
        r2=assemble_r2,
        outname=str(assemble_outname),
        coord=args.coord,
        rc=args.rc,
        log_path=log_path,
        failed=args.failed,
        swap_r1r2=args.swap_r1r2,
        dry_run=args.dry_run,
        outdir=run_dir,
    )

    assemble_pass = assemble_outname.parent / f"{assemble_outname.name}_assemble-pass.fastq"
    if not args.dry_run and not assemble_pass.exists():
        print(f"AssemblePairs output not found: {assemble_pass}", file=sys.stderr)
        return 3
    if args.mode == "merge-trim":
        run_cutadapt_single(
            in_fastq=assemble_pass,
            out_fastq=merge_then_trim_fastq,
            adapters=adapters,
            threads=args.threads,
            error_rate=args.cutadapt_error_rate,
            overlap=args.cutadapt_overlap,
            min_length=args.cutadapt_min_length,
            dry_run=args.dry_run,
        )
        if not args.dry_run and not merge_then_trim_fastq.exists():
            print(f"merge-then-trim output not found: {merge_then_trim_fastq}", file=sys.stderr)
            return 4
        fasta_source = merge_then_trim_fastq
    else:
        fasta_source = assemble_pass

    fasta_out = fasta_source.with_suffix(".fasta")
    if args.dry_run:
        print(f">> fastq_to_fasta {fasta_source} -> {fasta_out}")
    else:
        fastq_to_fasta(fasta_source, fasta_out)

    print("Done.")
    print(f"output folder: {run_dir}")
    print(f"merged fastq: {fasta_source}")
    print(f"fasta: {fasta_out}")
    print(f"log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
