#!/usr/bin/env python3
import argparse
import os
import subprocess as sp
from pathlib import Path

from csplotter.circos_config import CircosConfig
from csplotter.genbank import Genbank

__version__ = "0.1.0"


def main():
    """CSplotter main function for entrypoint"""
    # Get argument values
    args = get_args()
    ref_file: Path = args.ref_file
    outdir: Path = args.outdir
    thread_num: int = args.thread_num
    evalue: float = args.evalue

    run(ref_file, outdir, thread_num, evalue)


def run(ref_file: Path, outdir: Path, thread_num: int, evalue: float):
    """Run CSplotter workflow"""
    outdir.mkdir(exist_ok=True)

    config_dir = outdir / "circos_config"
    config_dir.mkdir(exist_ok=True)

    gbk = Genbank(ref_file)
    circos_config = CircosConfig(gbk, config_dir, window_size=5000, step_size=2000)

    # karyotype_file = outdir / "karyotype.txt"
    # circos_config.write_karyotype_file(karyotype_file)

    config_file = config_dir / "circos.conf"
    circos_config.write_config_file(config_file)

    sp.run(f"circos -conf {config_file}", shell=True)

    circos_png_file, circos_svg_file = "circos.png", "circos.svg"
    (outdir / circos_png_file).unlink(missing_ok=True)
    os.rename(circos_png_file, outdir / circos_png_file)
    os.rename(circos_svg_file, outdir / circos_svg_file)


def get_args():
    """Get arguments

    Returns:
        argparse.Namespace: Argument values
    """
    description = "Classify prokaryote protein sequences into COG functional category"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-r",
        "--ref_file",
        required=True,
        type=Path,
        help="Reference genbank file (*.gb|*.gbk|*.gbff)",
        metavar="",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="",
    )
    cpu_num = os.cpu_count()
    default_thread_num = 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
    parser.add_argument(
        "-t",
        "--thread_num",
        type=int,
        help=f"MMseqs threads parameter (Default: {default_thread_num})",
        default=default_thread_num,
        metavar="",
    )
    default_evalue = 1e-2
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        help=f"MMseqs e-value parameter (Default: {default_evalue})",
        default=default_evalue,
        metavar="",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"COGclassifier: v{__version__}",
        help="Print version information",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
