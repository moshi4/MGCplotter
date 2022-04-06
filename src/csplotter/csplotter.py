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

    gbk = Genbank(ref_file)
    circos_config = CircosConfig(gbk, outdir)

    # karyotype_file = outdir / "karyotype.txt"
    # circos_config.write_karyotype_file(karyotype_file)

    config_file = outdir / "circos.conf"
    circos_config.write_config_file(config_file)

    sp.run(f"circos -conf {config_file}", shell=True)


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
