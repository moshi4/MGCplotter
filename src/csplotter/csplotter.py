#!/usr/bin/env python3
import argparse
import os
import subprocess as sp
from pathlib import Path

from matplotlib import colors

from csplotter import config
from csplotter.circos_config import CircosConfig
from csplotter.genbank import Genbank

__version__ = "0.1.0"


def main():
    """CSplotter main function for entrypoint"""
    # Get argument values
    run(**get_args().__dict__)


def run(
    ref_file: Path,
    outdir: Path,
    thread_num: int,
    evalue: float,
    ref_feature_r: float = 0.05,
    conserved_seq_r: float = 0.05,
    gc_content_r: float = 0.15,
    gc_skew_r: float = 0.15,
    forward_cds_color: str = "red",
    reverse_cds_color: str = "blue",
    rrna_color: str = "green",
    trna_color: str = "magenta",
    gc_content_p_color: str = "black",
    gc_content_n_color: str = "grey",
    gc_skew_p_color: str = "olive",
    gc_skew_n_color: str = "purple",
):
    """Run CSplotter workflow"""
    outdir.mkdir(exist_ok=True)

    config_dir = outdir / "circos_config"
    config_dir.mkdir(exist_ok=True)

    ref_gbk = Genbank(ref_file)
    circos_config = CircosConfig(
        ref_gbk=ref_gbk,
        config_dir=config_dir,
        img_dir=outdir,
        window_size=5000,
        step_size=2000,
        ref_feature_r=ref_feature_r,
        conserved_seq_r=conserved_seq_r,
        gc_content_r=gc_content_r,
        gc_skew_r=gc_skew_r,
        forward_cds_color=to_hex(forward_cds_color),
        reverse_cds_color=to_hex(reverse_cds_color),
        rrna_color=to_hex(rrna_color),
        trna_color=to_hex(trna_color),
        gc_content_p_color=to_hex(gc_content_p_color),
        gc_content_n_color=to_hex(gc_content_n_color),
        gc_skew_p_color=to_hex(gc_skew_p_color),
        gc_skew_n_color=to_hex(gc_skew_n_color),
    )

    config_file = config_dir / "circos.conf"
    circos_config.write_config_file(config_file)

    sp.run(f"circos -conf {config_file}", shell=True)


def to_hex(color_like_str: str) -> str:
    """Convert color like string to hexcolor code

    Args:
        color_like_str (str): Color like string

    Returns:
        str: hexcolor code
    """
    return colors.to_hex(color_like_str).lstrip("#")


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns:
        argparse.Namespace: Argument values
    """
    desc = "Circos plot tool for CS(Conserved Sequence)"
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        "-r",
        "--ref_file",
        required=True,
        type=Path,
        help="Reference genbank file (*.gb|*.gbk|*.gbff)",
        metavar="R",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="O",
    )
    cpu_num = os.cpu_count()
    default_thread_num = 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
    parser.add_argument(
        "-t",
        "--thread_num",
        type=int,
        help=f"Threads number parameter (Default: {default_thread_num})",
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
    # Track radius control arguments
    for k, v in config.radius_args_dict.items():
        parser.add_argument(
            f"--{k}",
            type=float,
            help=f"{v.desc} (Default: {v.default})",
            default=v.default,
            metavar="",
        )
    # Color control arguments
    for k, v in config.color_args_dict.items():
        parser.add_argument(
            f"--{k}",
            type=str,
            help=f"{v.desc} (Default: '{v.default}')",
            default=v.default,
            metavar="",
        )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"v{__version__}",
        help="Print version information",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
