#!/usr/bin/env python3
import argparse
import os
import subprocess as sp
from collections import defaultdict
from pathlib import Path
from typing import Dict

import pandas as pd
from cogclassifier import cogclassifier
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
    # Radius
    forward_cds_r: float = 0.07,
    reverse_cds_r: float = 0.07,
    rrna_r: float = 0.07,
    trna_r: float = 0.07,
    conserved_seq_r: float = 0.05,
    gc_content_r: float = 0.15,
    gc_skew_r: float = 0.15,
    # Color
    forward_cds_color: str = "red",
    reverse_cds_color: str = "blue",
    rrna_color: str = "green",
    trna_color: str = "magenta",
    gc_content_p_color: str = "black",
    gc_content_n_color: str = "grey",
    gc_skew_p_color: str = "olive",
    gc_skew_n_color: str = "purple",
):
    """Run MGCplotter workflow"""
    # Setup directory
    config_dir = outdir / "circos_config"
    outdir.mkdir(exist_ok=True)
    config_dir.mkdir(exist_ok=True)

    # TODO: Search conserved sequence

    # Setup Circos config
    ref_gbk = Genbank(ref_file)
    circos_config = CircosConfig(
        ref_gbk=ref_gbk,
        config_dir=config_dir,
        img_dir=outdir,
        # Radius
        forward_cds_r=forward_cds_r,
        reverse_cds_r=reverse_cds_r,
        rrna_r=rrna_r,
        trna_r=trna_r,
        conserved_seq_r=conserved_seq_r,
        gc_content_r=gc_content_r,
        gc_skew_r=gc_skew_r,
        # Color
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

    # Run COGclassifier
    ref_cds_fasta_file = outdir / "reference_cds.faa"
    cog_dir = outdir / "cogclassifier"
    cog_dl_dir = cog_dir / "cog_download"
    cog_classifier_result_file = cog_dir / "classifier_result.tsv"
    ref_gbk.write_cds_fasta(ref_cds_fasta_file)
    if not cog_classifier_result_file.exists():
        cogclassifier.run(ref_cds_fasta_file, cog_dir, cog_dl_dir, thread_num, 1e-2)

    # Assign COG color to reference CDS
    location_id2color = get_location_id2color(
        cog_classifier_result_file, config.cog_letter2color
    )
    rewrite_circos_cds_color(circos_config._forward_cds_file, location_id2color)
    rewrite_circos_cds_color(circos_config._reverse_cds_file, location_id2color)

    # Run Circos
    sp.run(f"circos -conf {config_file}", shell=True)


def to_hex(color_like_str: str) -> str:
    """Convert color like string to hexcolor code

    Args:
        color_like_str (str): Color like string

    Returns:
        str: hexcolor code
    """
    return colors.to_hex(color_like_str).lstrip("#")


def get_location_id2color(
    cog_classifier_result_file: Path,
    cog_letter2color: Dict[str, str],
) -> Dict[str, str]:
    """Get CDS location ID & Color dict

    Args:
        cog_classifier_result_file (Path): COGclassifier result file
        cog_letter2color (Dict[str, str]): COG letter & Color dict

    Returns:
        Dict[str, str]: CDS location ID & COG Color dict

    Notes:
        CDS location ID = "start end strand" (e.g. "300 1000 +")
    """
    df = pd.read_csv(cog_classifier_result_file, delimiter="\t")
    location_id2color = defaultdict(str)
    for query_id, cog_letter in zip(df["QUERY_ID"], df["COG_LETTER"]):
        location_id = query_id.split("|")[1].replace("_", " ")
        location_id2color[location_id] = cog_letter2color[cog_letter]
    return location_id2color


def rewrite_circos_cds_color(
    circos_cds_file: Path, location_id2color: Dict[str, str]
) -> None:
    """Rewrite Circos CDS color to COG classification color

    Args:
        circos_cds_file (Path): Circos CDS file
        location_id2color (Dict[str, str]): CDS location ID & COG Color dict
    """
    contents = ""
    with open(circos_cds_file) as f:
        for line in f.read().splitlines():
            location_id = " ".join(line.split(" ")[1:4])
            color = location_id2color.get(location_id, None)
            if color is None:
                color = config.cog_letter2color["-"]
                contents += " ".join(line.split(" ")[0:4]) + f" color={to_hex(color)}\n"
            else:
                contents += " ".join(line.split(" ")[0:4]) + f" color={to_hex(color)}\n"
    with open(circos_cds_file, "w") as f:
        f.write(contents)


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
