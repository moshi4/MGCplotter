#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess as sp
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import matplotlib as mpl
import pandas as pd
from cogclassifier import cogclassifier

from mgcplotter import config
from mgcplotter.circos_config import CircosConfig
from mgcplotter.genbank import Genbank

__version__ = "0.1.0"


def main():
    """MGCplotter main function for entrypoint"""
    # Copy circos 'etc' config directory to current directory
    # Required if Circos is installed in unusual location
    circos_etc_dir = Path(__file__).parent / "etc"
    shutil.copytree(circos_etc_dir, "etc", dirs_exist_ok=True)

    run(**get_args().__dict__)

    # Delete 'etc' config directory after run
    shutil.rmtree("etc", ignore_errors=True)


def run(
    ref_file: Path,
    outdir: Path,
    query_list: List[Path],
    thread_num: int,
    evalue: float,
    force: bool,
    ticks_labelsize: int = 35,
    # Radius
    forward_cds_r: float = 0.07,
    reverse_cds_r: float = 0.07,
    rrna_r: float = 0.07,
    trna_r: float = 0.07,
    conserved_seq_r: float = 0.04,
    gc_content_r: float = 0.15,
    gc_skew_r: float = 0.15,
    # Color
    forward_cds_color: str = "red",
    reverse_cds_color: str = "blue",
    rrna_color: str = "green",
    trna_color: str = "magenta",
    conserved_seq_color: str = "chocolate",
    gc_content_p_color: str = "black",
    gc_content_n_color: str = "grey",
    gc_skew_p_color: str = "olive",
    gc_skew_n_color: str = "purple",
) -> None:
    """Run MGCplotter workflow"""
    # Setup directory
    config_dir = outdir / "circos_config"
    config_rbh_dir = config_dir / "rbh_results"
    rbh_dir = outdir / "rbh_search"
    outdir.mkdir(exist_ok=True)
    config_dir.mkdir(exist_ok=True)
    rbh_dir.mkdir(exist_ok=True)
    config_rbh_dir.mkdir(exist_ok=True)

    # Search conserved sequence by MMseqs RBH method
    rbh_result_files: List[Path] = []
    ref_gbk = Genbank(ref_file)
    ref_fasta_file = outdir / "reference_cds.faa"
    ref_gbk.write_cds_fasta(ref_fasta_file)
    for query_file in query_list:
        query_fasta_file = rbh_dir / query_file.with_suffix(".faa").name
        if query_file.suffix in config.fasta_suffixs:
            shutil.copy(query_file, query_fasta_file)
        elif query_file.suffix in config.gbk_suffixs:
            Genbank(query_file).write_cds_fasta(query_fasta_file)
        else:
            continue
        ref_name = ref_file.with_suffix("").name
        rbh_result_file = rbh_dir / (
            ref_name + "_vs_" + query_file.with_suffix(".tsv").name
        )
        if force or not rbh_result_file.exists():
            run_mmseqs_rbh_search(
                query_fasta_file, ref_fasta_file, rbh_result_file, evalue, thread_num
            )
        rbh_result_files.append(rbh_result_file)

    # Setup Circos config
    circos_config = CircosConfig(
        ref_gbk=ref_gbk,
        config_dir=config_dir,
        img_dir=outdir,
        ticks_labelsize=ticks_labelsize,
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
        conserved_seq_color=to_hex(conserved_seq_color),
        gc_content_p_color=to_hex(gc_content_p_color),
        gc_content_n_color=to_hex(gc_content_n_color),
        gc_skew_p_color=to_hex(gc_skew_p_color),
        gc_skew_n_color=to_hex(gc_skew_n_color),
    )
    for rbh_result_file in rbh_result_files:
        rbh_config_file = config_rbh_dir / rbh_result_file.with_suffix(".txt").name
        circos_config.add_rbh_config(rbh_result_file, rbh_config_file)
    config_file = config_dir / "circos.conf"
    circos_config.write_config_file(config_file)

    # Run COGclassifier
    cog_dir = outdir / "cogclassifier"
    cog_classifier_result_file = cog_dir / "classifier_result.tsv"
    if force or not cog_classifier_result_file.exists():
        cogclassifier.run(ref_fasta_file, cog_dir, thread_num=thread_num, evalue=1e-2)

    # Assign COG color to reference CDS
    location_id2color = get_location_id2color(
        cog_classifier_result_file, config.cog_letter2color
    )
    rewrite_circos_cds_color(circos_config._f_cds_file, location_id2color)
    rewrite_circos_cds_color(circos_config._r_cds_file, location_id2color)

    # Run Circos
    sp.run(f"circos -conf {config_file}", shell=True)


def run_mmseqs_rbh_search(
    query_fasta_file: Path,
    ref_fasta_file: Path,
    rbh_result_file: Path,
    evalue: float = 1e-3,
    thread_num: int = 1,
) -> None:
    """Run MMseqs rbh search"""
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = (
            f"mmseqs easy-rbh {query_fasta_file} {ref_fasta_file} "
            + f"{rbh_result_file} {tmpdir} -e {evalue} --threads {thread_num}"
        )
        sp.run(cmd, shell=True)


def to_hex(color_like_str: str) -> str:
    """Convert color like string to hexcolor code

    Args:
        color_like_str (str): Color like string

    Returns:
        str: hexcolor code
    """
    return mpl.colors.to_hex(color_like_str).lstrip("#")


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
            color = location_id2color.get(location_id, config.cog_letter2color["-"])
            contents += " ".join(line.split(" ")[0:4]) + f" color={to_hex(color)}\n"
    with open(circos_cds_file, "w") as f:
        f.write(contents)


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns:
        argparse.Namespace: Argument values
    """
    desc = "Microbial Genome Circular plotting tool"
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
    parser.add_argument(
        "--query_list",
        nargs="+",
        type=Path,
        help="Query fasta or genbank files (*.fa|*.faa|*.fasta|*.gb|*.gbk|*.gbff)",
        default=[],
        metavar="",
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
    parser.add_argument(
        "-f",
        "--force",
        help="Forcibly overwrite previous result (Default: OFF)",
        action="store_true",
    )
    default_ticks_labelsize = 35
    parser.add_argument(
        "--ticks_labelsize",
        type=float,
        help=f"Ticks label size (Default: {default_ticks_labelsize})",
        default=default_ticks_labelsize,
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
    args = parser.parse_args()

    # Argument value validation
    err_info = ""
    for k, v in args.__dict__.items():
        if k in config.color_args_dict.keys():
            if not mpl.colors.is_color_like(v):
                err_info += f"'--{k} {v}' is invalid. Not color like string.\n"
        elif k in config.radius_args_dict.keys():
            if not 0 <= v <= 1:
                err_info += f"'--{k} {v}' is invalid. Value must be 0 <= value <= 1.\n"
    if err_info != "":
        parser.error("\n" + err_info)

    return args


if __name__ == "__main__":
    main()
