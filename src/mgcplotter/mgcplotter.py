#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import subprocess as sp
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

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
    mmseqs_evalue: float,
    cog_evalue: float,
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
    assign_cog_color: bool = False,
    cog_color_json: Optional[Path] = None,
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
    outdir.mkdir(exist_ok=True)
    config_dir = outdir / "circos_config"
    config_dir.mkdir(exist_ok=True)
    config_rbh_dir = config_dir / "rbh_results"
    config_rbh_dir.mkdir(exist_ok=True)
    rbh_dir = outdir / "rbh_search"
    rbh_dir.mkdir(exist_ok=True)

    # Search conserved sequence by MMseqs RBH method
    ref_gbk = Genbank(ref_file)
    ref_faa_file = outdir / "reference_cds.faa"
    ref_gbk.write_cds_fasta(ref_faa_file)
    rbh_result_files: List[Path] = []
    for idx, query_file in enumerate(query_list, 1):
        query_num = len(query_list)
        if idx == 1:
            em_print(f"Search Conserved Sequence ({query_num} Query vs Reference)")
        # Setup query CDS faa file
        query_faa_file = rbh_dir / query_file.with_suffix(".faa").name
        if query_file.suffix in config.fasta_suffixs:
            shutil.copy(query_file, query_faa_file)
        elif query_file.suffix in config.gbk_suffixs and not query_faa_file.exists():
            Genbank(query_file).write_cds_fasta(query_faa_file)
        # Run MMseqs RBH search
        query_name = query_file.with_suffix("").name
        ref_name = ref_file.with_suffix("").name
        target_info = f"{query_name} vs {ref_name}[reference]"
        rbh_result_file = rbh_dir / f"{query_name}_vs_reference_rbh.tsv"
        if force or not rbh_result_file.exists():
            print(f"# Run MMseqs RBH search ({target_info})")
            run_mmseqs_rbh_search(
                query_faa_file, ref_faa_file, rbh_result_file, mmseqs_evalue, thread_num
            )
        else:
            print(f"# Reuse previous MMseqs RBH search result ({target_info})")
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
        forward_cds_color=forward_cds_color,
        reverse_cds_color=reverse_cds_color,
        rrna_color=rrna_color,
        trna_color=trna_color,
        conserved_seq_color=conserved_seq_color,
        gc_content_p_color=gc_content_p_color,
        gc_content_n_color=gc_content_n_color,
        gc_skew_p_color=gc_skew_p_color,
        gc_skew_n_color=gc_skew_n_color,
    )
    for rbh_result_file in rbh_result_files:
        rbh_config_file = config_rbh_dir / rbh_result_file.with_suffix(".txt").name
        circos_config.add_rbh_config(rbh_result_file, rbh_config_file)
    config_file = config_dir / "circos.conf"
    circos_config.write_config_file(config_file)

    # Run COGclassifier for Functional Classification of Reference CDSs
    if assign_cog_color:
        cog_dir = outdir / "cogclassifier"
        cog_classifier_result_file = cog_dir / "classifier_result.tsv"
        em_print("Run COGclassifier for Functional Classification of Reference CDSs")
        if force or not cog_classifier_result_file.exists():
            cogclassifier.run(
                ref_faa_file, cog_dir, thread_num=thread_num, evalue=cog_evalue
            )
        else:
            print("# Reuse previous COGclassifier result")

        # Assign COG color to reference CDS
        if cog_color_json is not None:
            with open(cog_color_json) as f:
                config.cog_letter2color = json.load(f)
        location_id2color = get_location_id2color(
            cog_classifier_result_file, config.cog_letter2color
        )
        rewrite_circos_cds_color(circos_config._f_cds_file, location_id2color)
        rewrite_circos_cds_color(circos_config._r_cds_file, location_id2color)

    # Run Circos
    em_print("Run Circos")
    cmd = f"circos -conf {config_file}"
    print(f"$ {cmd}\n")
    sp.run(cmd, shell=True)


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
            f"mmseqs easy-rbh {query_fasta_file} {ref_fasta_file} {rbh_result_file} "
            + f"{tmpdir} -e {evalue} --threads {thread_num} > /dev/null"
        )
        print(f"$ {cmd}\n")
        sp.run(cmd, shell=True)


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
            hexcolor = mpl.colors.to_hex(color).lstrip("#")
            contents += " ".join(line.split(" ")[0:4]) + f" color={hexcolor}\n"
    with open(circos_cds_file, "w") as f:
        f.write(contents)


def em_print(content: str) -> None:
    """Emphasis print content

    Args:
        content (str): Print content
    """
    print(f"\n{'*' * 80}\n* {content}\n{'*'* 80}\n")


def generate_cog_color_template() -> None:
    """Generate COG color template json file"""
    color_template_json_file = "cog_color_template.json"
    with open(color_template_json_file, "w") as f:
        json.dump(config.cog_letter2color, f, indent=2)
    print(f"Generate COG color template json file ('{color_template_json_file}')")


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
        help=f"Query fasta or genbank files ({'|'.join(config.valid_query_suffixs)})",
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
    default_mmseqs_evalue = 1e-5
    parser.add_argument(
        "--mmseqs_evalue",
        type=float,
        help=f"MMseqs e-value parameter (Default: {default_mmseqs_evalue})",
        default=default_mmseqs_evalue,
        metavar="",
    )
    default_cog_evalue = 1e-2
    parser.add_argument(
        "--cog_evalue",
        type=float,
        help=f"COGclassifier e-value parameter (Default: {default_cog_evalue})",
        default=default_cog_evalue,
        metavar="",
    )
    parser.add_argument(
        "-f",
        "--force",
        help="Forcibly overwrite previous calculation result (Default: OFF)",
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
    parser.add_argument(
        "--assign_cog_color",
        help="Assign COG classification color to CDS (Default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "--cog_color_json",
        type=Path,
        help="User-defined COG classification color json file",
        default=None,
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

    # Argument value validation check
    err_info = ""
    for f in args.query_list:
        if f.suffix not in config.valid_query_suffixs:
            err_info += f"'{f.suffix}' is invalid file suffix ({f}).\n"
    for k, v in args.__dict__.items():
        if k in config.color_args_dict.keys():
            if not mpl.colors.is_color_like(v):
                err_info += f"'--{k} {v}' is invalid color like string.\n"
        elif k in config.radius_args_dict.keys():
            if not 0 <= v <= 0.3:
                err_info += f"'--{k} {v}' is invalid value range (0 <= value <= 0.3).\n"

    if args.cog_color_json is not None:
        if not args.cog_color_json.exists():
            err_info += f"--cog_color_json: File not found '{args.cog_color_json}'.\n"
        else:
            with open(args.cog_color_json) as f:
                cog_color_json_dict: Dict[str, str] = json.load(f)
            for k, v in cog_color_json_dict.items():
                if k not in config.cog_letter2color:
                    err_info += f"--cog_color_json: '{k}' is not COG letter.\n"
                if not mpl.colors.is_color_like(v):
                    err_info += f"--cog_color_json: '{v}' is not color like string.\n"

    if err_info != "":
        parser.error("\n" + err_info)

    return args


if __name__ == "__main__":
    main()
