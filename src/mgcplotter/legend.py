from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import matplotlib as mpl
import matplotlib.pyplot as plt

from mgcplotter.circos_config import CircosConfig


@dataclass
class Legend:
    """Legend DataClass"""

    color: str
    desc: str
    marker: str


def plot_legend(
    legends: List[Legend], legend_outfile: Path, title: str = "", ncol: int = 1
) -> None:
    """Plot legend

    Args:
        legends (List[Legend]): Legend list
        legend_outfile (Path): Legend output file
        title (str): Legend title
        ncol (str): Number of Column
    """
    # Setup matplotlib params for only plot legend
    mpl.rcParams["font.family"] = "monospace"
    for pos in ["left", "right", "top", "bottom"]:
        plt.gca().spines[pos].set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    mpl.rcParams["svg.fonttype"] = "none"

    handles = []
    for legend in legends:
        marker, color = legend.marker, legend.color
        handles.append(
            plt.plot([], [], marker=marker, color=color, linestyle="none")[0]
        )
    descs = [legend.desc for legend in legends]
    legend = plt.legend(
        handles,
        descs,
        frameon=False,
        title=title,
        handletextpad=0,
        ncol=ncol,
        columnspacing=1,
    )
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(legend_outfile, dpi=500, bbox_inches=bbox)


def plot_track_legend(circos_config: CircosConfig, legend_outfile: Path) -> None:
    """Plot track legend

    Args:
        circos_config (CircosConfig): CircosConfig object
        legend_outfile (Path): Legend output file
    """
    cc = circos_config
    legends = []
    if cc.f_cds_r != 0:
        legends.append(Legend(f"#{cc.f_cds_color}", "Forward CDS", "s"))
    if cc.r_cds_r != 0:
        legends.append(Legend(f"#{cc.r_cds_color}", "Reverse CDS", "s"))
    if cc.rrna_r != 0:
        legends.append(Legend(f"#{cc.rrna_color}", "rRNA", "s"))
    if cc.trna_r != 0:
        legends.append(Legend(f"#{cc.trna_color}", "tRNA", "s"))
    if cc.conserved_cds_r != 0 and len(cc._rbh_config_files) != 0:
        legends.append(Legend(f"#{cc.conserved_cds_color}", "Conserved CDS", "s"))
    if cc.gc_content_r != 0:
        legends.append(Legend(f"#{cc.gc_content_p_color}", "GC Content (+)", "^"))
        legends.append(Legend(f"#{cc.gc_content_n_color}", "GC Content (-)", "v"))
    if cc.gc_skew_r != 0:
        legends.append(Legend(f"#{cc.gc_skew_p_color}", "GC Skew (+)", "^"))
        legends.append(Legend(f"#{cc.gc_skew_n_color}", "GC Skew (-)", "v"))

    title = "Track Contents"
    plot_legend(legends, legend_outfile, title)


def plot_cog_letter_legend(
    cog_letter2color: Dict[str, str],
    legend_outfile: Path,
) -> None:
    """Plot COG functional classification letter legend

    Args:
        cog_letter2color (Dict[str, str]): COG letter & color dict
        legend_outfile (Path): Legend outpu file
    """
    legends = []
    for cog_letter, color in cog_letter2color.items():
        legends.append(Legend(color, cog_letter, "s"))

    plot_legend(legends, legend_outfile, ncol=6)


def plot_cog_def_legend(
    cog_letter2color: Dict[str, str],
    cog_letter2desc: Dict[str, str],
    legend_outfile: Path,
) -> None:
    """Plot COG functional classification definition legend

    Args:
        cog_letter2color (Dict[str, str]): COG letter & color dict
        cog_letter2desc (Dict[str, str]): COG letter & description dict
        legend_outfile (Path): Legend outpu file
    """
    legends = []
    for cog_letter, color in cog_letter2color.items():
        desc = f"{cog_letter} : {cog_letter2desc[cog_letter]}"
        legends.append(Legend(color, desc, "s"))

    plot_legend(legends, legend_outfile)
