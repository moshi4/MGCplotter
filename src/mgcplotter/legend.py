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


def plot_legend(legends: List[Legend], legend_outfile: Path) -> None:
    """Plot legend

    Args:
        legends (List[Legend]): Legend list
        legend_outfile (Path): Legend output file
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
    legend = plt.legend(handles, descs, frameon=False)
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
    if cc.conserved_seq_r != 0:
        for idx, f in enumerate(cc._rbh_config_files, 1):
            desc = f"Query{idx:02d}: {f.with_suffix('').name}"
            legends.append(Legend(f"#{cc.conserved_seq_color}", desc, "s"))
    if cc.gc_content_r != 0:
        legends.append(Legend(f"#{cc.gc_content_p_color}", "GC Content (+)", "^"))
        legends.append(Legend(f"#{cc.gc_content_n_color}", "GC Content (-)", "v"))
    if cc.gc_skew_r != 0:
        legends.append(Legend(f"#{cc.gc_skew_p_color}", "GC Skew (+)", "^"))
        legends.append(Legend(f"#{cc.gc_skew_n_color}", "GC Skew (-)", "v"))

    plot_legend(legends, legend_outfile)


def plot_cog_legend(
    cog_letter2color: Dict[str, str],
    cog_letter2desc: Dict[str, str],
    legend_outfile: Path,
) -> None:
    """Plot COG classification legend

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
