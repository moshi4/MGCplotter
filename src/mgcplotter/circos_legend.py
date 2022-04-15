from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import matplotlib as mpl
import matplotlib.pyplot as plt

from mgcplotter.circos_config import CircosConfig

mpl.rcParams["font.family"] = "monospace"
mpl.rcParams["svg.fonttype"] = "none"


@dataclass
class Legend:
    """Legend DataClass"""

    color: str
    desc: str
    marker: str


class CircosLegend:
    """Circos Legend Class"""

    def __init__(
        self,
        circos_config: CircosConfig,
        cog_letter2color: Dict[str, str],
        cog_letter2desc: Dict[str, str],
        outdir: Path,
        dpi: int = 500,
    ):
        """Constructor"""
        self.circos_config = circos_config
        self.cog_letter2color = cog_letter2color
        self.cog_letter2desc = cog_letter2desc
        self.outdir = outdir
        self.dpi = dpi
        self.outdir.mkdir(exist_ok=True)

        # Legend files
        self._track_contents_png_file = self.outdir / "track_contents.png"
        self._conserved_cds_ident_png_file = self.outdir / "conserved_cds_identity.png"
        self._cog_letter_png_file = self.outdir / "cog_letter.png"
        self._cog_def_png_file = self.outdir / "cog_definition.png"

    def plot_all_legends(self, svg: bool = True) -> None:
        """Plot circos all legends

        Args:
            svg (bool, optional): Output SVG or not
        """
        self.plot_track_contents(self._track_contents_png_file)
        self.plot_cog_letter(self._cog_letter_png_file)
        self.plot_cog_def(self._cog_def_png_file)
        self.plot_conserved_cds_ident(self._conserved_cds_ident_png_file)
        if svg:
            self.plot_track_contents(self._as_svg(self._track_contents_png_file))
            self.plot_cog_letter(self._as_svg(self._cog_letter_png_file))
            self.plot_cog_def(self._as_svg(self._cog_def_png_file))
            self.plot_conserved_cds_ident(
                self._as_svg(self._conserved_cds_ident_png_file)
            )

    def plot_track_contents(self, outfile: Path) -> None:
        """Plot track legend

        Args:
            outfile (Path): Legend output file
        """
        cc = self.circos_config
        legends = []
        if cc.f_cds_r != 0:
            legends.append(Legend(f"#{cc.f_cds_color}", "Forward CDS", "s"))
        if cc.r_cds_r != 0:
            legends.append(Legend(f"#{cc.r_cds_color}", "Reverse CDS", "s"))
        if cc.rrna_r != 0:
            legends.append(Legend(f"#{cc.rrna_color}", "rRNA", "s"))
        if cc.trna_r != 0:
            legends.append(Legend(f"#{cc.trna_color}", "tRNA", "s"))
        if cc.conserved_cds_r != 0 and len(cc._conserved_cds_files) != 0:
            legends.append(Legend(f"#{cc.conserved_cds_color}", "Conserved CDS", "s"))
        if cc.gc_content_r != 0:
            legends.append(Legend(f"#{cc.gc_content_p_color}", "GC Content (+)", "^"))
            legends.append(Legend(f"#{cc.gc_content_n_color}", "GC Content (-)", "v"))
        if cc.gc_skew_r != 0:
            legends.append(Legend(f"#{cc.gc_skew_p_color}", "GC Skew (+)", "^"))
            legends.append(Legend(f"#{cc.gc_skew_n_color}", "GC Skew (-)", "v"))

        self._plot_legend(legends, outfile, title="Track Contents")

    def plot_cog_letter(self, outfile: Path) -> None:
        """Plot COG functional classification letter legend

        Args:
            outfile (Path): Legend output file
        """
        legends = []
        for cog_letter, color in self.cog_letter2color.items():
            legends.append(Legend(color, cog_letter, "s"))
        self._plot_legend(legends, outfile, ncol=6)

    def plot_cog_def(self, outfile) -> None:
        """Plot COG functional classification definition legend

        Args:
            outfile (Path): Legend output file
        """
        legends = []
        for cog_letter, color in self.cog_letter2color.items():
            desc = f"{cog_letter} : {self.cog_letter2desc[cog_letter]}"
            legends.append(Legend(color, desc, "s"))
        self._plot_legend(legends, outfile)

    def plot_conserved_cds_ident(self, outfile: Path) -> None:
        """Plot conserved cds identity legend

        Args:
            outfile (Path): Legend output file

        Notes:
            Reference URL for implementation:
            https://matplotlib.org/2.0.2/examples/api/colorbar_only.html
        """
        fig = plt.figure(figsize=(5, 1))
        ax = fig.add_axes([0.05, 0.5, 0.9, 0.2])  # left, bottom, width, height
        color = "#" + self.circos_config.conserved_cds_color
        cmap = mpl.colors.LinearSegmentedColormap.from_list("cmap", ("white", color))
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        cb = mpl.colorbar.ColorbarBase(  # type: ignore
            ax, cmap=cmap, norm=norm, orientation="horizontal"
        )
        cb.set_label(label="Identity(%)", loc="center")
        cb.ax.invert_xaxis()

        fig.savefig(outfile, dpi=self.dpi)
        fig.clear()

    def _plot_legend(
        self,
        legends: List[Legend],
        legend_outfile: Path,
        title: str = "",
        ncol: int = 1,
    ) -> None:
        """Plot legend

        Args:
            legends (List[Legend]): Legend list
            legend_outfile (Path): Legend output file
            title (str): Legend title
            ncol (str): Number of Column

        Notes:
            Reference URL for implementation:
            https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib
        """
        # Setup matplotlib params for only plot legend (Disable 'spine', 'axis')
        for pos in ["left", "right", "top", "bottom"]:
            plt.gca().spines[pos].set_visible(False)
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)

        handles = []
        for legend in legends:
            marker, color = legend.marker, legend.color
            handles.append(
                plt.plot([], [], marker=marker, color=color, linestyle="none")[0]
            )
        descs = [legend.desc for legend in legends]
        mpl_legend = plt.legend(
            handles,
            descs,
            frameon=False,
            title=title,
            handletextpad=0,
            ncol=ncol,
            columnspacing=1,
        )
        fig = mpl_legend.figure
        fig.canvas.draw()
        bbox = mpl_legend.get_window_extent().transformed(
            fig.dpi_scale_trans.inverted()
        )
        fig.savefig(legend_outfile, dpi=self.dpi, bbox_inches=bbox)
        fig.clear()

    def _as_svg(self, file: Path) -> Path:
        """Convert filename extension to '.svg'"""
        return file.with_suffix(".svg")
