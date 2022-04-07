from pathlib import Path
from typing import List, Optional

from csplotter.genbank import Genbank


class CircosConfig:
    """Circos Config Class"""

    def __init__(
        self,
        ref_gbk: Genbank,
        outdir: Path,
        window_size: int = 5000,
        step_size: int = 2000,
    ):
        """Constructor"""
        self.ref_gbk = ref_gbk
        self.outdir = outdir
        self.window_size = window_size
        self.step_size = step_size

        self._config_file = outdir / "circos.conf"
        self._ideogram_file = outdir / "ideogram.conf"
        self._ticks_file = outdir / "ticks.conf"
        self._karyotype_file = outdir / "karyotype.txt"
        self._forward_cds_file = outdir / "forward_cds.txt"
        self._reverse_cds_file = outdir / "reverse_cds.txt"
        self._gc_skew_file = outdir / "gc_skew.txt"
        self._gc_content_file = outdir / "gc_content.txt"

    def write_config_file(self, config_outfile: Path) -> Path:
        """Write circos config file"""
        self._write_karyotype_file()
        self._write_ideogram_conf()
        self._write_ticks_conf()
        forward_cds_conf = self._add_feature(feature_types=["CDS"], target_strand=1)
        gc_content_conf = self._add_gc_content()
        gc_skew_conf = self._add_gc_skew()
        config_contents = (
            f"karyotype = {self._karyotype_file}\n"
            + "chromosomes_units = 1000000\n"
            + f"<<include {self._ideogram_file}>>\n"
            + f"<<include {self._ticks_file}>>\n"
            + "<plots>\n"
            + forward_cds_conf
            + gc_content_conf
            + gc_skew_conf
            + "</plots>\n"
            + "<image>\n"
            + "<<include image.conf>>\n"
            + "</image>\n"
            + "<<include colors_fonts_patterns.conf>> \n"
            + "<<include housekeeping.conf>> \n"
        )

        with open(config_outfile, "w") as f:
            f.write(config_contents)

        return self._config_file

    def _write_ideogram_conf(self) -> None:
        contents = (
            ""
            + "<ideogram>\n"
            + "<spacing>\n"
            + "default = 0.005r\n"
            + "</spacing>\n"
            + "# Ideogram position, fill and outline\n"
            + "radius           = 0.90r\n"
            + "thickness        = 20p\n"
            + "fill             = yes\n"
            + "stroke_color     = dgrey\n"
            + "stroke_thickness = 2p\n"
            + "# Minimum definition for ideogram labels.\n"
            + "show_label       = yes\n"
            + "# see etc/fonts.conf for list of font names\n"
            + "label_font       = default \n"
            + "label_radius     = 1r + 75p\n"
            + "label_size       = 30\n"
            + "label_parallel   = yes\n"
            + "</ideogram>\n"
        )
        with open(self._ideogram_file, "w") as f:
            f.write(contents)

    def _write_ticks_conf(self) -> None:
        contents = (
            ""
            + "show_ticks          = yes\n"
            + "show_tick_labels    = yes\n"
            + "<ticks>\n"
            + "radius           = 1r\n"
            + "color            = black\n"
            + "thickness        = 2p\n"
            + "multiplier       = 1e-6\n"
            + "format           = %d\n"
            + "<tick>\n"
            + "spacing        = 5u\n"
            + "size           = 10p\n"
            + "</tick>\n"
            + "<tick>\n"
            + "spacing        = 25u\n"
            + "size           = 15p\n"
            + "show_label     = yes\n"
            + "label_size     = 20p\n"
            + "label_offset   = 10p\n"
            + "format         = %d\n"
            + "</tick>\n"
            + "</ticks>\n"
        )
        with open(self._ticks_file, "w") as f:
            f.write(contents)

    def _add_feature(
        self, feature_types: List[str], target_strand: Optional[int] = None
    ) -> str:
        self._write_feature_file(feature_types, target_strand)
        contents = (
            f"##### {'-'.join(feature_types)} Features #####\n"
            + "<plot>\n"
            + "type = tile\n"
            + f"file = {self._forward_cds_file}\n"
            + "r0 = 0.95r\n"
            + "r1 = 1.00r\n"
            + "orientation = out\n"
            + "layers = 1\n"
            + "margin = 0.01u\n"
            + "thickness = 50\n"
            + "padding = 1\n"
            + "stroke_thickness = .5\n"
            + "stroke_color = black\n"
            + "layers_overflow = collapse\n"
            + "</plot>\n"
        )
        return contents

    def _write_feature_file(
        self, feature_types: List[str], target_strand: Optional[int]
    ) -> None:
        features = self.ref_gbk.extract_all_features(feature_types, target_strand)
        contents = ""
        for f in features:
            start, end, strand = f.location.start, f.location.end, f.strand
            strand = "+" if strand == 1 else "-"
            contents += f"main {start} {end} {strand} color=ff0000\n"
        with open(self._forward_cds_file, "w") as f:
            f.write(contents)

    def _add_gc_skew(self) -> str:
        abs_max_value = self._write_gc_skew_file()
        contents = (
            "##### GC skew #####\n"
            + "<plot>\n"
            + "type = line\n"
            + f"file = {self._gc_skew_file}\n"
            + "extend_bin = yes\n"
            + "thickness = 0\n"
            + "r0 = 0.2r\n"
            + "r1 = 0.35r\n"
            + "orientation = out\n"
            + f"min = -{abs_max_value}\n"
            + f"max = {abs_max_value}\n"
            + "</plot>\n"
        )
        return contents

    def _write_gc_skew_file(self) -> float:
        gc_skew_values = self.ref_gbk.gc_skew(self.window_size, self.step_size)
        contents = ""
        for i, gc_skew in enumerate(gc_skew_values):
            pos = i * self.step_size
            color = "blue" if gc_skew > 0 else "orange"
            contents += f"main {pos} {pos} {gc_skew} fill_color={color}\n"
        with open(self._gc_skew_file, "w") as f:
            f.write(contents)
        return max(abs(v) for v in gc_skew_values)

    def _add_gc_content(self) -> str:
        abs_max_value = self._write_gc_content_file()
        contents = (
            "##### GC content #####\n"
            + "<plot>\n"
            + "type = line\n"
            + f"file = {self._gc_content_file}\n"
            + "extend_bin = yes\n"
            + "thickness = 0\n"
            + "r0 = 0.35r\n"
            + "r1 = 0.5r\n"
            + "orientation = out\n"
            + f"min = -{abs_max_value}\n"
            + f"max = {abs_max_value}\n"
            + "</plot>\n"
        )
        return contents

    def _write_gc_content_file(self) -> float:
        gc_content_values = self.ref_gbk.gc_content(self.window_size, self.step_size)
        gc_content_values = [v - self.ref_gbk.average_gc for v in gc_content_values]
        contents = ""
        for i, gc_content in enumerate(gc_content_values):
            pos = i * self.step_size
            color = "black" if gc_content > 0 else "grey"
            contents += f"main {pos} {pos} {gc_content} fill_color={color}\n"
        with open(self._gc_content_file, "w") as f:
            f.write(contents)
        return max(abs(v) for v in gc_content_values)

    def _write_karyotype_file(self) -> None:
        genome_length = self.ref_gbk.genome_length
        with open(self._karyotype_file, "w") as f:
            f.write(f"chr - main 1 0 {genome_length} grey")
