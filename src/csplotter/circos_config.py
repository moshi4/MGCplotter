from pathlib import Path
from typing import List, Optional

from csplotter.genbank import Genbank


class CircosConfig:
    """Circos Config Class"""

    def __init__(
        self,
        ref_gbk: Genbank,
        config_dir: Path,
        img_dir: Path,
        # Radius
        forward_cds_r=0.07,
        reverse_cds_r=0.07,
        rrna_r=0.07,
        trna_r=0.07,
        conserved_seq_r=0.05,
        gc_content_r=0.15,
        gc_skew_r=0.15,
        # Color
        forward_cds_color: str = "ff0000",  # red
        reverse_cds_color: str = "0000ff",  # blue
        rrna_color: str = "008000",  # green
        trna_color: str = "800080",  # magenta
        gc_content_p_color: str = "000000",  # black
        gc_content_n_color: str = "808080",  # grey
        gc_skew_p_color: str = "808000",  # olive
        gc_skew_n_color: str = "800080",  # purple
    ):
        """Constructor"""
        self.ref_gbk = ref_gbk
        self.config_dir = config_dir
        self.img_dir = img_dir
        # Radius
        self.forward_cds_r = forward_cds_r
        self.reverse_cds_r = reverse_cds_r
        self.rrna_r = rrna_r
        self.trna_r = trna_r
        self.conserved_seq_r = conserved_seq_r
        self.gc_content_r = gc_content_r
        self.gc_skew_r = gc_skew_r
        # Color
        self.forward_cds_color = forward_cds_color
        self.reverse_cds_color = reverse_cds_color
        self.rrna_color = rrna_color
        self.trna_color = trna_color
        self.gc_content_p_color = gc_content_p_color
        self.gc_content_n_color = gc_content_n_color
        self.gc_skew_p_color = gc_skew_p_color
        self.gc_skew_n_color = gc_skew_n_color

        # Circos config file
        self._config_file = config_dir / "circos.conf"
        self._ideogram_file = config_dir / "ideogram.conf"
        self._ticks_file = config_dir / "ticks.conf"
        self._karyotype_file = config_dir / "karyotype.txt"
        self._forward_cds_file = config_dir / "feature_forward_cds.txt"
        self._reverse_cds_file = config_dir / "feature_reverse_cds.txt"
        self._rrna_file = config_dir / "feature_rRNA.txt"
        self._trna_file = config_dir / "feature_tRNA.txt"
        self._gc_skew_file = config_dir / "gc_skew.txt"
        self._gc_content_file = config_dir / "gc_content.txt"

        self._r_counter = 1.0

    def write_config_file(self, config_outfile: Path) -> None:
        """Write Circos config file

        Args:
            config_outfile (Path): Circos config file
        """
        # Ideogram config
        self._write_ideogram_conf()
        # Ticks config
        self._write_ticks_conf()
        # Karyotype txt
        self._write_karyotype_file()
        # Feature config
        feature_conf = ""
        feature_conf += self._add_feature(
            self._forward_cds_file,
            ["CDS"],
            1,
            self.forward_cds_color,
            self.forward_cds_r,
        )
        feature_conf += self._add_feature(
            self._reverse_cds_file,
            ["CDS"],
            -1,
            self.reverse_cds_color,
            self.reverse_cds_r,
        )
        feature_conf += self._add_feature(
            self._rrna_file, ["rRNA"], None, self.rrna_color, self.rrna_r
        )
        feature_conf += self._add_feature(
            self._trna_file, ["tRNA"], None, self.trna_color, self.trna_r
        )
        # GC content config
        gc_content_conf = self._add_gc_content()
        # GC skew config
        gc_skew_conf = self._add_gc_skew()

        # Circos overall config
        config_contents = self._concat_lines(
            [
                "karyotype = {0}".format(self._karyotype_file),
                "chromosomes_units = 1000000",
                "<<include {0}>>".format(self._ideogram_file),
                "<<include {0}>>".format(self._ticks_file),
                "<plots>",
                feature_conf.rstrip("\n"),
                gc_content_conf.rstrip("\n"),
                gc_skew_conf.rstrip("\n"),
                "</plots>",
                "<image>",
                "<<include image.conf>>",
                "dir* = {0}".format(self.img_dir),
                "</image>",
                "<<include colors_fonts_patterns.conf>> ",
                "<<include housekeeping.conf>> ",
            ]
        )
        with open(config_outfile, "w") as f:
            f.write(config_contents)

    ###########################################################################
    # Karyotype config
    ###########################################################################
    def _write_karyotype_file(self) -> None:
        """Write karyotype txt"""
        genome_length = self.ref_gbk.genome_length
        with open(self._karyotype_file, "w") as f:
            f.write(f"chr - main 1 0 {genome_length} grey")

    ###########################################################################
    # Ideogram config
    ###########################################################################
    def _write_ideogram_conf(self) -> None:
        """Write Circos Ideogram config"""
        contents = self._concat_lines(
            [
                "<ideogram>",
                "<spacing>",
                "default = 0.005r",
                "</spacing>",
                "radius           = 0.85r",
                "thickness        = 10p",
                "fill             = yes",
                "stroke_color     = dgrey",
                "stroke_thickness = 2p",
                "show_label       = no",
                "label_font       = default",
                "label_radius     = 1r + 75p",
                "label_size       = 30",
                "label_parallel   = yes",
                "</ideogram>",
            ]
        )
        with open(self._ideogram_file, "w") as f:
            f.write(contents)

    ###########################################################################
    # Ticks config
    ###########################################################################
    def _write_ticks_conf(self) -> None:
        """Write Circos ticks config"""
        contents = self._concat_lines(
            [
                "show_ticks       = yes",
                "show_tick_labels = yes",
                "<ticks>",
                "radius      = 1r",
                "color       = black",
                "thickness   = 2p",
                "multiplier  = {0}".format(self.ticks_multiplier),
                "orientation = out",
                "format      = {0} {1}".format(self.ticks_format, self.ticks_unit),
                "<tick>",
                "spacing      = 1u",
                "size         = 20p",
                "show_label   = yes",
                "label_size   = 40p",
                "label_offset = 20p",
                "</tick>",
                "<tick>",
                "spacing      = 0.2u",
                "size         = 10p",
                "show_label   = yes",
                "label_size   = 40p",
                "label_offset = 10p",
                "format       = {0}".format(self.ticks_format),
                "</tick>",
                "</ticks>",
            ]
        )
        with open(self._ticks_file, "w") as f:
            f.write(contents)

    ###########################################################################
    # Feature(CDS, rRNA, tRNA) config
    ###########################################################################
    def _add_feature(
        self,
        feature_file: Path,
        feature_types: List[str],
        target_strand: Optional[int] = None,
        color: str = "grey",
        feature_r: float = 0.07,
    ) -> str:
        """Add Feature track

        Args:
            feature_file (Path): Feature file to write
            feature_types (List[str]): Feature types (e.g. 'CDS', 'rRNA', 'tRNA')
            target_strand (Optional[int]): Strand ('1', '-1', 'None')
            color (str): Feature color to be drawn
            feature_r (float): Feature radius size
        """
        self._write_feature_file(feature_file, feature_types, target_strand, color)
        contents = self._concat_lines(
            [
                f"##### {'-'.join(feature_types)} Features #####",
                "<plot>",
                "type             = tile",
                "file             = {0}".format(feature_file),
                "r1               = {0:.3f}r".format(self._r_counter),
                "r0               = {0:.3f}r".format(self._r_counter - feature_r),
                "orientation      = out",
                "layers           = 1",
                "margin           = 0.01u",
                "thickness        = {0}".format(feature_r * 1000),
                "padding          = 1",
                "stroke_color     = black",
                "stroke_thickness = 0",
                "layers_overflow  = collapse",
                "</plot>",
            ]
        )
        self._r_counter -= feature_r
        return contents

    def _write_feature_file(
        self,
        feature_file: Path,
        feature_types: List[str],
        target_strand: Optional[int],
        color: str,
    ) -> None:
        """Write feature file

        Args:
            feature_file (Path): Feature file to write
            feature_types (List[str]): Feature types (e.g. 'CDS', 'rRNA', 'tRNA')
            target_strand (Optional[int]): Strand ('1', '-1', 'None')
            color (str): Feature color to be drawn
        """
        features = self.ref_gbk.extract_all_features(feature_types, target_strand)
        contents = ""
        for f in features:
            start, end, strand = f.location.start, f.location.end, f.strand
            strand = "+" if strand == 1 else "-"
            contents += f"main {start} {end} {strand} color={color}\n"
        with open(feature_file, "w") as f:
            f.write(contents)

    ###########################################################################
    # GC content config
    ###########################################################################
    def _add_gc_content(self) -> str:
        """Add GC Content track"""
        self._r_counter = 0.6 if self._r_counter > 0.6 else self._r_counter
        abs_max_value = self._write_gc_content_file()
        contents = self._concat_lines(
            [
                "##### GC Content #####",
                "<plot>",
                "type        = histogram",
                "file        = {0}".format(self._gc_content_file),
                "r1          = {0:.3f}r".format(self._r_counter),
                "r0          = {0:.3f}r".format(self._r_counter - self.gc_content_r),
                "min         = {0:.3f}".format(-abs_max_value),
                "max         = {0:.3f}".format(abs_max_value),
                "thickness   = 0",
                "orientation = out",
                "</plot>",
            ]
        )
        self._r_counter -= self.gc_content_r
        return contents

    def _write_gc_content_file(self) -> float:
        """Write GC Content file"""
        gc_content_values = self.ref_gbk.gc_content(self.window_size, self.step_size)
        gc_content_values = [v - self.ref_gbk.average_gc for v in gc_content_values]
        contents = ""
        for i, value in enumerate(gc_content_values):
            pos = i * self.step_size
            color = self.gc_content_p_color if value > 0 else self.gc_content_n_color
            contents += f"main {pos} {pos} {value} fill_color={color}\n"
        with open(self._gc_content_file, "w") as f:
            f.write(contents)
        return max(abs(v) for v in gc_content_values)

    ###########################################################################
    # GC skew config
    ###########################################################################
    def _add_gc_skew(self) -> str:
        """Add GC Skew track"""
        self._r_counter = 0.6 if self._r_counter > 0.6 else self._r_counter
        abs_max_value = self._write_gc_skew_file()
        contents = self._concat_lines(
            [
                "##### GC Skew #####",
                "<plot>",
                "type        = histogram",
                "file        = {0}".format(self._gc_skew_file),
                "r1          = {0:.3f}r".format(self._r_counter),
                "r0          = {0:.3f}r".format(self._r_counter - self.gc_skew_r),
                "min         = {0:.3f}".format(-abs_max_value),
                "max         = {0:.3f}".format(abs_max_value),
                "thickness   = 0",
                "orientation = out",
                "</plot>",
            ]
        )
        self._r_counter -= self.gc_skew_r
        return contents

    def _write_gc_skew_file(self) -> float:
        """Write GC Skew file"""
        gc_skew_values = self.ref_gbk.gc_skew(self.window_size, self.step_size)
        contents = ""
        for i, value in enumerate(gc_skew_values):
            pos = i * self.step_size
            color = self.gc_skew_p_color if value > 0 else self.gc_skew_n_color
            contents += f"main {pos} {pos} {value} fill_color={color}\n"
        with open(self._gc_skew_file, "w") as f:
            f.write(contents)
        return max(abs(v) for v in gc_skew_values)

    ###########################################################################
    # Properties
    ###########################################################################
    @property
    def window_size(self) -> int:
        """Window size for GC content & skew calculation"""
        return int(self.ref_gbk.genome_length / 1000)

    @property
    def step_size(self) -> int:
        """Step size for GC content & skew calculation"""
        return int(self.window_size * 0.4)

    @property
    def ticks_format(self) -> str:
        """Ticks format"""
        return "%.1f"

    @property
    def ticks_multiplier(self) -> float:
        """Ticks multiplier"""
        return 1e-6

    @property
    def ticks_unit(self) -> str:
        """Ticks unit ('Mb' or 'Kb')"""
        return "Mb"

    ###########################################################################
    # Util functions
    ###########################################################################
    def _concat_lines(self, lines: List[str]) -> str:
        """Concatenate lines

        Args:
            lines (List[str]): Target lines

        Returns:
            str: Concatenated lines string
        """
        return "\n".join(lines) + "\n"
