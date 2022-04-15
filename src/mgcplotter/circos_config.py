from pathlib import Path
from typing import List, Optional

import matplotlib as mpl
import pandas as pd

from mgcplotter.genbank import Genbank


class CircosConfig:
    """Circos Config Class"""

    def __init__(
        self,
        ref_gbk: Genbank,
        config_dir: Path,
        img_dir: Path,
        ticks_labelsize=35,
        # Radius
        forward_cds_r=0.07,
        reverse_cds_r=0.07,
        rrna_r=0.07,
        trna_r=0.07,
        conserved_cds_r=0.04,
        gc_content_r=0.15,
        gc_skew_r=0.15,
        # Color
        forward_cds_color: str = "red",
        reverse_cds_color: str = "blue",
        rrna_color: str = "green",
        trna_color: str = "magenta",
        conserved_cds_color: str = "chocolate",
        gc_content_p_color: str = "black",
        gc_content_n_color: str = "grey",
        gc_skew_p_color: str = "olive",
        gc_skew_n_color: str = "purple",
    ):
        """Constructor"""
        self.ref_gbk = ref_gbk
        self.config_dir = config_dir
        self.img_dir = img_dir
        self.ticks_labelsize = ticks_labelsize
        # Radius
        self.f_cds_r = forward_cds_r
        self.r_cds_r = reverse_cds_r
        self.rrna_r = rrna_r
        self.trna_r = trna_r
        self.conserved_cds_r = conserved_cds_r
        self.gc_content_r = gc_content_r
        self.gc_skew_r = gc_skew_r
        self.separate_r = 0.005
        # Color
        self.f_cds_color = self._to_hex(forward_cds_color)
        self.r_cds_color = self._to_hex(reverse_cds_color)
        self.rrna_color = self._to_hex(rrna_color)
        self.trna_color = self._to_hex(trna_color)
        self.conserved_cds_color = self._to_hex(conserved_cds_color)
        self.gc_content_p_color = self._to_hex(gc_content_p_color)
        self.gc_content_n_color = self._to_hex(gc_content_n_color)
        self.gc_skew_p_color = self._to_hex(gc_skew_p_color)
        self.gc_skew_n_color = self._to_hex(gc_skew_n_color)
        self.separate_color = self._to_hex("grey")

        # Circos config file
        self._config_file = config_dir / "circos.conf"
        self._ideogram_file = config_dir / "ideogram.conf"
        self._ticks_file = config_dir / "ticks.conf"
        self._karyotype_file = config_dir / "karyotype.txt"
        # Features config files
        self._ref_features_dir = config_dir / "reference_features"
        self._ref_features_dir.mkdir(exist_ok=True)
        self._f_cds_file = self._ref_features_dir / "forward_cds.txt"
        self._r_cds_file = self._ref_features_dir / "reverse_cds.txt"
        self._rrna_file = self._ref_features_dir / "rRNA.txt"
        self._trna_file = self._ref_features_dir / "tRNA.txt"
        # GC content & GC skew config file
        self._gc_content_file = config_dir / "gc_content.txt"
        self._gc_skew_file = config_dir / "gc_skew.txt"
        # Separate config files
        self._separate_file = config_dir / "separate.txt"
        # Conserved CDS config files
        self._conserved_cds_dir = config_dir / "conserved_cds"
        self._conserved_cds_dir.mkdir(exist_ok=True)
        self._conserved_cds_files: List[Path] = []

        self._r = 1.0
        self._track_config = ""

    def write_config_file(self, config_outfile: Path) -> None:
        """Write Circos config file

        Args:
            config_outfile (Path): Circos config file
        """
        self._write_ideogram_conf()
        self._write_ticks_conf()
        self._write_karyotype_file()
        self._add_feature_track(
            self._f_cds_file, "CDS", 1, self.f_cds_color, self.f_cds_r
        )
        self._add_feature_track(
            self._r_cds_file, "CDS", -1, self.r_cds_color, self.r_cds_r
        )
        self._add_feature_track(
            self._rrna_file, "rRNA", None, self.rrna_color, self.rrna_r
        )
        self._add_feature_track(
            self._trna_file, "tRNA", None, self.trna_color, self.trna_r
        )
        if len(self._conserved_cds_files) != 0 and self.conserved_cds_r != 0:
            self._r -= 0.01
            self._add_separate_track()
            for conserved_cds_file in self._conserved_cds_files:
                self._add_conservd_cds_track(conserved_cds_file)
            self._add_separate_track()
            self._r -= 0.01
        self._add_gc_content_track()
        self._add_gc_skew_track()

        # Circos overall config
        config_contents = self._concat_lines(
            [
                "karyotype = {0}".format(self._karyotype_file),
                "chromosomes_units = {0}".format(self._chromosome_units),
                "<<include {0}>>".format(self._ideogram_file),
                "<<include {0}>>".format(self._ticks_file),
                "<plots>",
                self._track_config.rstrip("\n"),
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
        contents = f"chr - main 1 0 {self._genome_length} grey\n"
        colors = ["lgrey", "dgrey"]
        base_len = 0
        for idx, contig_seq in enumerate(self.ref_gbk.contig_seqs):
            start, end, color = base_len, base_len + len(contig_seq), colors[idx % 2]
            contents += f"band main band{idx+1} band{idx+1} {start} {end} {color}\n"
            base_len += len(contig_seq)
        with open(self._karyotype_file, "w") as f:
            f.write(contents)

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
                "radius           = 0.80r",
                "thickness        = 15p",
                "fill             = yes",
                "stroke_color     = dgrey",
                "show_bands       = yes",
                "fill_bands       = yes",
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
                "multiplier  = {0}".format(self._ticks_multiplier),
                "orientation = out",
                "format      = {0} {1}".format(self._ticks_format, self._ticks_unit),
                "<tick>",
                "spacing      = {0:.2f}u".format(self._largeticks_spacing),
                "show_label   = yes",
                "label_size   = {0}p".format(self.ticks_labelsize),
                "label_offset = 10p",
                "size         = 25p",
                "</tick>",
                "<tick>",
                "spacing      = {0:.2f}u".format(self._smallticks_spacing),
                "show_label   = no",
                "size         = 15p",
                "</tick>",
                "</ticks>",
            ]
        )
        with open(self._ticks_file, "w") as f:
            f.write(contents)

    ###########################################################################
    # Add Feature(CDS, rRNA, tRNA) track
    ###########################################################################
    def _add_feature_track(
        self,
        feature_file: Path,
        feature_type: str,
        target_strand: Optional[int] = None,
        color: str = "grey",
        feature_r: float = 0.07,
    ) -> None:
        """Add Feature track

        Args:
            feature_file (Path): Feature file to write
            feature_type (str): Feature type (e.g. 'CDS', 'rRNA', 'tRNA')
            target_strand (Optional[int]): Strand ('1', '-1', 'None')
            color (str): Feature color to be drawn
            feature_r (float): Feature radius size
        """
        self._write_feature_file(feature_file, feature_type, target_strand, color)
        self._track_config += self._concat_lines(
            [
                f"##### {feature_type} Feature Track #####",
                "<plot>",
                "type             = tile",
                "file             = {0}".format(feature_file),
                "r1               = {0:.3f}r".format(self._r),
                "r0               = {0:.3f}r".format(self._r - feature_r),
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
        self._r -= feature_r

    def _write_feature_file(
        self,
        feature_file: Path,
        feature_type: str,
        target_strand: Optional[int],
        color: str,
    ) -> None:
        """Write feature file

        Args:
            feature_file (Path): Feature file to write
            feature_type (str): Feature type (e.g. 'CDS', 'rRNA', 'tRNA')
            target_strand (Optional[int]): Strand ('1', '-1', 'None')
            color (str): Feature color to be drawn
        """
        features = self.ref_gbk.extract_all_features(feature_type, target_strand)
        contents = ""
        for f in features:
            start, end, strand = f.location.start, f.location.end, f.strand
            strand = "+" if strand == 1 else "-"
            contents += f"main {start} {end} {strand} color={color}\n"
        with open(feature_file, "w") as f:
            f.write(contents)

    ###########################################################################
    # Add Separate track
    ###########################################################################
    def _add_separate_track(self) -> None:
        """Add separate track"""
        self._write_separate_file()
        self._track_config += self._concat_lines(
            [
                "##### Separate Track #####",
                "<plot>",
                "type             = tile",
                "file             = {0}".format(self._separate_file),
                "r1               = {0:.3f}r".format(self._r),
                "r0               = {0:.3f}r".format(self._r - self.separate_r),
                "orientation      = center",
                "layers           = 1",
                "margin           = 0.01u",
                "thickness        = {0}".format(self.separate_r * 1000),
                "padding          = 1",
                "layers_overflow  = collapse",
                "</plot>",
            ]
        )
        self._r -= self.separate_r

    def _write_separate_file(self) -> None:
        with open(self._separate_file, "w") as f:
            f.write(f"main 0 {self._genome_length} + color={self.separate_color}")

    ###########################################################################
    # Add Conserved CDS track
    ###########################################################################
    def _add_conservd_cds_track(self, conserved_cds_config_file: Path) -> None:
        """Add Conserved CDS track"""
        self._track_config += self._concat_lines(
            [
                "##### Conserved CDS Track #####",
                "<plot>",
                "type             = tile",
                "file             = {0}".format(conserved_cds_config_file),
                "r1               = {0:.3f}r".format(self._r),
                "r0               = {0:.3f}r".format(self._r - self.conserved_cds_r),
                "orientation      = center",
                "layers           = 1",
                "margin           = 0.01u",
                "thickness        = {0}".format(self.conserved_cds_r * 1000),
                "padding          = 1",
                "stroke_color     = black",
                "stroke_thickness = 0",
                "layers_overflow  = collapse",
                "</plot>",
            ]
        )
        self._r -= self.conserved_cds_r

    ###########################################################################
    # Add GC content track
    ###########################################################################
    def _add_gc_content_track(self) -> None:
        """Add GC Content track"""
        self._r = self._boundary if self._r > self._boundary else self._r
        abs_max_value = self._write_gc_content_file()
        self._track_config += self._concat_lines(
            [
                "##### GC Content Track #####",
                "<plot>",
                "type        = histogram",
                "file        = {0}".format(self._gc_content_file),
                "r1          = {0:.3f}r".format(self._r),
                "r0          = {0:.3f}r".format(self._r - self.gc_content_r),
                "min         = {0:.3f}".format(-abs_max_value),
                "max         = {0:.3f}".format(abs_max_value),
                "thickness   = 0",
                "orientation = out",
                "</plot>",
            ]
        )
        self._r -= self.gc_content_r

    def _write_gc_content_file(self) -> float:
        """Write GC Content file"""
        gc_content_values = self.ref_gbk.gc_content(self._window_size, self._step_size)
        gc_content_values = [v - self.ref_gbk.average_gc for v in gc_content_values]
        contents = ""
        for i, value in enumerate(gc_content_values):
            start = i * self._step_size
            end = start + self._step_size
            end = self._genome_length if end > self._genome_length else end
            color = self.gc_content_p_color if value > 0 else self.gc_content_n_color
            contents += f"main {start} {end} {value} fill_color={color}\n"
        with open(self._gc_content_file, "w") as f:
            f.write(contents)
        return max(abs(v) for v in gc_content_values)

    ###########################################################################
    # Add GC skew track
    ###########################################################################
    def _add_gc_skew_track(self) -> None:
        """Add GC Skew track"""
        self._r = self._boundary if self._r > self._boundary else self._r
        abs_max_value = self._write_gc_skew_file()
        self._track_config += self._concat_lines(
            [
                "##### GC Skew Track #####",
                "<plot>",
                "type        = histogram",
                "file        = {0}".format(self._gc_skew_file),
                "r1          = {0:.3f}r".format(self._r),
                "r0          = {0:.3f}r".format(self._r - self.gc_skew_r),
                "min         = {0:.3f}".format(-abs_max_value),
                "max         = {0:.3f}".format(abs_max_value),
                "thickness   = 0",
                "orientation = out",
                "</plot>",
            ]
        )
        self._r -= self.gc_skew_r

    def _write_gc_skew_file(self) -> float:
        """Write GC Skew file"""
        gc_skew_values = self.ref_gbk.gc_skew(self._window_size, self._step_size)
        contents = ""
        for i, value in enumerate(gc_skew_values):
            start = i * self._step_size
            end = start + self._step_size
            end = self._genome_length if end > self._genome_length else end
            color = self.gc_skew_p_color if value > 0 else self.gc_skew_n_color
            contents += f"main {start} {end} {value} fill_color={color}\n"
        with open(self._gc_skew_file, "w") as f:
            f.write(contents)
        return max(abs(v) for v in gc_skew_values)

    ###########################################################################
    # Properties
    ###########################################################################
    @property
    def _genome_length(self) -> int:
        """Genome length"""
        return self.ref_gbk.genome_length

    @property
    def _window_size(self) -> int:
        """Window size for GC content & GC skew calculation"""
        return int(self._genome_length / 1000)

    @property
    def _step_size(self) -> int:
        """Step size for GC content & GC skew calculation"""
        return int(self._window_size * 0.4)

    @property
    def _chromosome_units(self) -> int:
        """Chromosome units"""
        return 10 ** (len(str(self._genome_length)) - 1)

    @property
    def _ticks_format(self) -> str:
        """Ticks format"""
        if self._chromosome_units >= 10**6:
            return "%.1f"
        else:
            return "%d"

    @property
    def _ticks_multiplier(self) -> float:
        """Ticks multiplier"""
        if self._chromosome_units >= 10**6:
            return 1e-6
        else:
            return 1e-3

    @property
    def _ticks_unit(self) -> str:
        """Ticks unit"""
        if self._chromosome_units >= 10**6:
            return "Mb"
        else:
            return "Kb"

    @property
    def _largeticks_spacing(self) -> float:
        """Largeticks spacing"""
        if self._genome_length / self._chromosome_units >= 2:
            return 0.5
        else:
            return 0.1

    @property
    def _smallticks_spacing(self) -> float:
        """Smallticks spacing"""
        if self._genome_length / self._chromosome_units >= 2:
            return 0.1
        else:
            return 0.01

    @property
    def _boundary(self) -> float:
        """Boundary radius for GC content & GC skew track"""
        return 0.7

    ###########################################################################
    # Add Conserved CDS config function
    ###########################################################################
    def add_conserved_cds_config(self, rbh_result_file: Path) -> None:
        """Add conserved CDS config from MMseqs RBH result

        Args:
            rbh_result_file (Path): MMseqs RBH result file
        """
        df = self._load_rbh_result(rbh_result_file)
        contents = ""
        for query, ident in zip(df["TARGET"], df["FIDENT"]):
            start, end, strand = str(query).split("|")[1].split("_")
            color = self._get_interpolated_color(self.conserved_cds_color, ident)
            contents += f"main {start} {end} {strand} color={color}\n"

        filename = rbh_result_file.with_suffix(".txt").name
        conserved_cds_config_file = self._conserved_cds_dir / filename
        with open(conserved_cds_config_file, "w") as f:
            f.write(contents)
        self._conserved_cds_files.append(conserved_cds_config_file)

    def _load_rbh_result(self, rbh_result_file: Path) -> pd.DataFrame:
        header_names = (
            "QUERY,TARGET,FIDENT,ALNLEN,MISMATCH,GAPOPEN,"
            + "QSTART,QEND,TSTART,TEND,EVALUE,BITS"
        ).split(",")
        df = pd.read_table(rbh_result_file, header=None, names=header_names)
        return df.drop_duplicates(subset="TARGET").sort_values("TARGET")

    def _get_interpolated_color(
        self, hexcolor: str, interpolate_value: float, vmin: float = 0.0
    ) -> str:
        """Get interpolate color from float value (Interpolate: Target color -> white)

        Args:
            hexcolor (str): Target hexcolor
            interpolate_value (float): Interpolate float value (0.0 - 1.0)
            vmin (float): Minimum under value for interpolation

        Returns:
            str: Interpolated color
        """
        hexcolor = hexcolor if hexcolor.startswith("#") else f"#{hexcolor}"
        cmap = mpl.colors.LinearSegmentedColormap.from_list("cmap", ("white", hexcolor))
        norm = mpl.colors.Normalize(vmin=vmin, vmax=1.0)
        norm_value = norm(interpolate_value)
        return mpl.colors.to_hex(cmap(norm_value)).lstrip("#")

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

    def _to_hex(self, color_like_str: str) -> str:
        """Convert color like string to hexcolor code

        Args:
            color_like_str (str): Color like string

        Returns:
            str: hexcolor code

        Notes:
            Circos cannot accept '#' character in hexcolor code
        """
        return mpl.colors.to_hex(color_like_str).lstrip("#")
