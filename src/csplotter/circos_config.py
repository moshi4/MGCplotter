from pathlib import Path

from csplotter.genbank import Genbank


class CircosConfig:
    """Circos Config Class"""

    def __init__(
        self,
        ref_gbk: Genbank,
        outdir: Path,
        window_size: int = 5000,
        step_size: int = 200,
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
        self._gc_skew_file = outdir / "gc_skew.txt"

        self._karyotype_filename = "karyotype.txt"

    def write_config_file(self, config_outfile: Path) -> Path:
        """Write circos config file"""
        self._write_karyotype_file()
        gc_skew_contents = self._add_gc_skew()
        config_contents = (
            f"karyotype = {self._karyotype_file}\n"
            + "chromosomes_units = 1000000\n"
            + f"<<include {self._ideogram_file}>>\n"
            + f"<<include {self._ticks_file}>>\n"
            + gc_skew_contents
            + "<image>\n"
            + "<<include image.conf>>\n"
            + "</image>\n"
            + "<<include colors_fonts_patterns.conf>> \n"
            + "<<include housekeeping.conf>> \n"
        )

        with open(config_outfile, "w") as f:
            f.write(config_contents)

        return self._config_file

    def _add_gc_skew(self) -> str:
        self._write_gc_skew_file()
        contents = (
            "<plots>\n"
            + "<plot>\n"
            + "type = histogram\n"
            + f"file = {self._gc_skew_file}\n"
            + "extend_bin = yes\n"
            + "thickness = 0\n"
            + "r0 = 0.6r\n"
            + "r1 = 1.0r\n"
            + "orientation = out\n"
            + "min = -0.2\n"
            + "max = 0.2\n"
            + "</plot>\n"
            + "</plots>\n"
        )
        return contents

    def _write_gc_skew_file(self) -> None:
        gc_skew_values = self.ref_gbk.gc_skew(self.window_size, self.step_size)
        contents = ""
        for i, gc_skew in enumerate(gc_skew_values):
            pos = i * self.step_size
            color = "blue" if gc_skew > 0 else "orange"
            contents += f"main {pos} {pos} {gc_skew} fill_color={color}\n"
        with open(self._gc_skew_file, "w") as f:
            f.write(contents)

    def _write_karyotype_file(self) -> None:
        genome_length = self.ref_gbk.genome_length
        with open(self._karyotype_file, "w") as f:
            f.write(f"chr - main 1 0 {genome_length} grey")
