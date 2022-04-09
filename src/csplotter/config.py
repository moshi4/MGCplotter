import colorsys
from dataclasses import dataclass
from typing import Union

from matplotlib import colors


@dataclass
class Arg:
    """Arg DataClass"""

    default: Union[int, float, str]
    desc: str


radius_args_dict = {
    "forward_cds_r": Arg(0.07, "Forward CDS track radius size"),
    "reverse_cds_r": Arg(0.07, "Reverse CDS track radius size"),
    "rrna_r": Arg(0.07, "rRNA track radius size"),
    "trna_r": Arg(0.07, "tRNA track radius size"),
    "conserved_seq_r": Arg(0.05, "Conserved seq track radius size"),
    "gc_content_r": Arg(0.15, "GC content track radius size"),
    "gc_skew_r": Arg(0.15, "GC skew track radius size"),
}

color_args_dict = {
    "forward_cds_color": Arg("red", "Forward CDS color"),
    "reverse_cds_color": Arg("blue", "Reverse CDS color"),
    "rrna_color": Arg("green", "rRNA color"),
    "trna_color": Arg("magenta", "tRNA color"),
    "gc_content_p_color": Arg("black", "GC content color for plus value from average"),
    "gc_content_n_color": Arg("grey", "GC content color for minus value from average"),
    "gc_skew_p_color": Arg("olive", "GC skew color for positive value"),
    "gc_skew_n_color": Arg("purple", "GC skew color for negative value"),
}

cog_letter2color = {
    "J": "#FCCCFC",  # Translation, ribosomal structure and biogenesis
    "A": "#FCDCFC",  # RNA processing and modification
    "K": "#FCDCEC",  # Transcription
    "L": "#FCDCDC",  # Replication, recombination and repair
    "B": "#FCDCCC",  # Chromatin structure and dynamics
    "D": "#FCFCDC",  # Cell cycle control, cell division, chromosome partitioning
    "Y": "#FCFCCC",  # Nuclear structure
    "V": "#FCFCBC",  # Defense mechanisms
    "T": "#FCFCAC",  # Signal transduction mechanisms
    "M": "#ECFCAC",  # Cell wall/membrane/envelope biogenesis
    "N": "#DCFCAC",  # Cell motility
    "Z": "#CCFCAC",  # Cytoskeleton
    "W": "#BCFCAC",  # Extracellular structures
    "U": "#ACFCAC",  # Intracellular trafficking, secretion, and vesicular transport
    "O": "#9CFCAC",  # Posttranslational modification, protein turnover, chaperones
    "X": "#9CFC9C",  # Mobilome: prophages, transposons
    "C": "#BCFCFC",  # Energy production and conversion
    "G": "#CCFCFC",  # Carbohydrate transport and metabolism
    "E": "#DCFCFC",  # Amino acid transport and metabolism
    "F": "#DCECFC",  # Nucleotide transport and metabolism
    "H": "#DCDCFC",  # Coenzyme transport and metabolism
    "I": "#DCCCFC",  # Lipid transport and metabolism
    "P": "#CCCCFC",  # Inorganic ion transport and metabolism
    "Q": "#BCCCFC",  # Secondary metabolites biosynthesis, transport and catabolism
    "R": "#E0E0E0",  # General function prediction only
    "S": "#CCCCCC",  # Function unknown
    "-": "#B8B8B8",  # No COG classified
}


def change_luminance(hexcolor: str, value: float) -> str:
    """Change color luminance

    Args:
        hexcolor (str): Hexcolor code
        value (float): Change value for luminance

    Returns:
        str: Changed hexcolor
    """
    rgb = colors.to_rgb(hexcolor)
    hls = colorsys.rgb_to_hls(*rgb)
    hue, luminance, saturation = hls
    new_luminance = luminance + value
    if new_luminance > 1:
        new_luminance = 1.0
    elif new_luminance < 0:
        new_luminance = 0.0
    new_hls = [hue, new_luminance, saturation]
    new_rgb = colorsys.hls_to_rgb(*new_hls)
    return colors.to_hex(new_rgb)


for letter, color in cog_letter2color.items():
    cog_letter2color[letter] = change_luminance(color, -0.3)
