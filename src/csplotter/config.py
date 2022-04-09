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

for letter, color in cog_letter2color.items():
    original_rgb = colors.to_rgb(color)
    original_hsv = colors.rgb_to_hsv(original_rgb)
    converted_hsv = [original_hsv[0], original_hsv[1] + 0.5, original_hsv[2]]
    converted_rgb = colors.hsv_to_rgb(converted_hsv)
    converted_hex = colors.to_hex(converted_rgb)
    cog_letter2color[letter] = converted_hex
