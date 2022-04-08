from dataclasses import dataclass
from typing import Union


@dataclass
class Arg:
    """Arg DataClass"""

    default: Union[int, float, str]
    desc: str


radius_args_dict = {
    "ref_feature_r": Arg(0.07, "Reference species feature track radius size"),
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
