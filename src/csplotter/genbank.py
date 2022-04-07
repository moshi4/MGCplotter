from functools import cached_property
from pathlib import Path
from typing import List, Union

from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import Seq, SeqFeature
from Bio.SeqRecord import SeqRecord


class Genbank:
    """Genbank Class"""

    def __init__(
        self,
        gbk_file: Union[str, Path],
        name: str = "",
    ):
        """Constructor

        Args:
            gbk_file (Union[str, StringIO, Path]): Genbank file
            name (str, optional): Name
        """
        self.gbk_file: Path = Path(gbk_file)
        self.name: str = name if name != "" else self.gbk_file.with_suffix("").name
        self._records: List[SeqRecord] = list(SeqIO.parse(gbk_file, "genbank"))
        self._first_record: SeqRecord = self._records[0]

    @property
    def genome_length(self) -> int:
        """Genome sequence length of Genbank first record"""
        return len(self._first_record.seq)

    def gc_skew(self, window_size: int = 5000, step_size: int = 2000) -> List[float]:
        """GC Skew"""
        gc_skew_values = []
        seq = self._first_record.seq
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos

            subseq = seq[start_pos:end_pos]
            g = subseq.count("G") + subseq.count("g")
            c = subseq.count("C") + subseq.count("c")
            try:
                skew = (g - c) / float(g + c)
            except ZeroDivisionError:
                skew = 0.0
            gc_skew_values.append(skew)
        return gc_skew_values

    def gc_content(self, window_size: int = 5000, step_size: int = 2000) -> List[float]:
        """GC Content"""
        gc_content_values = []
        seq = self._first_record.seq
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos

            subseq = seq[start_pos:end_pos]
            gc_content_values.append(SeqUtils.GC(subseq))
        return gc_content_values

    @cached_property
    def average_gc(self) -> float:
        """Average GC content"""
        return SeqUtils.GC(self._first_record.seq)

    def extract_all_features(
        self,
        only_first_record: bool = False,
        feature_types: List[str] = ["CDS"],
    ) -> List[SeqFeature]:
        """Extract all features

        Args:
            only_first_record (bool): Exatract from only first record or not
            feature_types (List[str]): Feature types to extract
        Returns:
            List[SeqFeature]: All features
        """
        if only_first_record:
            return [f for f in self._first_record.features if f.type in feature_types]
        else:
            features = []
            for record in self._records:
                for f in record.features:
                    if f.type in feature_types:
                        features.append(f)
            return features

    def write_cds_protein_features_fasta(
        self,
        fasta_outfile: Union[str, Path],
        only_first_record: bool = False,
    ):
        """Write CDS protein features fasta"""
        features = self.extract_all_features(only_first_record, ["CDS"])
        cds_seq_records: List[SeqRecord] = []
        for idx, feature in enumerate(features, 1):
            qualifiers = feature.qualifiers
            protein_id = qualifiers.get("protein_id", [None])[0]
            product = qualifiers.get("product", [""])[0]
            translation = qualifiers.get("translation", [None])[0]

            if translation is None:
                continue
            else:
                cds_seq = Seq(translation)

            seq_id = f"GENE{idx:06d}"
            seq_id = seq_id if protein_id is None else protein_id

            cds_seq_record = SeqRecord(seq=cds_seq, id=seq_id, description=product)
            cds_seq_records.append(cds_seq_record)

        SeqIO.write(cds_seq_records, fasta_outfile, "fasta-2line")

    def write_genome_fasta(
        self,
        outfile: Union[str, Path],
    ) -> None:
        """Write genome fasta file

        Args:
            outfile (Union[str, Path]): Output genome fasta file
        """
        write_seq = self._first_record.seq
        with open(outfile, "w") as f:
            f.write(f">{self.name}\n{write_seq}\n")
