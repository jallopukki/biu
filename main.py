from pathlib import Path
from typing import Union, Iterator
from abc import ABC, abstractmethod
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


########### An abstract sequences class


class BiologicalSequence(ABC):
    """An abstract class for nucleic acids"""

    def __init__(self, sequence: str):
        self._sequence = sequence.upper()
        if not self.check_alphabet():
            raise ValueError(f"Invalid alphabet for {self.__class__.__name__}")

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, key: Union[int, slice]) -> Union[str, "BiologicalSequence"]:
        if isinstance(key, slice):
            return self.__class__(self._sequence[key])
        return self._sequence[key]

    def __str__(self) -> str:
        return self._sequence

    def check_alphabet(self) -> bool:
        return self._validate_alphabet(self._sequence)

    @abstractmethod
    def _validate_alphabet(self, seq: str) -> bool:
        """DNA/RNA alphabet check (implemented in the corresponding subclasses)"""
        pass


class NucleicAcidSequence(BiologicalSequence):
    """Base class for nucleic acids (DNA/RNA)"""

    def __init__(self, sequence: str):
        if type(self) == NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is an abstract class and cannot be instantiated directly"
            )
        super().__init__(sequence)

    def complement(self) -> "NucleicAcidSequence":
        """Returns the complementary sequence (polymorphism)"""
        complement_table = self._get_complement_table()
        complemented = self._sequence.translate(complement_table)
        return self.__class__(complemented)

    def reverse(self) -> "NucleicAcidSequence":
        """Returns the reversed sequence"""
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """Returns the reverse complementary sequence"""
        return self.complement().reverse()

    @abstractmethod
    def _get_complement_table(self) -> dict:
        """Complementarity table (applies to DNA/RNA in subclasses)"""
        pass

    @abstractmethod
    def _validate_alphabet(self, seq: str) -> bool:
        """Alphabet validation (applies to DNA/RNA in subclasses)"""
        pass


class DNASequence(NucleicAcidSequence):
    """Class for DNA"""

    VALID_LETTERS = set("ATGC")
    """Set of valid DNA letters"""

    def _validate_alphabet(self, seq: str) -> bool:
        """DNA letters validation"""
        if not seq:
            return True
        return all(char in self.VALID_LETTERS for char in seq) and "U" not in seq

    def _get_complement_table(self) -> dict:
        """Make complementary DNA sequence"""
        return str.maketrans("ATGC", "TACG")

    def transcribe(self) -> "RNASequence":
        """Transcribe DNA to RNA (T -> U)."""
        rna_seq = self._sequence.replace("T", "U")
        return RNASequence(rna_seq)


class RNASequence(NucleicAcidSequence):
    """Class for RNA"""

    VALID_LETTERS = set("AUGC")
    """Set of valid RNA letters"""

    def _validate_alphabet(self, seq: str) -> bool:
        """RNA letters validation"""
        if not seq:
            return True
        return all(char in self.VALID_LETTERS for char in seq) and "T" not in seq

    def _get_complement_table(self) -> dict:
        """Make complementary RNA sequence"""
        return str.maketrans("AUGC", "UACG")


class AminoAcidSequence(BiologicalSequence):
    """Class for amino acid (AA) sequence"""

    VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
    """Set of valid AA letters"""

    def _validate_alphabet(self, seq: str) -> bool:
        """AA letters validation"""
        if not seq:
            return True
        return all(char in self.VALID_AMINO_ACIDS for char in seq)

    def count_amino_acids(self) -> int:
        """Count number of AA in protein sequence"""
        return len(self._sequence)


###########  Fastq filtration using biopython


def normalize_bounds(bounds):
    """Normalize bounds to a tuple of floats (min, max)."""
    if isinstance(bounds, (int, float)):
        return (0.0, float(bounds))
    return (float(bounds[0]), float(bounds[1]))


def filter_fastq(
    input_fastq: str,
    output_fastq: str = "filtered/output_fastq.fastq",
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: float = 0.0,
) -> dict:
    """Fastq reads filtration GC, length and quality (biopython)
    Returns:
    dict: {read_id: (sequence_str, List[int])}
        where quality scores are Phred values from Biopython"""
    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)

    filtered_reads = {}
    filtered_records = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        seq_str = str(record.seq).upper()
        length = len(record.seq)
        gc_content = gc_fraction(record.seq) * 100
        qualities = record.letter_annotations.get("phred_quality", [])
        avg_quality = sum(qualities) / len(qualities) if qualities else 0

        if (
            gc_min <= gc_content <= gc_max
            and len_min <= length <= len_max
            and avg_quality >= quality_threshold
        ):
            filtered_reads[record.id] = (seq_str, qualities)
            filtered_records.append(record)

    if output_fastq and filtered_records:
        os.makedirs(os.path.dirname(output_fastq), exist_ok=True)
        SeqIO.write(filtered_records, output_fastq, "fastq")

    return filtered_reads


########### Fasta и Blast


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> str:
    """Converting multi-line Fasta file to a single-line Fasta file"""
    input_path = Path(input_fasta)

    if output_fasta is None:
        output_fasta = input_path.with_stem(input_path.stem + "_oneline")

    with open(input_path, "r") as fin, open(output_fasta, "w") as fout:
        current_header = None
        current_seq = []

        for line in fin:
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith(">"):
                if current_header is not None:
                    fout.write(f">{current_header}\n{''.join(current_seq)}\n")
                current_header = stripped[1:]
                current_seq = []
            else:
                current_seq.append(stripped)

        if current_header is not None:
            fout.write(f">{current_header}\n{''.join(current_seq)}\n")

    return str(output_fasta)


def parse_blast_output(input_file: str, output_file: str):
    """Parses BLAST text output and extracts unique protein descriptions"""
    descriptions = set()
    with open(input_file, "r") as f:
        lines = f.readlines()

    in_results = False
    skip_next = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("Sequences producing significant alignments:"):
            in_results = True
            skip_next = True
            continue

        if not in_results:
            continue

        if skip_next:
            skip_next = False
            continue

        if not stripped or stripped.startswith(
            (
                "Lambda",
                "Effective search",
                "Query #",
                "RID:",
                "Job Title",
                "Program:",
                "Database:",
            )
        ):
            in_results = False
            continue

        if "Scientific" in stripped and "Common" in stripped and "Max" in stripped:
            continue

        parts = [part for part in line.split("  ") if part.strip()]
        if parts:
            description = parts[0].strip()
            if description.startswith("TPA: "):
                description = description[5:]
            descriptions.add(description)

    with open(output_file, "w") as f:
        for desc in sorted(descriptions):
            f.write(desc + "\n")


###########  main

if __name__ == "__main__":
    DATA_DIR = Path("data")
    RESULTS_DIR = Path("filtered")
    RESULTS_DIR.mkdir(exist_ok=True)

    # DNA
    dna = DNASequence("ATGCATGC")
    print(f"Original DNA sequence: {dna}")
    print(f"DNA length (letters): {len(dna)}")
    print(f"Complementary DNA sequence: {dna.complement()}")
    print(f"Reverse DNA sequencing: {dna.reverse()}")
    print(f"Reverse complementary DNA sequence: {dna.reverse_complement()}")
    print(f"Slice [2:6]: {dna[2:6]}")
    print(f"Validation of the DNA alphabet: {dna.check_alphabet()}")

    # DNA to RNA
    rna = dna.transcribe()
    print(f"Original DNA sequence: {dna}")
    print(f"Сomplementary RNA sequence: {rna}")
    print(f"Сhecking the result type: {type(rna).__name__}")

    # RNA
    rna_seq = RNASequence("AUGCAUGC")
    print(f"Original RNA sequence: {rna_seq}")
    print(f"Complementary RNA sequence: {rna_seq.complement()}")
    print(f"Reverse complementary RNA sequence: {rna_seq.reverse_complement()}")

    # AA
    protein = AminoAcidSequence("ACDEFGHIK")
    print(f"Original protein sequence: {protein}")
    print(f"Protein length: {len(protein)}")
    print(f"AA number in protein: {protein.count_amino_acids()}")

    # Check if NucleicAcidSequence cannot be created
    print("\nAttempt to create NucleicAcidSequence:")
    try:
        na = NucleicAcidSequence("ATGC")
    except (TypeError, NotImplementedError) as e:
        print(f"The error was thrown correctly: {e}")

    # Fastq filtration using Biopython
    input_path = "data/example_fastq.fastq"
    if Path(input_path).exists():
        print("\nStart filtration with Biopython...")
        result = filter_fastq(
            input_fastq=input_path,
            gc_bounds=(0, 100),
            length_bounds=(0, 10000),
            quality_threshold=0.0,
        )
        print(f"{len(result)} reads were filtered.")
        if result:
            print(f"FASTA is saved: filtered/output_fastq.fastq")
        else:
            print("Empty result — check input file or parameters.")
    else:
        print(f"FASTQ file is missing: {input_path}")

    # FASTA and BLAST
    input_fasta = DATA_DIR / "example_multiline_fasta.fasta"
    output_fasta = RESULTS_DIR / "oneline.fasta"
    input_file = DATA_DIR / "example_blast_results.txt"
    output_file = RESULTS_DIR / "blast_results.txt"

    if input_fasta.exists():
        out_fasta = convert_multiline_fasta_to_oneline(
            str(input_fasta), str(output_fasta)
        )
        print(f"FASTA is formed and saved: {out_fasta}")
    else:
        print(f"FASTA file is missing: {input_fasta}")

    if input_file.exists():
        parse_blast_output(str(input_file), str(output_file))
        print(f"BLAST rersults are saved: {output_file}")
    else:
        print(f"BLAST file is missing: {input_file}")
