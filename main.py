from modules.run_dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
)
from modules.filter_fastq import mean_quality


def run_dna_rna_tools(*args):
    *seqs, operation = args
    result = []
    """Tools for working with RNA/DNA sequences."""

    for seq in seqs:
        if operation == "is_nucleic_acid":
            """Validation RNA or DNA. Return bool result"""
            valid_letters = set("ATGCUatgcu")
            return all(char in valid_letters for char in seq) and not (
                "T" in seq.upper() and "U" in seq.upper()
            )
        elif operation == "transcribe":
            result.append(transcribe(seq))
        elif operation == "reverse":
            result.append(reverse(seq))
        elif operation == "complement":
            result.append(complement(seq))
        elif operation == "reverse_complement":
            result.append(complement(reverse(seq)))

    if len(result) == 1:
        return result[0]
    else:
        return result


def filter_fastq(
    seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0.0
):
    """Filters FASTQ sequences based on GC content, length, and quality."""

    def normalize_bounds(bounds):
        if isinstance(bounds, (int, float)):
            return (0.0, float(bounds))
        return (float(bounds[0]), float(bounds[1]))

    gc_min, gc_max = normalize_bounds(gc_bounds)
    len_min, len_max = normalize_bounds(length_bounds)
    filtered = {}

    for name, (seq, quality) in seqs.items():
        if not seq:
            continue

        """Calculate GC (%)."""
        gc_percentage = (
            (seq.upper().count("G") + seq.upper().count("C")) / len(seq) * 100
        )
        length = len(seq)
        avg_qual = mean_quality(quality)

        if (
            gc_min <= gc_percentage <= gc_max
            and len_min <= length <= len_max
            and avg_qual >= quality_threshold
        ):
            filtered[name] = (seq, quality)

    return filtered
