from pathlib import Path
from typing import List, Union, Dict, Tuple
import sys
import os


sys.path.append("./modules")

from modules.run_dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
    is_nucleic_acid,
)

from modules.filter_fastq import (
    filter_fastq as _filter_fastq,
    write_fastq,
)

DATA_DIR = Path("data")
FILTERED_DIR = Path("filtered")
DEFAULT_INPUT = DATA_DIR / "example_fastq.fastq"
DEFAULT_OUTPUT = FILTERED_DIR / "output_fastq.fastq"


def run_dna_rna_tools(*args: str) -> Union[bool, str, List[Union[str, bool]]]:
    """
    Tools for working with RNA/DNA sequences.

    Args:
        *args: One or more sequences, followed by an operation name.
               Supported operations: 'is_nucleic_acid', 'transcribe',
               'reverse', 'complement', 'reverse_complement'.

    Returns:
        bool: if operation is 'is_nucleic_acid' and one sequence given.
        str or List[str/bool]: result for one or multiple sequences.
    """
    *seqs, operation = args

    operations = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    func = operations[operation]
    result = [func(seq) for seq in seqs]

    return result[0] if len(result) == 1 else result


def filter_fastq(
    input_fastq: str,
    output_fastq: str = "filtered/output_fastq.fastq",
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0.0,
) -> dict:
    """
    Filter FASTQ reads by GC content, length, and quality.
    Saves to output_fastq (default: 'filtered/output_fastq.fastq').
    Returns dict of filtered reads: {read_name: (sequence, quality)},
    where read_name includes the '@' symbol.
    """
    filtered_reads = _filter_fastq(
        input_fastq=input_fastq,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )

    if output_fastq is not None:
        os.makedirs(os.path.dirname(output_fastq), exist_ok=True)
        write_fastq(output_fastq, filtered_reads)

    return filtered_reads


if __name__ == "__main__":
    print("Start filtration...")
    input_path = "data/example_fastq.fastq"

    result = filter_fastq(
        input_fastq=input_path,
        gc_bounds=(0, 100),
        length_bounds=(0, 10000),
        quality_threshold=0.0,
    )

    print(f"✅ {len(result)} reads were filtered.")
    if result:
        print(f"File saved: filtered/output_fastq.fastq")
    else:
        print("⚠️  Empty result — check input file or parameters.")
