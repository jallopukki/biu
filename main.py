from pathlib import Path
from typing import List, Union, Dict, Tuple
import sys

sys.path.append("./modules")


from modules.run_dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
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
        str: if single sequence and List[str] if multiple sequences.
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
    resperations[operation]
    result = [func(seq) for seq in seqs]

    return result[0] if len(result) == 1 else result


def filter_fastq(
    input_fastq: str,
    output_fastq: str = None,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0.0,
) -> dict:
    """
    Filter FASTQ reads by GC content, length, and quality.
    Reads from input_fastq, applies filters, optionally writes to output_fastq.

    Returns:
        dict: filtered reads as {read_name: (sequence, quality)}
    """
    filtered_reads = _filter_fastq(
        input_fastq=input_fastq,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )

    if output_fastq is not None:
        write_fastq(output_fastq, filtered_reads)

    return filtered_reads

    # этот блок мне сделал жипити, спас меня от мучений запуска.
    # В терминале классно выглядит, можно оставить, а!? (на рус) пжлст


if __name__ == "__main__":
    print("Запуск фильтрации...")
    input_path = "data/example_fastq.fastq"
    output_path = "filtered/output_fastq.fastq"

    result = filter_fastq(
        input_fastq=input_path,
        output_fastq=output_path,
        gc_bounds=(0, 100),
        length_bounds=(0, 10000),
        quality_threshold=0.0,
    )

    print(f"✅ Отфильтровано {len(result)} ридов.")
    if result:
        print(f"Файл сохранён: {output_path}")
    else:
        print("⚠️  Результат пуст — проверьте входной файл или параметры.")
