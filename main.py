"""CLI entry point for FASTQ filtering and DNA/RNA sequence tools."""
import argparse
import logging
import sys
from pathlib import Path

sys.path.append("./modules")

from modules.run_dna_rna_tools import (
    transcribe, reverse, complement, reverse_complement, is_nucleic_acid,
)
from modules.filter_fastq import (
    filter_fastq as _filter_fastq,
    write_fastq,
)

DATA_DIR = Path("data")
FILTERED_DIR = Path("filtered")


def run_dna_rna_tools(*args: str):
    """Apply DNA/RNA operations to one or more sequences."""
    *seqs, operation = args
    ops = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    result = [ops[operation](seq) for seq in seqs]
    return result[0] if len(result) == 1 else result


def filter_fastq(
    input_fastq: str,
    output_fastq: str | None = "filtered/output_fastq.fastq",
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0.0,
) -> dict:
    """Filter FASTQ reads and optionally write results to file.

    Returns dict of filtered reads: {read_name: (sequence, quality)}.
    """
    filtered = _filter_fastq(input_fastq, gc_bounds, length_bounds, quality_threshold)
    if output_fastq:
        Path(output_fastq).parent.mkdir(parents=True, exist_ok=True)
        write_fastq(output_fastq, filtered)
    return filtered


def setup_logger(log_file: str = "filter.log"):
    """Configure logging to file and stdout with INFO level."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def parse_args():
    """Parse command-line arguments for FASTQ filtering."""
    p = argparse.ArgumentParser(description="FastQ Filter CLI")
    p.add_argument("input", help="Input FASTQ file")
    p.add_argument("-o", "--output", default="filtered/output_fastq.fastq", help="Output FASTQ path")
    p.add_argument("--gc-min", type=float, default=0.0, help="Min GC content (percent)")
    p.add_argument("--gc-max", type=float, default=100.0, help="Max GC content (percent)")
    p.add_argument("--len-min", type=int, default=0, help="Min read length")
    p.add_argument("--len-max", type=int, default=2**32, help="Max read length")
    p.add_argument("--quality", type=float, default=0.0, help="Min average Phred quality")
    return p.parse_args()


def main():
    """CLI entry point: parse args, run filtration, handle errors."""
    setup_logger()
    args = parse_args()

    if not Path(args.input).exists():
        logging.error(f"Input file not found: {args.input}")
        sys.exit(1)

    logging.info(f"Starting filtration: {args.input}")
    try:
        result = filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=(args.gc_min, args.gc_max),
            length_bounds=(args.len_min, args.len_max),
            quality_threshold=args.quality,
        )
        logging.info(f"Done. {len(result)} reads saved to {args.output}")
    except Exception as e:
        logging.error(f"Filtration failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
