"""Unit tests for FASTQ filtering utilities and core logic."""
import sys
from pathlib import Path

# Add project root to sys.path to resolve 'modules' imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest
from modules.filter_fastq import (
    mean_quality,
    calculate_gc,
    normalize_bounds,
    read_fastq,
    write_fastq,
    filter_fastq,
)


@pytest.fixture
def sample_fastq(tmp_path):
    """Create a temporary FASTQ file with two reads of different quality."""
    p = tmp_path / "sample.fastq"
    # 'I' = Phred 40, '?' = Phred 30
    p.write_text("@r1\nATGC\n+\nIIII\n@r2\nAAAA\n+\n????\n")
    return str(p)


class TestCoreUtils:
    """Tests for standalone utility functions."""

    def test_mean_quality(self):
        """Verify Phred+33 quality score calculation."""
        assert mean_quality("IIII") == 40.0

    def test_gc_calculation(self):
        """Check GC content percentage (50% for 'ATGC')."""
        assert calculate_gc("ATGC") == 50.0

    def test_bounds_normalization(self):
        """Validate conversion of bounds to (min, max) float tuples."""
        assert normalize_bounds((0, 100)) == (0.0, 100.0)
        assert normalize_bounds(50) == (0.0, 50.0)


class TestFilterLogic:
    """Tests for read filtering based on quality, length, and content."""

    def test_filters_by_quality(self, sample_fastq):
        """Ensure reads below the quality threshold are excluded."""
        # Threshold 35.0 filters out '?' (Phred 30)
        res = filter_fastq(sample_fastq, quality_threshold=35.0)
        assert "r1" in res and "r2" not in res

    def test_filters_by_length(self, sample_fastq):
        """Verify length bounds exclude reads outside the range."""
        res = filter_fastq(sample_fastq, length_bounds=(5, 10))
        assert len(res) == 0

    def test_empty_sequence_skip(self, tmp_path):
        """Check that empty sequences are safely ignored."""
        p = tmp_path / "empty.fastq"
        p.write_text("@empty\n\n+\n\n")
        res = filter_fastq(str(p))
        assert len(res) == 0


class TestIOAndErrors:
    """Tests for file I/O operations and error handling."""

    def test_write_read_roundtrip(self, tmp_path):
        """Verify data integrity after writing and reading a FASTQ file."""
        data = {"seq1": ("ATGC", "IIII")}
        p = tmp_path / "rw.fastq"
        write_fastq(str(p), data)
        assert read_fastq(str(p)) == data

    def test_missing_file_error(self):
        """Ensure FileNotFoundError is raised for invalid input paths."""
        with pytest.raises(FileNotFoundError):
            filter_fastq("nonexistent.fastq")
