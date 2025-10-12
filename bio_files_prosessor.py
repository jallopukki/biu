from pathlib import Path

DATA_DIR = Path("data")
RESULTS_DIR = Path("filtered")

input_fasta = DATA_DIR / "example_multiline_fasta.fasta"
output_fasta = RESULTS_DIR / "oneline.fasta"

input_file = DATA_DIR / "example_blast_results.txt"
output_file = RESULTS_DIR / "blast_results.txt"


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> str:
    """
    Converts a multi-line FASTA file to a one-line-per-sequence FASTA file.

    Parameters:
        input_fasta (str): Path to the input FASTA file (may have multi-line sequences).
        output_fasta (str, optional): Path to the output file.
                                      If not provided, uses '<stem>_oneline.fasta'.

    Returns:
        str: Path to the output file.
    """
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
                # write previous seq
                if current_header is not None:
                    fout.write(f">{current_header}\n{''.join(current_seq)}\n")
                current_header = stripped[1:]  # убираем '>'
                current_seq = []
            else:
                current_seq.append(stripped)

        # write last seq
        if current_header is not None:
            fout.write(f">{current_header}\n{''.join(current_seq)}\n")

    return str(output_fasta)


def parse_blast_output(input_file: str, output_file: str):
    """
    Parses BLAST text output (format 0) and extracts unique protein descriptions
    from the 'Sequences producing significant alignments' section.
    """
    descriptions = set()
    with open(input_file, "r") as f:
        lines = f.readlines()

    in_results = False
    skip_next = False  # skip str with heads

    for line in lines:
        stripped = line.strip()

        # start of section with results
        if stripped.startswith("Sequences producing significant alignments:"):
            in_results = True
            skip_next = True  # следующая строка — заголовок таблицы
            continue

        if not in_results:
            continue

        # skip head of the table
        if skip_next:
            skip_next = False
            continue

        # exit from the section
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

        # skip something like heads
        if "Scientific" in stripped and "Common" in stripped and "Max" in stripped:
            continue

        # this is start of required data
        # split by 2 spaces
        parts = [part for part in line.split("  ") if part.strip()]
        if parts:
            description = parts[0].strip()
            # del "TPA: on start if it is
            if description.startswith("TPA: "):
                description = description[5:]
            descriptions.add(description)

    # save sorted unique descriptions
    with open(output_file, "w") as f:
        for desc in sorted(descriptions):
            f.write(desc + "\n")


if __name__ == "__main__":
    RESULTS_DIR.mkdir(exist_ok=True)

    # а здесь я уже сама-сама))
    out_fasta = convert_multiline_fasta_to_oneline(str(input_fasta), str(output_fasta))
    print("Запускаю конвертацию одна последовательность - одна строка...")
    print(f"✅ FASTA сформирован и сохранен в: {out_fasta}")

    if input_file.exists():
        parse_blast_output(str(input_file), str(output_file))
        print("Отбираю первые строки из столбца Description...")
        print("Сохраняю набор полученных белков в новый файл в алфавитном порядке...")
        print(f"✅ результаты BLAST сохранены в: {output_file}")
    else:
        print(f"⚠️  отсутствует BLAST файл: {input_file}")
