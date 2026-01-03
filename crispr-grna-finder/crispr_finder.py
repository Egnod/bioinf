from pathlib import Path

import typer

app = typer.Typer()


def find_ngg_pam(sequence: str) -> list[int]:
    candidates_positions = []

    last_idx = -1

    while True:
        last_idx = sequence.find("GG", last_idx + 1)

        if last_idx > 0:
            candidates_positions.append(last_idx - 1)
        elif last_idx == -1:
            break

    return candidates_positions


def extract_grna_candidates(sequence: str, pam_positions: list[int]) -> list[dict]:
    grna_candidates = []

    for pam_pos in pam_positions:
        if pam_pos >= 20:
            grna_seq = sequence[pam_pos - 20 : pam_pos]
            grna_candidates.append({"sequence": grna_seq, "pam_position": pam_pos})

    return grna_candidates


def calculate_gc_content(sequence: str) -> float:
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0.0


def has_poly_repeats(sequence: str, min_length: int = 4) -> bool:
    for base in "ACGT":
        if base * min_length in sequence:
            return True
    return False


def _analyze_sequence(sequence: str, is_strand: bool) -> list[dict]:
    pam_positions = find_ngg_pam(sequence)
    grna_candidates = extract_grna_candidates(sequence, pam_positions)

    analyzed_candidates = []
    for candidate in grna_candidates:
        gc_content = calculate_gc_content(candidate["sequence"])
        poly_repeats = has_poly_repeats(candidate["sequence"])

        analyzed_candidates.append(
            {
                "sequence": candidate["sequence"],
                "pam_position": candidate["pam_position"],
                "gc_content": gc_content,
                "has_poly_repeats": poly_repeats,
                "is_strand": is_strand,
            }
        )

    return analyzed_candidates


def analyze_sequence(sequence: str) -> list[dict]:
    candidates = []

    candidates.extend(_analyze_sequence(sequence, is_strand=True))
    rev_comp_sequence = reverse_complement(sequence)
    candidates.extend(_analyze_sequence(rev_comp_sequence, is_strand=False))

    return candidates


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.translate(complement)[::-1]


def filter_candidates(
    candidates: list[dict],
    gc_min: float = 40.0,
    gc_max: float = 60.0,
    allow_repeats: bool = False,
) -> list[dict]:
    filtered = []
    for candidate in candidates:
        if gc_min <= candidate["gc_content"] <= gc_max:
            if allow_repeats or not candidate["has_poly_repeats"]:
                filtered.append(candidate)
    return filtered


def parse_fasta(file_path: Path) -> dict[str, str]:
    """Парсит FASTA файл, возвращает {название: последовательность}"""
    if not file_path.exists():
        typer.echo(f"Файл не найден: {file_path}")
        raise typer.Exit(code=1)

    raw_data = file_path.read_text()
    results = {}

    sequence = ""
    name = ""

    for line in raw_data.splitlines():
        if line.startswith(">"):
            if sequence and name:
                results[name] = sequence

            name = line[1:].strip()
            sequence = ""
        else:
            sequence += line.strip().upper()

    if sequence and name:
        results[name] = sequence

    return results


@app.command()
def find(
    fasta_path: Path = typer.Argument(..., help="Путь к FASTA файлу"),
    gc_min: float = typer.Option(40.0, help="Минимальный GC-состав"),
    gc_max: float = typer.Option(60.0, help="Максимальный GC-состав"),
    allow_repeats: bool = typer.Option(False, help="Разрешить повторы"),
):
    """Поиск CRISPR gRNA кандидатов в последовательности"""

    sequences = parse_fasta(fasta_path)

    for name, sequence in sequences.items():
        typer.echo(f"\n=== {name} ===")

        candidates = analyze_sequence(sequence)
        filtered = filter_candidates(candidates, gc_min, gc_max, allow_repeats)

        if not filtered:
            typer.echo("Кандидатов не найдено")
            continue

        for c in filtered:
            strand = "+" if c["is_strand"] else "-"
            typer.echo(
                f"  {c['sequence']} | pos: {c['pam_position']} | "
                f"GC: {c['gc_content']:.1f}% | strand: {strand}"
            )


if __name__ == "__main__":
    app()
