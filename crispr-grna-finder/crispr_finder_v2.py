from pathlib import Path

from Bio.Blast import NCBIWWW, NCBIXML
import typer

from crispr_finder import (
    parse_fasta,
    analyze_sequence,
    filter_candidates,
)

app = typer.Typer(help="CRISPR gRNA анализатор v2: скоринг и off-target")


def calculate_score(sequence: str, gc_content: float, has_repeats: bool) -> float:
    """
    Оценивает качество gRNA по шкале 0-100+.
    """
    score = 100.0

    score -= abs(gc_content - 50) * 1.5

    if has_repeats:
        score -= 20

    last_6 = sequence[-6:]
    gc_at_end = last_6.count("G") + last_6.count("C")
    if gc_at_end > 3:
        score -= (gc_at_end - 3) * 5

    if sequence[0] == "T":
        score -= 10  # U6 плохо работает с T
    elif sequence[0] == "G":
        score += 5  # U6 любит G

    return max(0, score)


def score_candidates(candidates: list[dict]) -> list[dict]:
    """
    Добавляет скор к списку кандидатов.
    """
    for candidate in candidates:
        candidate["score"] = calculate_score(
            candidate["sequence"],
            candidate["gc_content"],
            candidate["has_poly_repeats"],
        )
    return candidates


def check_off_targets(
    sequence: str,
    organism: str = "Homo sapiens",
    max_hits: int = 50,
) -> list[dict]:
    """
    Ищет потенциальные off-target сайты через NCBI BLAST.

    Args:
        sequence: последовательность gRNA (20 нт)
        organism: ограничить поиск организмом
        max_hits: максимум результатов

    Returns:
        Список словарей с информацией о каждом совпадении:
        - title: название последовательности
        - identity: процент идентичности
        - mismatches: количество несовпадений
        - e_value: статистическая значимость
    """
    query = sequence + "NGG"

    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
        sequence=query,
        entrez_query=f'"{organism}"[organism]',
        hitlist_size=max_hits,
        word_size=7,  # маленький word_size для коротких последовательностей
        expect=1000,  # высокий порог, чтобы поймать слабые совпадения
        megablast=False,  # обычный blastn, не megablast
    )

    off_targets = []

    for record in NCBIXML.parse(result_handle):
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                identity = (hsp.identities / hsp.align_length) * 100

                off_targets.append(
                    {
                        "title": alignment.title[:80],
                        "identity": round(identity, 1),
                        "mismatches": hsp.align_length - hsp.identities,
                        "e_value": hsp.expect,
                    }
                )

    result_handle.close()
    return off_targets


def assess_off_target_risk(off_targets: list[dict], target_gene: str) -> dict:
    """
    Классифицирует off-target по уровню риска.

    Args:
        off_targets: результаты BLAST
        target_gene: название целевого гена (для фильтрации)

    Returns:
        Словарь с категоризированными off-target и флагом is_safe
    """

    short_name = target_gene.split()[0].split("_")[0].lower()
    target_patterns = [
        short_name,
        "chemokine receptor 5",
        "cc chemokine receptor",
        "c-c motif chemokine receptor",
        "ccr-5",
        "nonfunctional cc chemokine",
        "cc chemokine re",
    ]

    high_risk = []
    medium_risk = []
    low_risk = []

    for ot in off_targets:
        title_lower = ot["title"].lower()

        if any(pattern in title_lower for pattern in target_patterns):
            continue

        if ot["identity"] >= 95:
            high_risk.append(ot)
        elif ot["identity"] >= 85:
            medium_risk.append(ot)
        else:
            low_risk.append(ot)

    return {
        "high_risk": high_risk,
        "medium_risk": medium_risk,
        "low_risk": low_risk,
        "high_risk_count": len(high_risk),
        "medium_risk_count": len(medium_risk),
        "low_risk_count": len(low_risk),
        "is_safe": len(high_risk) == 0,
    }


@app.command()
def analyze(
    fasta_path: Path = typer.Argument(..., help="Путь к FASTA файлу"),
    gc_min: float = typer.Option(40.0, help="Минимальный GC-состав (%)"),
    gc_max: float = typer.Option(60.0, help="Максимальный GC-состав (%)"),
    top_n: int = typer.Option(10, help="Показать топ N кандидатов"),
    check_offtarget: bool = typer.Option(
        False, "--check-offtarget", help="Проверить off-target через BLAST (медленно!)"
    ),
):
    """
    Анализ gRNA кандидатов с скорингом.
    """
    sequences = parse_fasta(fasta_path)

    for name, sequence in sequences.items():
        typer.echo(f"\n{'=' * 60}")
        typer.echo(f"Ген: {name}")
        typer.echo(f"Длина: {len(sequence)} нуклеотидов")
        typer.echo(f"{'=' * 60}")

        candidates = analyze_sequence(sequence)
        filtered = filter_candidates(candidates, gc_min, gc_max, allow_repeats=False)

        scored = score_candidates(filtered)
        scored.sort(key=lambda x: x["score"], reverse=True)

        if not scored:
            typer.echo("Кандидатов не найдено")
            continue

        typer.echo(f"\nНайдено кандидатов: {len(scored)}")
        typer.echo(f"Топ-{min(top_n, len(scored))} по скору:\n")

        for i, c in enumerate(scored[:top_n], 1):
            strand = "+" if c["is_strand"] else "-"
            typer.echo(
                f"{i:2}. {c['sequence']} | "
                f"score: {c['score']:.1f} | "
                f"GC: {c['gc_content']:.1f}% | "
                f"strand: {strand}"
            )

            if check_offtarget:
                typer.echo("    Проверяю off-target (займёт ~30-60 сек)...")
                try:
                    off_targets = check_off_targets(c["sequence"])
                    risk = assess_off_target_risk(off_targets, name)

                    if risk["is_safe"]:
                        typer.echo("    ✓ Высокий риск off-target не обнаружен")
                    else:
                        typer.echo(
                            f"    ⚠ Риск: {risk['high_risk_count']} высокий, "
                            f"{risk['medium_risk_count']} средний"
                        )

                        for ot in risk["high_risk"][:3]:
                            typer.echo(f"      - {ot['title']}: {ot['identity']}%")

                except Exception as e:
                    typer.echo(f"    ✗ Ошибка BLAST: {e}")


@app.command()
def score(
    sequence: str = typer.Argument(
        ..., help="Последовательность gRNA (20 нуклеотидов)"
    ),
):
    """
    Рассчитать скор для одной gRNA.
    """
    if len(sequence) != 20:
        typer.echo(f"Ошибка: gRNA должна быть 20 нуклеотидов, получено {len(sequence)}")
        raise typer.Exit(code=1)

    sequence = sequence.upper()

    gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100

    has_repeats = any(base * 4 in sequence for base in "ACGT")

    score_value = calculate_score(sequence, gc, has_repeats)

    typer.echo(f"Последовательность: {sequence}")
    typer.echo(f"GC-состав: {gc:.1f}%")
    typer.echo(f"Повторы: {'да' if has_repeats else 'нет'}")
    typer.echo(f"Скор: {score_value:.1f}")


@app.command()
def blast(
    sequence: str = typer.Argument(
        ..., help="Последовательность gRNA (20 нуклеотидов)"
    ),
    organism: str = typer.Option("Homo sapiens", help="Организм для поиска"),
):
    """
    Проверить одну gRNA на off-target через BLAST.
    """
    if len(sequence) != 20:
        typer.echo(f"Ошибка: gRNA должна быть 20 нуклеотидов, получено {len(sequence)}")
        raise typer.Exit(code=1)

    sequence = sequence.upper()

    typer.echo(f"Ищу off-target для {sequence} в {organism}...")
    typer.echo("(это займёт 30-60 секунд)\n")

    try:
        off_targets = check_off_targets(sequence, organism)

        if not off_targets:
            typer.echo("Совпадений не найдено")
            return

        typer.echo(f"Найдено совпадений: {len(off_targets)}\n")

        for i, ot in enumerate(off_targets[:10], 1):
            typer.echo(f"{i:2}. {ot['title']}")
            typer.echo(
                f"    Identity: {ot['identity']}% | Mismatches: {ot['mismatches']}"
            )

    except Exception as e:
        typer.echo(f"Ошибка BLAST: {e}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
