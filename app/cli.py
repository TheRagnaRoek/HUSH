"""HUSH command-line interface.

CLI-first entry point for the pipeline. Currently exposes the Input Layer via
'hush input'; siRNA design, MSA/allele, and shRNA conversion commands are
added as those stages come online.

Run 'hush --help' (after 'pip install -e .') or 'python -m app.cli --help'.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from core.errors import HushError
from core.inputs import build_design_input, load_sequences
from core.sequence import BioSequence, DesignInput

app = typer.Typer(
    add_completion=False,
    help="HUSH - Human genome Universal Silencing Hairpin-designer.",
    no_args_is_help=True,
)


def _echo_sequence(seq: BioSequence, index: int | None = None) -> None:
    prefix = f"  [{index}] " if index is not None else "  "
    header = f"{prefix}{seq.id or '<unnamed>'}"
    typer.echo(typer.style(header, bold=True))
    if seq.description:
        typer.echo(f"      {seq.description}")
    line = f"      type={seq.moltype.value}  length={seq.length}"
    if seq.is_nucleotide and seq.length:
        line += f"  GC={seq.gc_content:.1%}"
    typer.echo(line)
    preview = seq.seq if seq.length <= 60 else f"{seq.seq[:57]}..."
    typer.echo(f"      {preview}")


def _summarize_design_input(design: DesignInput) -> None:
    typer.echo(typer.style("Target", fg=typer.colors.CYAN, bold=True))
    _echo_sequence(design.target)

    if design.has_panel:
        assert design.panel is not None
        typer.echo("")
        typer.echo(
            typer.style(
                f"Panel: {design.panel.name}  ({len(design.panel)} sequences)",
                fg=typer.colors.CYAN,
                bold=True,
            )
        )
        uniform = "uniform" if design.panel.is_uniform else "MIXED types"
        typer.echo(
            f"  molecule types: "
            f"{', '.join(sorted(t.value for t in design.panel.moltypes))} "
            f"({uniform})"
        )
        for i, seq in enumerate(design.panel):
            _echo_sequence(seq, index=i)

    # Foreshadow the downstream requirement without failing the inspection.
    typer.echo("")
    if design.target.is_nucleotide:
        typer.echo(
            typer.style(
                "[OK] Target is a nucleotide sequence - ready for siRNA design.",
                fg=typer.colors.GREEN,
            )
        )
    else:
        typer.echo(
            typer.style(
                f"[!] Target is {design.target.moltype.value}; siRNA/shRNA design "
                "needs the coding transcript (mRNA/CDS). Supply the nucleotide "
                "sequence before designing.",
                fg=typer.colors.YELLOW,
            )
        )


@app.command("input")
def input_cmd(
    target: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        readable=True,
        help="Target sequence file (FASTA or PDB).",
    ),
    panel: Optional[Path] = typer.Option(
        None,
        "--panel",
        "-p",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Optional related-sequence panel (e.g. an HLA allele family FASTA).",
    ),
    target_format: Optional[str] = typer.Option(
        None, "--target-format", help="Override target format: fasta | pdb | cif."
    ),
    panel_format: Optional[str] = typer.Option(
        None, "--panel-format", help="Override panel format: fasta | pdb | cif."
    ),
) -> None:
    """Load and inspect target (and optional panel) sequences.

    This exercises the Input Layer: it parses the file(s), validates residues,
    auto-detects molecule type, and reports what the downstream pipeline will
    receive.
    """
    try:
        design = build_design_input(
            target,
            panel,
            target_fmt=target_format,
            panel_fmt=panel_format,
        )
    except HushError as exc:
        typer.secho(f"Input error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    _summarize_design_input(design)


@app.command("inspect")
def inspect_cmd(
    path: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="Any FASTA/PDB file to list sequences from.",
    ),
    fmt: Optional[str] = typer.Option(
        None, "--format", "-f", help="Override format: fasta | pdb | cif."
    ),
) -> None:
    """List every sequence in a file (no target/panel roles assigned)."""
    try:
        seqs = load_sequences(path, fmt=fmt)
    except HushError as exc:
        typer.secho(f"Input error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    typer.echo(f"{len(seqs)} sequence(s) in {path.name}:")
    for i, seq in enumerate(seqs):
        _echo_sequence(seq, index=i)


if __name__ == "__main__":
    app()
