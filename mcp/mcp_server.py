"""
This module exposes a basic MCP server making
hush design
hush convert
hush reference
available to AI agents. This MCP server exposes all possible parameters and executes
commands as subprocesses using the HUSH executable.
"""
import os
import subprocess
from pathlib import Path
from typing import Optional
from mcp.server.mcpserver import MCPServer

mcp = MCPServer("HUSH siRNA Designer")
HUSH_EXECUTABLE = "hush"

@mcp.tool()
def design_sirna(
    target: str,
    length: int = 19,
    overhang: str = "natural",
    rules: str = "all",
    mode: str = "any",
    seed: bool = True,
    seed_tm: Optional[float] = None,
    offtarget: bool = True,
    offtarget_db: Optional[str] = None,
    offtarget_min_mismatch: Optional[int] = None,
    offtarget_max_mismatch: Optional[int] = None,
    offtarget_self_accession: Optional[str] = None,
    offtarget_self_fasta: Optional[str] = None,
    offtarget_resolve_self: bool = True,
    offtarget_self_coverage: Optional[float] = None,
    gc: bool = False,
    gc_min: float = 30.0,
    gc_max: float = 52.0,
    allow_homopolymers: bool = False,
    max_run: Optional[int] = None,
    asymmetry: bool = False,
    hairpin: bool = False,
    hairpin_tm: Optional[float] = None,
    asym_window: Optional[int] = None,
    asym_temp: Optional[float] = None,
    params: str = "xia",
    panel: Optional[str] = None,
    panel_format: Optional[str] = None,
    pan_allele: bool = False,
    allele_specific: Optional[str] = None,
    target_specific: bool = False,
    hit: Optional[str] = None,
    spare: Optional[str] = None,
    min_spare_mismatches: Optional[int] = None,
    sort: str = "auto",
    explain: bool = False,
    limit: int = 20,
    reynolds_threshold: Optional[int] = None,
    amarzguioui_threshold: Optional[int] = None,
    output: Optional[str] = None
) -> str:
    """Design siRNA candidates against a target mRNA/CDS.

    The target must be a FASTA or PDB file accessible to the server.
    Optionally supply an allele-family panel for allele-specific,
    pan-allele, or custom hit/spare targeting.

    Returns the formatted candidate table and any warnings thrown by HUSH.
    """

    target_path = Path(target).expanduser().resolve()

    if not target_path.is_file():
        return f"Target file not found: {target}"

    cmd = [
        HUSH_EXECUTABLE, "design", str(target_path),
        "--length", str(length),
        "--overhang", overhang,
        "--rules", rules,
        "--mode", mode,
        "--seed" if seed else "--no-seed",
        "--offtarget" if offtarget else "--no-offtarget",
        "--offtarget-resolve-self" if offtarget_resolve_self
        else "--no-offtarget-resolve-self",
        "--gc" if gc else "--no-gc",
        "--allow-homopolymers" if allow_homopolymers
        else "--forbid-homopolymers",
        "--asymmetry" if asymmetry else "--no-asymmetry",
        "--hairpin" if hairpin else "--no-hairpin",
        "--params", params,
        "--sort", sort,
        "--limit", str(limit),
    ]

    optional_args = {
        "--seed-tm": seed_tm,
        "--offtarget-db": offtarget_db,
        "--offtarget-min-mismatch": offtarget_min_mismatch,
        "--offtarget-max-mismatch": offtarget_max_mismatch,
        "--offtarget-self-accession": offtarget_self_accession,
        "--offtarget-self-fasta": offtarget_self_fasta,
        "--offtarget-self-coverage": offtarget_self_coverage,
        "--gc-min": gc_min,
        "--gc-max": gc_max,
        "--max-run": max_run,
        "--hairpin-tm": hairpin_tm,
        "--asym-window": asym_window,
        "--asym-temp": asym_temp,
        "--panel": panel,
        "--panel-format": panel_format,
        "--allele-specific": allele_specific,
        "--hit": hit,
        "--spare": spare,
        "--min-spare-mismatches": min_spare_mismatches,
        "--reynolds-threshold": reynolds_threshold,
        "--amarzguioui-threshold": amarzguioui_threshold,
    }

    for flag, value in optional_args.items():
        if value is not None:
            cmd.extend([flag, str(value)])

    # Handle optional flags
    if pan_allele:
        cmd.append("--pan-allele")

    if target_specific:
        cmd.append("--target-specific")

    if explain:
        cmd.append("--explain")

    # Write to file if target is specified
    if output is not None:
        cmd.extend([
            "--output",
            str(Path(output).expanduser().resolve())
        ])

    try:
        # Run with dry formatting <-> LLM
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
            env={**__import__("os").environ, "NO_COLOR": "1"},
        )

        output = "\n".join(
            part for part in [result.stdout, result.stderr] if part
        )

        if result.returncode != 0:
            return f"Design failed (exit {result.returncode}):\n{output}"

        return output or "Design completed without output."

    except subprocess.TimeoutExpired:
        return "Design timed out after 300 seconds."

    except OSError as exc:
        return f"Could not execute hush: {exc}"



@mcp.tool()
def convert_shrna(
    selected: str,
    loop: str = "jensen",
    terminator: bool = True,
    terminator_length: Optional[int] = None,
    drop_premature_terminators: bool = True,
    limit: int = 20,
    output: Optional[str] = None,
) -> str:
    """Convert selected siRNA candidates into stem-loop shRNA constructs.

    Reads a TSV file produced by 'hush design -o' and converts
    each selected siRNA candidate into a sense-loop-antisense
    shRNA hairpin.

    Parameters:
        selected: Path to the selected-candidate TSV file.
        loop: Loop scaffold: jensen (default), psuper, or trc.
        terminator: Append a Pol III poly-T terminator.
        terminator_length: Length of the poly-T terminator.
            Uses the Hush default when omitted.
        drop_premature_terminators: Exclude constructs containing
            premature Pol III terminator sequences.
        limit: Maximum number of constructs to display.
            Use 0 to display all.
        output: Optional path to save all constructs as TSV.

    Returns:
        Formatted shRNA constructs with DNA oligonucleotide
        sequences, folding validation, and QC information.
        The oligos are vector-agnostic and do not include
        cloning-specific sticky ends.
    """

    selected_path = Path(selected).expanduser().resolve()

    if not selected_path.is_file():
        return f"Selected-candidate TSV not found: {selected}"

    if loop not in ("jensen", "psuper", "trc"):
        return f"Invalid loop scaffold: {loop}"

    cmd = [
        HUSH_EXECUTABLE, "convert", str(selected_path),
        "--loop", loop,
        "--terminator" if terminator else "--no-terminator",
        (
            "--drop-premature-terminators"
            if drop_premature_terminators
            else "--keep-premature-terminators"
        ),
        "--limit", str(limit),
    ]

    if terminator_length is not None:
        cmd.extend(["--terminator-length", str(terminator_length)])

    if output is not None:
        output_path = Path(output).expanduser().resolve()
        cmd.extend(["--output", str(output_path)])

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
            env={**__import__("os").environ, "NO_COLOR": "1"},
        )

        output_text = "\n".join(
            part for part in [result.stdout, result.stderr] if part
        )

        if result.returncode != 0:
            return (
                f"Conversion failed (exit {result.returncode}):\n"
                f"{output_text}"
            )

        return output_text or "Conversion completed without output."

    except subprocess.TimeoutExpired:
        return "Conversion timed out after 300 seconds."

    except OSError as exc:
        return f"Could not execute hush: {exc}"



@mcp.tool()
def build_reference(
    outdir: Optional[str] = None,
    base_url: Optional[str] = None,
    index_name: Optional[str] = None,
    keep_fasta: bool = False,
) -> str:
    """Download the human RefSeq transcriptome and build a Bowtie1 index.

    Builds the reference database used by HUSH for transcriptome-wide
    siRNA off-target checking.

    This operation downloads several gigabytes of data and may take
    considerable time. Requires Bowtie1 to be installed.

    Parameters:
        outdir: Directory for the reference index and manifest.
            Defaults to HUSH_REFERENCE_DIR or data/reference.
        base_url: RefSeq mRNA_Prot directory URL.
            Can point to an archived release for reproducibility.
        index_name: Basename for the Bowtie index files.
        keep_fasta: Retain the concatenated transcriptome FASTA
            after building the index.

    Returns:
        Build output, including the location of the reference index.
        Once built, HUSH automatically uses the index for
        off-target filtering during siRNA design.
    """

    cmd = [HUSH_EXECUTABLE, "reference", "build"]

    if outdir is not None:
        outdir_path = Path(outdir).expanduser().resolve()
        cmd.extend(["--outdir", str(outdir_path)])

    if base_url is not None:
        cmd.extend(["--base-url", base_url])

    if index_name is not None:
        cmd.extend(["--index-name", index_name])

    cmd.append("--keep-fasta" if keep_fasta else "--no-keep-fasta")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,
            check=False,
            env={**__import__("os").environ, "NO_COLOR": "1"},
        )

        output = "\n".join(
            part for part in [result.stdout, result.stderr] if part
        )

        if result.returncode != 0:
            return (
                f"Reference build failed (exit {result.returncode}):\n"
                f"{output}"
            )

        return output or "Reference index built successfully."

    except subprocess.TimeoutExpired:
        return (
            "Reference build exceeded the 3600-second timeout. "
            "Check the output directory for partial downloads or "
            "incomplete index files."
        )

    except OSError as exc:
        return f"Could not execute HUSH reference builder: {exc}"



if __name__ == "__main__":
    HUSH_EXECUTABLE = os.environ.get('HUSH_EXECUTABLE', 'hush')
    mcp.run(transport="stdio")
