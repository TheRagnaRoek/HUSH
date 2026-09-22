"""The interactive HUSH shell (``hush shell``).

A persistent instance of HUSH: load the target (and optional panel) once, tune
the design criteria with 'set', and run 'design'/'convert' as often as you like
without re-parsing or re-loading. State lives in a single HushSession, so the
panel alignment is computed once and reused across designs.

The shell speaks the same option words as the one-shot CLI (--no-seed, --gc,
--allele-specific, ...) and renders results through the very same functions
(app.cli.render_design_result / render_convert_result), so a design in the shell
reads identically to 'hush design'. This module is presentation only; all the
pipeline work is delegated to the session.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import typer

from app.cli import (
    _summarize_design_input,
    _write_convert_tsv,
    _write_design_tsv,
    render_convert_result,
    render_design_result,
    render_discrimination_detail,
)
from app.parsing import (
    parse_mode,
    parse_overhang,
    parse_params,
    parse_rules,
    split_ids,
)
from core.errors import HushError
from service.params import AlleleTargetingSpec, DesignParams
from service.session import HushSession
from shrna.convert import DEFAULT_TERMINATOR_LENGTH, ShrnaOptions
from shrna.loops import DEFAULT_LOOP, LOOPS, get_loop


# --------------------------------------------------------------------------- #
# flag parsing (mirrors the one-shot CLI's design/convert options)
# --------------------------------------------------------------------------- #


class _ShellArgError(ValueError):
    """A flag-parsing failure inside the shell (never exits the process)."""


class _ShellArgParser(argparse.ArgumentParser):
    """An argparse parser that raises instead of calling sys.exit.

    argparse's default is to print usage and terminate the interpreter; inside a
    REPL that would kill the whole session, so both exits are turned into a
    catchable error the command loop reports as 'error: ...'.
    """

    def error(self, message: str) -> None:  # type: ignore[override]
        raise _ShellArgError(message)

    def exit(self, status: int = 0, message: str | None = None) -> None:  # type: ignore[override]
        if message:
            raise _ShellArgError(message.strip())
        raise _ShellArgError("")


def _build_design_parser() -> _ShellArgParser:
    # SUPPRESS defaults so an unset option is simply absent from the namespace;
    # only options the user actually typed are applied over the current criteria.
    p = _ShellArgParser(
        prog="set", add_help=False, argument_default=argparse.SUPPRESS
    )
    p.add_argument("--length", "-l", type=int)
    p.add_argument("--overhang")
    p.add_argument("--rules")
    p.add_argument("--mode")
    p.add_argument("--reynolds-threshold", type=int)
    p.add_argument("--amarzguioui-threshold", type=int)
    p.add_argument("--seed", action=argparse.BooleanOptionalAction)
    p.add_argument("--seed-tm", type=float)
    p.add_argument("--offtarget", action=argparse.BooleanOptionalAction)
    p.add_argument("--offtarget-db")
    p.add_argument("--offtarget-min-mismatch", type=int)
    p.add_argument("--offtarget-max-mismatch", type=int)
    p.add_argument("--offtarget-self-accession")
    p.add_argument("--bowtie-path")
    p.add_argument("--gc", action=argparse.BooleanOptionalAction)
    p.add_argument("--gc-min", type=float)
    p.add_argument("--gc-max", type=float)
    p.add_argument("--allow-homopolymers", dest="allow_homopolymers", action="store_true")
    p.add_argument("--forbid-homopolymers", dest="allow_homopolymers", action="store_false")
    p.add_argument("--max-run", type=int)
    p.add_argument("--asymmetry", action=argparse.BooleanOptionalAction)
    p.add_argument("--asym-window", type=int)
    p.add_argument("--asym-temp", type=float)
    p.add_argument("--hairpin", action=argparse.BooleanOptionalAction)
    p.add_argument("--hairpin-tm", type=float)
    p.add_argument("--params")
    p.add_argument("--pan-allele", action="store_true")
    p.add_argument("--allele-specific")
    p.add_argument("--target-specific", action="store_true")
    p.add_argument("--hit")
    p.add_argument("--spare")
    p.add_argument("--min-spare-mismatches", type=int)
    p.add_argument("--sort")
    return p


def _build_convert_parser() -> _ShellArgParser:
    p = _ShellArgParser(
        prog="convert", add_help=False, argument_default=argparse.SUPPRESS
    )
    p.add_argument("--loop")
    p.add_argument("--terminator", action=argparse.BooleanOptionalAction)
    p.add_argument("--terminator-length", type=int)
    return p


_DESIGN_PARSER = _build_design_parser()
_CONVERT_PARSER = _build_convert_parser()


def apply_design_flags(params: DesignParams, tokens: list[str]) -> DesignParams:
    """Return a copy of params updated by the design flags in tokens.

    Only options present in tokens change; everything else is inherited. Any
    targeting flag replaces the whole targeting intent (the intents are mutually
    exclusive, and coherence is checked when the design runs).
    """
    provided = vars(_DESIGN_PARSER.parse_args(tokens))
    updates: dict = {}

    def take(key: str, field: str, transform=lambda v: v) -> None:
        if key in provided:
            updates[field] = transform(provided[key])

    take("length", "length")
    take("overhang", "overhang", parse_overhang)
    if "rules" in provided:
        use_ut, use_rey, use_amz = parse_rules(provided["rules"])
        updates.update(
            use_ui_tei=use_ut, use_reynolds=use_rey, use_amarzguioui=use_amz
        )
    take("mode", "rule_mode", parse_mode)
    take("reynolds_threshold", "reynolds_threshold")
    take("amarzguioui_threshold", "amarzguioui_threshold")
    take("seed", "use_seed_filter")
    take("seed_tm", "seed_tm_threshold")
    take("offtarget", "use_offtarget_filter")
    take("offtarget_db", "offtarget_db", Path)
    take("offtarget_min_mismatch", "offtarget_min_mismatch")
    take("offtarget_max_mismatch", "offtarget_max_mismatch")
    take("offtarget_self_accession", "offtarget_self_accession")
    take("bowtie_path", "bowtie_path")
    take("gc", "use_gc_filter")
    take("gc_min", "gc_min", lambda v: v / 100.0)
    take("gc_max", "gc_max", lambda v: v / 100.0)
    take("allow_homopolymers", "use_homopolymer_filter", lambda v: not v)
    take("max_run", "max_homopolymer_run")
    take("asymmetry", "use_asymmetry_filter")
    take("asym_window", "asymmetry_window")
    take("asym_temp", "asymmetry_temp_c")
    take("hairpin", "use_hairpin_filter")
    take("hairpin_tm", "hairpin_tm_max")
    take("params", "thermo_params", parse_params)
    take("min_spare_mismatches", "min_spare_mismatches")
    take("sort", "sort")

    if {
        "pan_allele", "allele_specific", "target_specific", "hit", "spare"
    } & provided.keys():
        updates["targeting"] = AlleleTargetingSpec(
            pan_allele=provided.get("pan_allele", False),
            allele_specific=provided.get("allele_specific"),
            target_specific=provided.get("target_specific", False),
            hit=tuple(split_ids(provided["hit"])) if "hit" in provided else None,
            spare=tuple(split_ids(provided["spare"])) if "spare" in provided else None,
        )
    return replace(params, **updates)


def apply_convert_flags(tokens: list[str]) -> ShrnaOptions:
    """Build ShrnaOptions from the convert flags in tokens (unset ones default)."""
    provided = vars(_CONVERT_PARSER.parse_args(tokens))
    loop = get_loop(provided["loop"]) if "loop" in provided else LOOPS[DEFAULT_LOOP]
    return ShrnaOptions(
        loop=loop,
        include_terminator=provided.get("terminator", True),
        terminator_length=provided.get("terminator_length", DEFAULT_TERMINATOR_LENGTH),
    )


# --------------------------------------------------------------------------- #
# the shell
# --------------------------------------------------------------------------- #

_HELP = """\
Commands:
  load <path>           load the target mRNA/CDS (FASTA or PDB)
  panel <path>          load a related-sequence panel (turns on allele targeting)
  unpanel               drop the panel
  set [flags]           set persistent design criteria (same flags as 'hush design');
                        with no flags, print the current criteria
  reset                 reset the criteria to the best-practice defaults
  limit [n]             set how many rows design/convert print (0 = all)
  show                  show the loaded inputs, criteria, and last results
  inspect               preview the design input (roles, readiness)
  design [flags]        run design with the current criteria (flags apply for this run only;
                        add --explain for the full per-allele discrimination breakdown)
  assess [region]       show every spared allele's discrimination for the last design
                        (all selected, or just <region> e.g. 6-24)
  convert [flags]       convert the last design's selected candidates to shRNA
                        (flags: --loop, --terminator/--no-terminator, --terminator-length)
  save <dir>            write design.tsv (and shrna.tsv, if converted) into <dir>
  help                  show this help
  quit / exit           leave the shell
"""


class HushShell:
    """A read-eval-print loop over one HushSession."""

    def __init__(self, session: HushSession | None = None) -> None:
        self.session = session or HushSession()
        self.params = DesignParams()
        self.limit = 20
        # The last convert result, kept so 'save' writes the scaffold the user
        # actually chose. Cleared when a new design supersedes it.
        self.last_convert = None
        self._commands = {
            "load": self._cmd_load,
            "panel": self._cmd_panel,
            "unpanel": self._cmd_unpanel,
            "set": self._cmd_set,
            "reset": self._cmd_reset,
            "limit": self._cmd_limit,
            "show": self._cmd_show,
            "inspect": self._cmd_inspect,
            "design": self._cmd_design,
            "assess": self._cmd_assess,
            "convert": self._cmd_convert,
            "save": self._cmd_save,
            "help": self._cmd_help,
            "?": self._cmd_help,
            "quit": self._cmd_quit,
            "exit": self._cmd_quit,
            "q": self._cmd_quit,
        }

    # -- loop -------------------------------------------------------------- #

    def run(self, prompt: str = "hush> ") -> None:
        typer.echo("HUSH interactive shell. Type 'help' for commands, 'quit' to exit.")
        while True:
            try:
                line = input(prompt)
            except EOFError:
                typer.echo("")
                break
            except KeyboardInterrupt:
                # Ctrl-C abandons the current line, it does not leave the shell.
                typer.echo("")
                continue
            if not self.dispatch(line):
                break

    def dispatch(self, line: str) -> bool:
        """Run one command line. Returns False only for quit/exit."""
        line = line.strip()
        if not line or line.startswith("#"):
            return True
        cmd, _, rest = line.partition(" ")
        handler = self._commands.get(cmd.lower())
        if handler is None:
            typer.secho(f"unknown command {cmd!r}; type 'help'.", fg=typer.colors.RED)
            return True
        try:
            return handler(rest.strip())
        except (HushError, ValueError) as exc:
            typer.secho(f"error: {exc}", fg=typer.colors.RED)
            return True

    # -- commands ---------------------------------------------------------- #

    def _cmd_load(self, rest: str) -> bool:
        if not rest:
            raise ValueError("usage: load <path>")
        target = self.session.load_target(rest)
        note = "" if target.is_nucleotide else "  [!] needs a coding transcript to design"
        typer.echo(
            f"loaded target: {target.id or '<unnamed>'} "
            f"({target.length} nt, {target.moltype.value}){note}"
        )
        return True

    def _cmd_panel(self, rest: str) -> bool:
        if not rest:
            raise ValueError("usage: panel <path>")
        panel = self.session.load_panel(rest)
        typer.echo(f"loaded panel: {panel.name} ({len(panel)} sequences)")
        return True

    def _cmd_unpanel(self, rest: str) -> bool:
        self.session.clear_panel()
        typer.echo("panel cleared.")
        return True

    def _cmd_set(self, rest: str) -> bool:
        if not rest:
            self._print_criteria()
            return True
        self.params = apply_design_flags(self.params, rest.split())
        return True

    def _cmd_reset(self, rest: str) -> bool:
        self.params = DesignParams()
        typer.echo("criteria reset to the best-practice defaults.")
        return True

    def _cmd_limit(self, rest: str) -> bool:
        if not rest:
            typer.echo(f"limit: {self.limit}")
            return True
        try:
            self.limit = int(rest)
        except ValueError:
            raise ValueError(f"limit takes an integer, got {rest!r}.")
        return True

    def _cmd_show(self, rest: str) -> bool:
        target = self.session.target
        if target is None:
            typer.echo("target: (none)")
        else:
            typer.echo(
                f"target: {target.id or '<unnamed>'} "
                f"({target.length} nt, {target.moltype.value})"
            )
        panel = self.session.panel
        typer.echo(
            f"panel:  {panel.name} ({len(panel)} sequences)"
            if self.session.has_panel
            else "panel:  (none)"
        )
        typer.echo("")
        self._print_criteria()
        typer.echo("")
        last = self.session.last_design
        if last is None:
            typer.echo("last design: (none)")
        else:
            typer.echo(
                f"last design: {len(last.selected)} selected of "
                f"{len(last.evaluations)} candidate(s)"
            )
        if self.last_convert is not None:
            typer.echo(f"last convert: {len(self.last_convert.constructs)} construct(s)")
        return True

    def _cmd_inspect(self, rest: str) -> bool:
        _summarize_design_input(self.session.inspect().design_input)
        return True

    def _cmd_design(self, rest: str) -> bool:
        tokens = rest.split() if rest else []
        # --explain is a render flag, not a persistent design knob, so it is
        # stripped before the criteria parser (which does not know it) runs.
        explain = "--explain" in tokens
        tokens = [t for t in tokens if t != "--explain"]
        params = apply_design_flags(self.params, tokens) if tokens else self.params
        result = self.session.design(params)
        # A previous convert belonged to the old design; drop it so 'save' never
        # pairs a fresh design with stale shRNA.
        self.last_convert = None
        render_design_result(result, self.limit, explain=explain)
        return True

    def _cmd_assess(self, rest: str) -> bool:
        last = self.session.last_design
        if last is None:
            raise ValueError("no design yet; run 'design' first.")
        if (
            last.window_by_pos is None
            or last.targeting is None
            or not last.targeting.discriminating
        ):
            raise ValueError(
                "the last design spares nothing to assess; use an allele-specific, "
                "target-specific, or hit/spare request."
            )
        selected = last.selected
        limit = self.limit
        if rest:
            selected = [e for e in selected if e.candidate.region == rest]
            if not selected:
                raise ValueError(
                    f"no selected candidate with region {rest!r}; regions are the "
                    "left column of the design table (e.g. 6-24)."
                )
            limit = 0  # a named region is a single candidate; never truncate it
        render_discrimination_detail(selected, last.window_by_pos, limit)
        return True

    def _cmd_convert(self, rest: str) -> bool:
        options = apply_convert_flags(rest.split()) if rest else ShrnaOptions()
        result = self.session.convert(options)
        self.last_convert = result
        render_convert_result(result, self.limit)
        return True

    def _cmd_save(self, rest: str) -> bool:
        if not rest:
            raise ValueError("usage: save <dir>")
        out_dir = Path(rest)
        written = []
        last_design = self.session.last_design
        if last_design is not None:
            out_dir.mkdir(parents=True, exist_ok=True)
            path = out_dir / "design.tsv"
            _write_design_tsv(
                path,
                last_design.selected,
                last_design.criteria,
                last_design.window_by_pos,
            )
            written.append(path)
        # Only write shRNA when a convert has actually been run, so the saved
        # scaffold is the one the user chose (not a silent default).
        if self.last_convert is not None:
            path = out_dir / "shrna.tsv"
            _write_convert_tsv(path, self.last_convert.constructs)
            written.append(path)
        if not written:
            raise ValueError("nothing to save yet; run design first.")
        for path in written:
            typer.echo(f"wrote {path}")
        return True

    def _cmd_help(self, rest: str) -> bool:
        typer.echo(_HELP)
        return True

    def _cmd_quit(self, rest: str) -> bool:
        return False

    # -- helpers ----------------------------------------------------------- #

    def _print_criteria(self) -> None:
        """A compact readout of the current persistent criteria."""
        p = self.params
        rules = [
            name
            for name, on in (
                ("ui-tei", p.use_ui_tei),
                ("reynolds", p.use_reynolds),
                ("amarzguioui", p.use_amarzguioui),
            )
            if on
        ]
        typer.echo(typer.style("Criteria", bold=True))
        typer.echo(f"  length/overhang: {p.length} nt, {p.overhang.value}")
        typer.echo(
            f"  rules:           {', '.join(rules) or 'none'} "
            f"(mode {p.rule_mode.value.upper()})"
        )
        typer.echo(
            f"  seed Tm:         {'< %.1f C' % p.seed_tm_threshold if p.use_seed_filter else 'off'}"
        )
        off = "on (needs --offtarget-db)" if p.use_offtarget_filter else "off"
        if p.use_offtarget_filter and p.offtarget_db is not None:
            off = f"on ({p.offtarget_db})"
        typer.echo(f"  off-target:      {off}")
        typer.echo(
            f"  GC window:       {'%.0f-%.0f%%' % (p.gc_min * 100, p.gc_max * 100) if p.use_gc_filter else 'off'}"
        )
        typer.echo(
            f"  hairpin:         {'Tm <= %.1f C' % p.hairpin_tm_max if p.use_hairpin_filter else 'off'}"
        )
        typer.echo(
            f"  homopolymer:     {'reject runs > %d nt' % p.max_homopolymer_run if p.use_homopolymer_filter else 'off'}"
        )
        typer.echo(f"  asymmetry:       {'on' if p.use_asymmetry_filter else 'off'}")
        typer.echo(f"  targeting:       {_format_targeting(p.targeting)}")
        typer.echo(f"  sort:            {p.sort}")


def _format_targeting(spec: AlleleTargetingSpec | None) -> str:
    """One-line description of a targeting intent for the criteria readout."""
    if spec is None or (
        not spec.pan_allele
        and spec.allele_specific is None
        and not spec.target_specific
        and spec.hit is None
    ):
        return "pan-allele (panel only)"
    if spec.allele_specific is not None:
        return f"allele-specific: {spec.allele_specific}"
    if spec.target_specific:
        return "target-specific (spare the whole panel)"
    if spec.hit is not None:
        spare = f", spare {list(spec.spare)}" if spec.spare else ""
        return f"hit {list(spec.hit)}{spare}"
    return "pan-allele"
