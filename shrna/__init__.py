"""HUSH shRNA conversion module (pipeline stage 4).

Converts a validated siRNA candidate into a cloneable stem-loop shRNA. The
19-bp duplex core is used unchanged as the stem, so the guide preserved through
the design stage is exactly the guide the hairpin expresses. Constructs are
emitted in the sense-loop-antisense orientation as the bare hairpin insert
(top/bottom DNA oligos plus a Pol III poly-T terminator), vector-agnostic:
restriction sticky-end tails are added by the wet-lab for their chosen vector.

Each construct is scored rather than template-stamped - the hairpin is folded to
confirm it forms the intended stem-loop, and QC flags (Pol III 5'-A/G initiation,
loop terminator-safety, homopolymer runs) are attached.

Submodules:

* loops   - the established loop registry (Jensen / pSUPER / TRC), each with its
            sequence and citation.
* convert - the options, construct model, and the conversion + scoring logic.

Wired into the CLI as ``hush convert``, which reads the ``hush design`` TSV so
the two stages compose without re-running design. Scope and references are in
docs/shrna_design_note.md; the miR-30 / miR-E scaffold is deliberately out of
scope.
"""
