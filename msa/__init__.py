"""HUSH allele module (stage 3): align a panel to the target and evaluate windows.

Panel members are projected onto target coordinates, so a window lines up with
the siRNA candidate at the same position. See docs/architecture.md.

Submodules: alignment (result types), align_star (backend), variation
(per-position comparison), windows (hit and spare gates), targeting (request
model), discrimination (ranking), msa_errors.
"""
