#!/usr/bin/env python3
"""Diff a HUSH design TSV against a siDirect result export.

Aligns the two candidate sets by guide sequence and reports which candidates
both tools picked, and which HUSH picked that siDirect did not (with HUSH's
per-rule pass codes). siDirect files are read as raw ODF (a zipped XML bundle)
so no odfpy/pandas dependency is needed.

Usage:
    python tools/compare_sidirect.py --sidirect results.ods --hush hush.tsv

Note: HUSH-only candidates are expected until the transcriptome off-target
filter lands (see the WIP note printed with the diff); siDirect drops candidates
with near-perfect off-targets, which HUSH cannot yet detect.
"""

from __future__ import annotations

import argparse
import csv
import sys
import xml.etree.ElementTree as ET
import zipfile

_NS = {
    "table": "urn:oasis:names:tc:opendocument:xmlns:table:1.0",
    "text": "urn:oasis:names:tc:opendocument:xmlns:text:1.0",
}


def _q(tag: str) -> str:
    prefix, local = tag.split(":")
    return f"{{{_NS[prefix]}}}{local}"


def _read_ods_rows(path: str) -> list[list[str]]:
    """Return the sheet's non-empty rows as lists of cell strings."""
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("content.xml"))
    rows: list[list[str]] = []
    for row in root.iter(_q("table:table-row")):
        cells: list[str] = []
        for cell in row.findall(_q("table:table-cell")):
            # A cell can stand in for a run of identical cells; cap the expansion
            # so trailing-empty runs do not blow up the row.
            repeat = min(int(cell.get(_q("table:number-columns-repeated"), "1")), 50)
            text = "".join(p.text or "" for p in cell.iter(_q("text:p")))
            cells.extend([text] * repeat)
        while cells and cells[-1] == "":
            cells.pop()
        if any(c.strip() for c in cells):
            rows.append(cells)
    return rows


def _col(header: list[str], *needles: str) -> int:
    """Index of the first header cell containing every needle (case-insensitive)."""
    for i, cell in enumerate(header):
        low = cell.lower()
        if all(n in low for n in needles):
            return i
    raise KeyError(f"no siDirect column matching {needles!r} in {header!r}")


def load_sidirect(path: str) -> dict[str, dict[str, str]]:
    """Map guide sequence -> {position, functional, min_mismatch} from a siDirect ODS."""
    rows = _read_ods_rows(path)
    header_idx = next(
        i for i, r in enumerate(rows)
        if any("functional" in c.lower() for c in r)
    )
    header = rows[header_idx]
    i_pos = _col(header, "target position")
    i_guide = _col(header, "oligo", "guide")
    i_func = _col(header, "functional")
    try:
        i_mm = _col(header, "mismatch")
    except KeyError:
        i_mm = None
    out: dict[str, dict[str, str]] = {}
    for row in rows[header_idx + 1:]:
        if len(row) <= i_guide:
            continue
        guide = row[i_guide].strip().upper().replace("T", "U")
        if not guide:
            continue
        out[guide] = {
            "position": row[i_pos] if len(row) > i_pos else "",
            "functional": row[i_func] if len(row) > i_func else "",
            "min_mismatch": row[i_mm] if i_mm is not None and len(row) > i_mm else "",
        }
    return out


def load_hush(path: str) -> dict[str, dict[str, str]]:
    """Map guide sequence -> {region, rules} from a HUSH design TSV."""
    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="ascii", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            guide = row["guide_5to3"].strip().upper().replace("T", "U")
            out[guide] = {
                "region": row.get("region", ""),
                "rules": row.get("rules_passed", ""),
            }
    return out


def _start(region: str) -> int:
    try:
        return int(region.split("-")[0])
    except (ValueError, IndexError):
        return 0


def diff(sidirect: dict, hush: dict) -> int:
    sid_g, hush_g = set(sidirect), set(hush)
    matched = sorted(sid_g & hush_g, key=lambda g: _start(hush[g]["region"]))
    hush_only = sorted(hush_g - sid_g, key=lambda g: _start(hush[g]["region"]))
    sid_only = sorted(sid_g - hush_g, key=lambda g: _start(sidirect[g]["position"]))

    print(f"siDirect: {len(sidirect)}   HUSH: {len(hush)}\n")

    print(f"== MATCHED (both) : {len(matched)} ==")
    for g in matched:
        h, s = hush[g], sidirect[g]
        mm = f" minMM={s['min_mismatch']}" if s["min_mismatch"] else ""
        print(f"  HUSH {h['region']:<10} rules={h['rules']:<3} | "
              f"siDirect {s['position']:<10} func={s['functional']:<3}{mm}")

    print(f"\n== HUSH-ONLY (siDirect excluded) : {len(hush_only)} ==")
    for g in hush_only:
        print(f"  HUSH {hush[g]['region']:<10} rules={hush[g]['rules']:<3} guide={g}")
    if hush_only:
        print("  NOTE: expected until the transcriptome off-target filter is "
              "implemented (WIP);\n        siDirect drops near-perfect-off-target "
              "candidates that HUSH cannot yet see.")

    print(f"\n== siDirect-ONLY (HUSH missed) : {len(sid_only)} ==")
    for g in sid_only:
        print(f"  siDirect {sidirect[g]['position']:<10} "
              f"func={sidirect[g]['functional']:<3} guide={g}")
    if not sid_only:
        print("  (none - siDirect's picks are a subset of HUSH's)")

    # Non-zero exit only if HUSH missed something siDirect kept (a real regression).
    return 1 if sid_only else 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sidirect", required=True, help="siDirect result export (.ods).")
    ap.add_argument("--hush", required=True, help="HUSH design output (.tsv).")
    args = ap.parse_args(argv)
    return diff(load_sidirect(args.sidirect), load_hush(args.hush))


if __name__ == "__main__":
    sys.exit(main())
