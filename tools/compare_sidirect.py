#!/usr/bin/env python3
"""Diff a HUSH design TSV against a siDirect result export.

Matches candidates by guide sequence and reports matched, HUSH-only and
siDirect-only picks. The .ods file is read directly as zipped XML. Run HUSH
with the off-target check active for a like-for-like diff.

    python tools/compare_sidirect.py --sidirect results.ods --hush hush.tsv
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
    """Non-empty rows of the sheet as lists of cell strings."""
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("content.xml"))
    rows: list[list[str]] = []
    for row in root.iter(_q("table:table-row")):
        cells: list[str] = []
        for cell in row.findall(_q("table:table-cell")):
            # A cell may repeat; cap the expansion so empty runs stay small.
            repeat = min(int(cell.get(_q("table:number-columns-repeated"), "1")), 50)
            text = "".join(p.text or "" for p in cell.iter(_q("text:p")))
            cells.extend([text] * repeat)
        while cells and cells[-1] == "":
            cells.pop()
        if any(c.strip() for c in cells):
            rows.append(cells)
    return rows


def _col(header: list[str], *needles: str) -> int:
    """Index of the first header cell containing all needles, ignoring case."""
    for i, cell in enumerate(header):
        low = cell.lower()
        if all(n in low for n in needles):
            return i
    raise KeyError(f"no siDirect column matching {needles!r} in {header!r}")


def load_sidirect(path: str) -> dict[str, dict[str, str]]:
    """Guide sequence to {position, functional, min_mismatch} from a siDirect ODS."""
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
    """Guide sequence to {region, rules, min_mismatch} from a HUSH design TSV."""
    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="ascii", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            guide = row["guide_5to3"].strip().upper().replace("T", "U")
            out[guide] = {
                "region": row.get("region", ""),
                "rules": row.get("rules_passed", ""),
                "min_mismatch": row.get("offtarget_min_mismatch", "").strip(),
            }
    return out


def _start(region: str) -> int:
    try:
        return int(region.split("-")[0])
    except (ValueError, IndexError):
        return 0


def _first_int(text: str):
    """First integer in a cell ('3' or '>=3'), or None."""
    digits = ""
    for ch in text:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    return int(digits) if digits else None


def _mm_agree(hush_mm: str, sid_mm: str) -> bool:
    """True when both report the same nearest-off-target mismatch count."""
    a, b = _first_int(hush_mm), _first_int(sid_mm)
    return a is not None and a == b


def diff(sidirect: dict, hush: dict) -> int:
    sid_g, hush_g = set(sidirect), set(hush)
    matched = sorted(sid_g & hush_g, key=lambda g: _start(hush[g]["region"]))
    hush_only = sorted(hush_g - sid_g, key=lambda g: _start(hush[g]["region"]))
    sid_only = sorted(sid_g - hush_g, key=lambda g: _start(sidirect[g]["position"]))

    print(f"siDirect: {len(sidirect)}   HUSH: {len(hush)}\n")

    mm_disagree = 0
    print(f"== MATCHED (both) : {len(matched)} ==")
    for g in matched:
        h, s = hush[g], sidirect[g]
        sid_mm, hush_mm = s["min_mismatch"], h["min_mismatch"]
        if sid_mm and hush_mm:
            agree = _mm_agree(hush_mm, sid_mm)
            mm = f" minMM: HUSH={hush_mm} siDirect={sid_mm} {'ok' if agree else 'DIFF'}"
            mm_disagree += 0 if agree else 1
        else:
            mm = f" minMM={sid_mm}" if sid_mm else ""
        print(f"  HUSH {h['region']:<10} rules={h['rules']:<3} | "
              f"siDirect {s['position']:<10} func={s['functional']:<3}{mm}")
    if mm_disagree:
        print(f"  NOTE: {mm_disagree} matched candidate(s) disagree on min-mismatch; "
              "check the off-target DB/settings against siDirect's.")

    print(f"\n== HUSH-ONLY (siDirect excluded) : {len(hush_only)} ==")
    for g in hush_only:
        print(f"  HUSH {hush[g]['region']:<10} rules={hush[g]['rules']:<3} guide={g}")
    if hush_only:
        any_offtarget = any(hush[g]["min_mismatch"] for g in hush)
        if any_offtarget:
            print("  NOTE: with the off-target filter on, these should be siDirect's "
                  "on-target/rule differences,\n        not missing off-target "
                  "detection; investigate any remainder.")
        else:
            print("  NOTE: expected when HUSH ran unreferenced; build the "
                  "reference with\n        hush reference build, then re-run "
                  "the design to drop near-perfect-off-target candidates.")

    print(f"\n== siDirect-ONLY (HUSH missed) : {len(sid_only)} ==")
    for g in sid_only:
        print(f"  siDirect {sidirect[g]['position']:<10} "
              f"func={sidirect[g]['functional']:<3} guide={g}")
    if not sid_only:
        print("  (none - siDirect's picks are a subset of HUSH's)")

    # Exit 1 only when HUSH missed a siDirect pick.
    return 1 if sid_only else 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sidirect", required=True, help="siDirect result export (.ods).")
    ap.add_argument("--hush", required=True, help="HUSH design output (.tsv).")
    args = ap.parse_args(argv)
    return diff(load_sidirect(args.sidirect), load_hush(args.hush))


if __name__ == "__main__":
    sys.exit(main())
