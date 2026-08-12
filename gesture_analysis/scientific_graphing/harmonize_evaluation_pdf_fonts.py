"""Harmonize the visible font scale of StarDial evaluation PDFs.

The source figures use different page sizes and are inserted at different
LaTeX widths.  Consequently, equal source font sizes do not look equal in the
compiled paper.  This utility applies conservative per-figure scale factors to
PDF text operators so that the compiled labels visually follow Fig. 11(b),
while avoiding oversized legends and dense confusion-matrix annotations.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.generic import FloatObject


# Exact source font sizes are used to distinguish axes/ticks from legends and
# mathematical annotations.  Values are intentionally not forced to one hard
# size: half-width plots need larger source text than full-width plots.
FONT_MAPS: dict[str, dict[float, float]] = {
    "grid_avg_affected_satellites_vs_height.pdf": {
        66.667: 76.0, 80.0: 76.0,
        83.333: 90.0, 96.0: 90.0,
        70.0: 50.0,
    },
    "height_sensitivity_dual_axis.pdf": {
        66.667: 76.0, 80.0: 76.0,
        83.333: 90.0, 96.0: 90.0,
        70.0: 50.0,
    },
    "figure_satellite_range_difference_distribution.pdf": {
        50.0: 60.0, 70.0: 60.0,
        52.0: 40.0,
    },
    "figure_space_uniqueness_heatmap.pdf": {
        50.0: 60.0, 68.0: 60.0,
    },
    "figure_spatial_auth_pass_rate.pdf": {
        40.0: 48.0, 62.0: 48.0, 52.0: 48.0,
        38.0: 29.333, 30.0: 23.467,
    },
    "figure_temporal_auth_pass_rate.pdf": {
        40.0: 48.0, 62.0: 48.0, 52.0: 48.0,
        38.0: 29.333, 30.0: 23.467,
    },
    "feature_space_tsne.pdf": {
        46.667: 56.0, 63.333: 78.0, 86.667: 104.0,
    },
    "feature_space_pca.pdf": {
        51.667: 58.0, 66.667: 78.0, 90.0: 106.0,
    },
}

# The long vertical label in the two pass-rate panels cannot use the same
# enlarged source size as the shorter x-axis label without crossing the PDF
# page boundary.  Keep only that rotated label at a fit-safe size while the
# ticks and x-axis title retain the larger half-width-figure scale.
ROTATED_AXIS_FONT_SIZES = {
    "figure_spatial_auth_pass_rate.pdf": 40.0,
    "figure_temporal_auth_pass_rate.pdf": 40.0,
}


def _encoded_text_length(operand: object) -> int:
    values = operand if isinstance(operand, list) else [operand]
    return sum(len(str(value)) for value in values
               if not isinstance(value, (int, float)))


def resize_long_rotated_label(contents: object, target_size: float) -> int:
    """Resize the long vertical axis label without shrinking other labels."""
    active_font = None
    active_size = None
    rotated_stack = [False]
    changed = 0

    index = 0
    while index < len(contents.operations):
        operands, operator = contents.operations[index]
        if operator == b"q":
            rotated_stack.append(rotated_stack[-1])
        elif operator == b"Q":
            if len(rotated_stack) > 1:
                rotated_stack.pop()
        elif operator == b"cm" and len(operands) == 6:
            a, b, c, d = (float(value) for value in operands[:4])
            if abs(a) < 0.01 and abs(d) < 0.01 and abs(b) > 0.1 and abs(c) > 0.1:
                rotated_stack[-1] = True
        elif operator == b"Tf" and len(operands) >= 2:
            active_font = operands[0]
            active_size = float(operands[1])
        elif (operator in (b"Tj", b"TJ") and rotated_stack[-1]
              and operands and _encoded_text_length(operands[0]) >= 20
              and active_font is not None
              and (active_size is None or abs(active_size - target_size) >= 0.01)):
            contents.operations.insert(
                index, ([active_font, FloatObject(target_size)], b"Tf"))
            active_size = target_size
            changed += 1
            index += 1
        index += 1
    return changed


def remove_shadowed_font_operators(contents: object) -> int:
    """Remove a font selection immediately overwritten before any text."""
    changed = 0
    index = 0
    while index + 1 < len(contents.operations):
        if (contents.operations[index][1] == b"Tf"
                and contents.operations[index + 1][1] == b"Tf"):
            del contents.operations[index]
            changed += 1
            continue
        index += 1
    return changed


def replace_font_sizes(pdf_path: Path, mapping: dict[float, float],
                       rotated_axis_size: float | None = None) -> int:
    reader = PdfReader(str(pdf_path))
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    changed = 0

    for page in writer.pages:
        contents = page.get_contents()
        if contents is None:
            continue
        rotated_stack = [False]
        for operands, operator in contents.operations:
            if operator == b"q":
                rotated_stack.append(rotated_stack[-1])
            elif operator == b"Q":
                if len(rotated_stack) > 1:
                    rotated_stack.pop()
            elif operator == b"cm" and len(operands) == 6:
                a, b, c, d = (float(value) for value in operands[:4])
                if (abs(a) < 0.01 and abs(d) < 0.01
                        and abs(b) > 0.1 and abs(c) > 0.1):
                    rotated_stack[-1] = True
            if operator != b"Tf" or len(operands) < 2:
                continue
            current = float(operands[1])
            if (rotated_axis_size is not None and rotated_stack[-1]
                    and abs(current - rotated_axis_size) < 0.01):
                continue
            for source, target in mapping.items():
                if abs(current - source) < 0.01:
                    operands[1] = FloatObject(target)
                    changed += 1
                    break
        if rotated_axis_size is not None:
            changed += resize_long_rotated_label(contents, rotated_axis_size)
        changed += remove_shadowed_font_operators(contents)
        page.replace_contents(contents)

    temporary = pdf_path.with_name(f".{pdf_path.name}.font-harmonized.tmp")
    with temporary.open("wb") as stream:
        writer.write(stream)
    temporary.replace(pdf_path)
    return changed


def harmonize(paper_root: Path) -> None:
    locations = {
        "grid_avg_affected_satellites_vs_height.pdf": paper_root / "figures",
        "height_sensitivity_dual_axis.pdf": paper_root / "figures",
        "figure_satellite_range_difference_distribution.pdf": paper_root / "figures" / "dif_time_space",
        "figure_space_uniqueness_heatmap.pdf": paper_root / "figures" / "dif_time_space",
        "figure_spatial_auth_pass_rate.pdf": paper_root / "figures" / "dif_time_space",
        "figure_temporal_auth_pass_rate.pdf": paper_root / "figures" / "dif_time_space",
        "feature_space_tsne.pdf": paper_root / "figures" / "attacks",
        "feature_space_pca.pdf": paper_root / "figures" / "attacks",
    }
    for name, mapping in FONT_MAPS.items():
        path = locations[name] / name
        if not path.exists():
            raise FileNotFoundError(path)
        count = replace_font_sizes(
            path, mapping, ROTATED_AXIS_FONT_SIZES.get(name))
        print(f"{name}: updated {count} font operators")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paper_root", type=Path)
    harmonize(parser.parse_args().paper_root)
