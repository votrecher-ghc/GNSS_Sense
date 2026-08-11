"""Standardize all four panels in the outdoor-environment comparison.

The script preserves the existing confusion-matrix values and ROC curve while
placing every panel on the same 800-by-800 pt canvas.  It also applies one
font hierarchy to the regenerated confusion matrices, metrics panel, and ROC
panel.  The ROC curve vertices are extracted from the existing vector PDF, so
this formatting pass does not change the measured data.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from pypdf import PdfReader
from reportlab.pdfgen import canvas
from reportlab.pdfbase.pdfmetrics import stringWidth

from square_environment_metrics import render as render_metrics


PAGE = 800.0
AXIS_FONT = 50
TICK_FONT = 27
ANNOTATION_FONT = 10
COLORBAR_FONT = 23
CURVE_COLORS = {
    "Open field": (0.1490, 0.4118, 0.7412),
    "Near building": (0.7804, 0.4196, 0.1804),
    "Near trees": (0.2784, 0.5608, 0.2510),
}
EER_LABELS = {
    "Open field": "Open field (EER 2.79%)",
    "Near building": "Near building (EER 17.65%)",
    "Near trees": "Near trees (EER 18.83%)",
}
CLASSES = [
    "L-swipe", "R-swipe", "A", "C", "L", "M",
    "N", "V", "X", "Z", "Star", "Rect.",
]


def extract_matrix(pdf_path: Path) -> list[list[float]]:
    text = PdfReader(str(pdf_path)).pages[0].extract_text() or ""
    values = [float(value) for value in re.findall(r"\b\d+\.\d{2}\b", text)[:144]]
    if len(values) != 144:
        raise ValueError(f"Expected 144 matrix values in {pdf_path}, found {len(values)}")
    return [values[index:index + 12] for index in range(0, 144, 12)]


def blend_color(value: float, maximum: float = 0.9) -> tuple[float, float, float]:
    start = (0.98, 0.98, 0.98)
    end = (0.20, 0.36, 0.58)
    ratio = max(0.0, min(1.0, value / maximum))
    return tuple((1.0 - ratio) * a + ratio * b for a, b in zip(start, end))


def draw_confusion(output: Path, matrix: list[list[float]]) -> None:
    pdf = canvas.Canvas(str(output), pagesize=(PAGE, PAGE), pageCompression=1)
    pdf.setTitle("Environment confusion matrix")
    pdf.setFillColorRGB(1, 1, 1)
    pdf.rect(0, 0, PAGE, PAGE, stroke=0, fill=1)

    left, bottom, heatmap_size = 115.0, 150.0, 600.0
    cell = heatmap_size / len(CLASSES)

    for row, values in enumerate(matrix):
        for col, value in enumerate(values):
            x = left + col * cell
            y = bottom + (len(CLASSES) - row - 1) * cell
            pdf.setFillColorRGB(*blend_color(value))
            pdf.rect(x, y, cell, cell, stroke=0, fill=1)
            pdf.setFillColorRGB(1, 1, 1) if value > 0.55 else pdf.setFillColorRGB(0.10, 0.10, 0.10)
            pdf.setFont("Helvetica", ANNOTATION_FONT)
            pdf.drawCentredString(x + cell / 2, y + cell / 2 - 2, f"{value:.2f}")

    pdf.setStrokeColorRGB(0.15, 0.15, 0.15)
    pdf.setLineWidth(1.0)
    pdf.rect(left, bottom, heatmap_size, heatmap_size, stroke=1, fill=0)

    pdf.setFillColorRGB(0.12, 0.12, 0.12)
    pdf.setFont("Helvetica", TICK_FONT)
    for index, label in enumerate(CLASSES):
        center_x = left + (index + 0.5) * cell
        center_y = bottom + heatmap_size - (index + 0.5) * cell
        pdf.drawRightString(left - 10, center_y - 5, label)
        pdf.saveState()
        pdf.translate(center_x + 4, bottom - 11)
        pdf.rotate(-55)
        pdf.drawRightString(0, 0, label)
        pdf.restoreState()

    pdf.setFont("Helvetica", AXIS_FONT)
    pdf.drawCentredString(left + heatmap_size / 2, 25, "Predicted class")
    pdf.saveState()
    pdf.translate(35, bottom + heatmap_size / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, "True class")
    pdf.restoreState()

    colorbar_x, colorbar_w = 735.0, 20.0
    steps = 120
    for index in range(steps):
        value = 0.9 * index / (steps - 1)
        pdf.setFillColorRGB(*blend_color(value))
        pdf.rect(colorbar_x, bottom + heatmap_size * index / steps,
                 colorbar_w, heatmap_size / steps + 0.5, stroke=0, fill=1)
    pdf.setStrokeColorRGB(0.15, 0.15, 0.15)
    pdf.rect(colorbar_x, bottom, colorbar_w, heatmap_size, stroke=1, fill=0)
    pdf.setFillColorRGB(0.12, 0.12, 0.12)
    pdf.setFont("Helvetica", COLORBAR_FONT)
    for tick in (0.0, 0.3, 0.6, 0.9):
        y = bottom + heatmap_size * tick / 0.9
        pdf.drawString(colorbar_x + colorbar_w + 8, y - 4, f"{tick:.1f}")
    pdf.saveState()
    pdf.translate(790, bottom + heatmap_size / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, "Average score")
    pdf.restoreState()

    pdf.showPage()
    pdf.save()


def _matches_color(operands: list[object], target: tuple[float, float, float]) -> bool:
    if len(operands) != 3:
        return False
    return all(abs(float(actual) - expected) < 0.012
               for actual, expected in zip(operands, target))


def extract_roc_curves(pdf_path: Path) -> dict[str, list[tuple[float, float]]]:
    """Recover the three measured curves from either the original or restyled PDF."""
    operations = PdfReader(str(pdf_path)).pages[0].get_contents().operations
    raw_curves: dict[str, list[tuple[float, float]]] = {}

    for label, color in CURVE_COLORS.items():
        candidates: list[list[tuple[float, float]]] = []
        for index, (operands, operator) in enumerate(operations):
            if operator != b"RG" or not _matches_color(operands, color):
                continue
            points: list[tuple[float, float]] = []
            for path_operands, path_operator in operations[index + 1:]:
                if path_operator in (b"m", b"l") and len(path_operands) == 2:
                    points.append((float(path_operands[0]), float(path_operands[1])))
                elif path_operator in (b"S", b"s"):
                    break
            if points:
                candidates.append(points)
        if not candidates:
            raise ValueError(f"Could not recover {label} ROC curve from {pdf_path}")
        raw_curves[label] = max(candidates, key=len)

    # Each measured curve includes (0, 0) and (1, 1), allowing the axes to be
    # recovered without relying on the source figure's page size or margins.
    all_points = [point for points in raw_curves.values() for point in points]
    x_min = min(x for x, _ in all_points)
    x_max = max(x for x, _ in all_points)
    y_min = min(y for _, y in all_points)
    y_max = max(y for _, y in all_points)
    if x_max <= x_min or y_max <= y_min:
        raise ValueError(f"Degenerate ROC coordinates in {pdf_path}")

    # MATLAB's exported PDF uses a downward plot coordinate system; the
    # ReportLab restyled PDF uses the usual upward system.  Detect both so the
    # conversion is idempotent.
    first_curve = next(iter(raw_curves.values()))
    downward = first_curve[0][1] > first_curve[-1][1]
    normalized: dict[str, list[tuple[float, float]]] = {}
    for label, points in raw_curves.items():
        curve: list[tuple[float, float]] = []
        for x, y in points:
            nx = (x - x_min) / (x_max - x_min)
            ny = ((y_max - y) if downward else (y - y_min)) / (y_max - y_min)
            curve.append((max(0.0, min(1.0, nx)), max(0.0, min(1.0, ny))))
        normalized[label] = curve
    return normalized


def draw_roc(output: Path, curves: dict[str, list[tuple[float, float]]]) -> None:
    pdf = canvas.Canvas(str(output), pagesize=(PAGE, PAGE), pageCompression=1)
    pdf.setTitle("ROC curves across environments")
    pdf.setFillColorRGB(1, 1, 1)
    pdf.rect(0, 0, PAGE, PAGE, stroke=0, fill=1)

    left, bottom, size = 100.0, 100.0, 650.0
    top, right = bottom + size, left + size

    pdf.setFont("Helvetica", TICK_FONT)
    for index in range(11):
        value = index / 10
        x = left + value * size
        y = bottom + value * size
        pdf.setStrokeColorRGB(0.84, 0.84, 0.84)
        pdf.setLineWidth(0.7)
        pdf.line(x, bottom, x, top)
        pdf.line(left, y, right, y)
        label = f"{value:g}"
        pdf.setFillColorRGB(0.15, 0.15, 0.15)
        pdf.drawCentredString(x, bottom - 35, label)
        pdf.drawRightString(left - 12, y - 7, label)

    pdf.setStrokeColorRGB(0.70, 0.70, 0.70)
    pdf.setLineWidth(2.0)
    pdf.setDash(9, 7)
    pdf.line(left, bottom, right, top)
    pdf.setDash()

    for label, points in curves.items():
        path = pdf.beginPath()
        first_x, first_y = points[0]
        path.moveTo(left + first_x * size, bottom + first_y * size)
        for fpr, tpr in points[1:]:
            path.lineTo(left + fpr * size, bottom + tpr * size)
        pdf.setStrokeColorRGB(*CURVE_COLORS[label])
        pdf.setLineWidth(3.0)
        pdf.drawPath(path, stroke=1, fill=0)

    pdf.setStrokeColorRGB(0.15, 0.15, 0.15)
    pdf.setLineWidth(1.1)
    pdf.rect(left, bottom, size, size, stroke=1, fill=0)

    pdf.setFillColorRGB(0.15, 0.15, 0.15)
    pdf.setFont("Helvetica", AXIS_FONT)
    pdf.drawCentredString(left + size / 2, 25, "False positive rate")
    pdf.saveState()
    pdf.translate(43, bottom + size / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, "True positive rate")
    pdf.restoreState()

    legend_entries = ["Chance", *CURVE_COLORS.keys()]
    legend_labels = ["Chance", *(EER_LABELS[label] for label in CURVE_COLORS)]
    sample_width, text_gap = 42.0, 10.0
    legend_width = max(stringWidth(label, "Helvetica", COLORBAR_FONT)
                       for label in legend_labels) + sample_width + text_gap + 24
    row_height = 29.0
    legend_height = row_height * len(legend_entries) + 16
    legend_x = right - legend_width - 12
    legend_y = bottom + 12
    pdf.setFillColorRGB(1, 1, 1)
    pdf.setStrokeColorRGB(0.15, 0.15, 0.15)
    pdf.setLineWidth(1.0)
    pdf.rect(legend_x, legend_y, legend_width, legend_height, stroke=1, fill=1)
    pdf.setFont("Helvetica", COLORBAR_FONT)

    for index, (entry, legend_label) in enumerate(zip(legend_entries, legend_labels)):
        y = legend_y + legend_height - 18 - row_height * index
        x0 = legend_x + 11
        if entry == "Chance":
            pdf.setStrokeColorRGB(0.70, 0.70, 0.70)
            pdf.setDash(7, 5)
            pdf.setLineWidth(2.0)
        else:
            pdf.setStrokeColorRGB(*CURVE_COLORS[entry])
            pdf.setDash()
            pdf.setLineWidth(3.0)
        pdf.line(x0, y - 3, x0 + sample_width, y - 3)
        pdf.setDash()
        pdf.setFillColorRGB(0.10, 0.10, 0.10)
        pdf.drawString(x0 + sample_width + text_gap, y - 10, legend_label)

    pdf.showPage()
    pdf.save()


def standardize(directory: Path) -> None:
    building_path = directory / "scenario_confusion_matrix_near_building.pdf"
    trees_path = directory / "scenario_confusion_matrix_near_trees.pdf"
    metrics_path = directory / "scenario_authentication_metrics_bar.pdf"
    roc_path = directory / "scenario_authentication_roc.pdf"

    building_matrix = extract_matrix(building_path)
    trees_matrix = extract_matrix(trees_path)
    draw_confusion(building_path, building_matrix)
    draw_confusion(trees_path, trees_matrix)
    render_metrics(metrics_path)

    curves = extract_roc_curves(roc_path)
    temporary_roc = directory / ".scenario_authentication_roc.standardized.pdf"
    draw_roc(temporary_roc, curves)
    temporary_roc.replace(roc_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    standardize(parser.parse_args().directory)
