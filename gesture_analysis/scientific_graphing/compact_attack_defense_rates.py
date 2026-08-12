"""Render Fig. 16 on a compact, wide vector canvas.

The values below reproduce the current attack-defense PDF exactly.  This
layout pass changes only the aspect ratio, margins, and legend placement so
that the figure consumes less vertical space in a single paper column.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas


GROUPS = ["Legitimate", "Replay", "SDR Spoof", "Ghost/Injection"]
SUCCESS = [95.8333, 8.3333, 8.3333, 4.1667]
ERROR = [4.1667, 8.3333, 12.5000, 8.3333]
SUCCESS_COLOR = (0.1882, 0.4000, 0.6706)
ERROR_COLOR = (0.8196, 0.5412, 0.2118)

PAGE_W = 1056.0
PAGE_H = 575.0
LEFT = 116.0
BOTTOM = 146.0
PLOT_W = 890.0
PLOT_H = 345.0


def _draw_legend(pdf: canvas.Canvas) -> None:
    font_size = 30
    swatch_w = 31
    swatch_h = 18
    gap = 9
    item_gap = 28
    labels = ["Success metric (TAR / ASR)", "Error metric (FRR / FAR)"]
    colors = [SUCCESS_COLOR, ERROR_COLOR]
    widths = [swatch_w + gap + stringWidth(label, "Helvetica", font_size)
              for label in labels]
    total_w = sum(widths) + item_gap
    x = LEFT + PLOT_W - total_w
    y = 527.0

    pdf.setFont("Helvetica", font_size)
    for label, color, width in zip(labels, colors, widths):
        pdf.setFillColorRGB(*color)
        pdf.rect(x, y - 4, swatch_w, swatch_h, stroke=0, fill=1)
        pdf.setFillColorRGB(0.10, 0.10, 0.10)
        pdf.drawString(x + swatch_w + gap, y - 5, label)
        x += width + item_gap


def render(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(output), pagesize=(PAGE_W, PAGE_H), pageCompression=1)
    pdf.setTitle("Attack-defense rates")
    pdf.setFillColorRGB(1, 1, 1)
    pdf.rect(0, 0, PAGE_W, PAGE_H, stroke=0, fill=1)

    # Grid and y-axis ticks.
    pdf.setFont("Helvetica", 34)
    for tick in range(0, 101, 10):
        y = BOTTOM + PLOT_H * tick / 105.0
        pdf.setStrokeColorRGB(0.86, 0.86, 0.86)
        pdf.setLineWidth(0.7)
        pdf.line(LEFT, y, LEFT + PLOT_W, y)
        pdf.setFillColorRGB(0.12, 0.12, 0.12)
        pdf.drawRightString(LEFT - 13, y - 10, str(tick))

    centers = [LEFT + PLOT_W * (index + 0.5) / len(GROUPS)
               for index in range(len(GROUPS))]
    bar_w = 52.0
    gap = 10.0
    for center, success, error in zip(centers, SUCCESS, ERROR):
        success_x = center - gap / 2 - bar_w
        error_x = center + gap / 2
        pdf.setFillColorRGB(*SUCCESS_COLOR)
        pdf.rect(success_x, BOTTOM, bar_w, PLOT_H * success / 105.0,
                 stroke=0, fill=1)
        pdf.setFillColorRGB(*ERROR_COLOR)
        pdf.rect(error_x, BOTTOM, bar_w, PLOT_H * error / 105.0,
                 stroke=0, fill=1)

    pdf.setStrokeColorRGB(0.12, 0.12, 0.12)
    pdf.setLineWidth(1.1)
    pdf.rect(LEFT, BOTTOM, PLOT_W, PLOT_H, stroke=1, fill=0)

    # Rotated category labels retain the original 18-degree presentation.
    pdf.setFont("Helvetica", 38)
    pdf.setFillColorRGB(0.12, 0.12, 0.12)
    for center, label in zip(centers, GROUPS):
        pdf.saveState()
        pdf.translate(center, BOTTOM - 48)
        pdf.rotate(18)
        pdf.drawCentredString(0, -10, label)
        pdf.restoreState()

    pdf.setFont("Helvetica", 48)
    pdf.drawCentredString(LEFT + PLOT_W / 2, 22, "Sample class / attack type")
    pdf.saveState()
    pdf.translate(30, BOTTOM + PLOT_H / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, "Rate (%)")
    pdf.restoreState()

    _draw_legend(pdf)
    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    render(parser.parse_args().output)
