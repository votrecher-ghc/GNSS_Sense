"""Regenerate the environment-metrics panel on a square paper canvas."""

from pathlib import Path
import argparse

from reportlab.pdfgen import canvas


SCENARIOS = ["Open\nfield", "Near\nbuilding", "Near\ntrees"]
METRICS = {
    "Accuracy": [95.7176, 82.1181, 82.6100],
    "BAC": [97.5063, 82.5126, 81.6761],
    "F1-score": [79.5014, 43.6131, 43.5681],
    "EER": [2.7936, 17.6452, 18.8289],
}
COLORS = [(0.188, 0.400, 0.671), (0.259, 0.588, 0.471),
          (0.820, 0.541, 0.212), (0.620, 0.341, 0.341)]


def render(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    page = 800.0
    # Match the ROC panel's large square plotting footprint rather than only
    # matching its outer PDF page.
    left, bottom, size = 100.0, 100.0, 650.0
    top = bottom + size
    pdf = canvas.Canvas(str(output), pagesize=(page, page), pageCompression=1)
    pdf.setTitle("Authentication metrics across environments")
    pdf.setFillColorRGB(1, 1, 1)
    pdf.rect(0, 0, page, page, stroke=0, fill=1)

    # Grid, ticks, and square plotting area.
    pdf.setFont("Helvetica", 31)
    for tick in range(0, 101, 10):
        y = bottom + size * tick / 105.0
        pdf.setStrokeColorRGB(0.80, 0.80, 0.80)
        pdf.setLineWidth(0.7)
        pdf.line(left, y, left + size, y)
        label = str(tick)
        pdf.setFillColorRGB(0.15, 0.15, 0.15)
        pdf.drawRightString(left - 12, y - 6, label)

    centers = [left + size * (i + 0.5) / 3 for i in range(3)]
    pdf.setStrokeColorRGB(0.80, 0.80, 0.80)
    for x in centers:
        pdf.line(x, bottom, x, top)

    # Grouped bars; preserve the values and palette of the existing panel.
    bar_width, gap = 28.0, 5.0
    metric_items = list(METRICS.items())
    group_width = len(metric_items) * bar_width + (len(metric_items) - 1) * gap
    for group, center in enumerate(centers):
        x0 = center - group_width / 2
        for metric, ((_, values), color) in enumerate(zip(metric_items, COLORS)):
            height = size * values[group] / 105.0
            pdf.setFillColorRGB(*color)
            pdf.rect(x0 + metric * (bar_width + gap), bottom, bar_width, height,
                     stroke=0, fill=1)

    pdf.setStrokeColorRGB(0.15, 0.15, 0.15)
    pdf.setLineWidth(1.1)
    pdf.rect(left, bottom, size, size, stroke=1, fill=0)

    pdf.setFillColorRGB(0.15, 0.15, 0.15)
    pdf.setFont("Helvetica", 31)
    for center, label in zip(centers, SCENARIOS):
        first, second = label.split("\n")
        pdf.drawCentredString(center, bottom - 36, first)
        pdf.drawCentredString(center, bottom - 60, second)

    pdf.saveState()
    pdf.translate(43, bottom + size / 2)
    pdf.rotate(90)
    pdf.setFont("Helvetica", 50)
    pdf.drawCentredString(0, 0, "Score (%)")
    pdf.restoreState()

    # Compact legend inside the plot, matching the ROC panel's placement.
    legend_width, legend_height = 180.0, 116.0
    legend_x = left + size - legend_width - 13
    legend_y = top - legend_height - 13
    pdf.setFillColorRGB(1, 1, 1)
    pdf.setStrokeColorRGB(0.15, 0.15, 0.15)
    pdf.rect(legend_x, legend_y, legend_width, legend_height, stroke=1, fill=1)
    pdf.setFont("Helvetica", 27)
    for index, ((label, _), color) in enumerate(zip(metric_items, COLORS)):
        y = legend_y + legend_height - 25 - 25 * index
        pdf.setFillColorRGB(*color)
        pdf.rect(legend_x + 11, y - 8, 20, 12, stroke=0, fill=1)
        pdf.setFillColorRGB(0.10, 0.10, 0.10)
        pdf.drawString(legend_x + 39, y - 6, label)

    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    render(parser.parse_args().output)
