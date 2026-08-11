"""Reflow the existing 4x3 trajectory gallery into a compact 6x2 figure.

This utility is a fallback for machines without MATLAB.  The canonical figure
generator is ``export_paper_figures_data_driven.m``; it now produces the same
6x2 layout from the original trajectory arrays.  This script crops only the 12
plot regions from an already rendered high-resolution gallery and redraws all
titles, shared axes, tick labels, and the legend at a consistent paper size.

Dependencies: Pillow and ReportLab.
"""

from __future__ import annotations

import argparse
import io
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import Color, black, white
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen.canvas import Canvas


LABELS = (
    "Left swipe",
    "Right swipe",
    "A",
    "C",
    "L",
    "M",
    "N",
    "V",
    "X",
    "Z",
    "Star",
    "Rectangle",
)

GT_COLOR = (0.302, 0.4902, 0.6784)
REC_COLOR = (0.6588, 0.3098, 0.3216)


def _group_consecutive(indices: np.ndarray) -> list[int]:
    if indices.size == 0:
        return []
    groups: list[list[int]] = [[int(indices[0])]]
    for value in indices[1:]:
        value = int(value)
        if value == groups[-1][-1] + 1:
            groups[-1].append(value)
        else:
            groups.append([value])
    return [int(round(sum(group) / len(group))) for group in groups]


def detect_axes(source: Image.Image) -> list[tuple[int, int, int, int]]:
    """Detect the four column and three row axes borders in the old gallery."""

    gray = np.asarray(source.convert("L"))
    dark = gray < 90
    x_lines = _group_consecutive(np.flatnonzero(dark.sum(axis=0) > 0.25 * gray.shape[0]))
    y_lines = _group_consecutive(np.flatnonzero(dark.sum(axis=1) > 0.25 * gray.shape[1]))
    if len(x_lines) < 8 or len(y_lines) < 6:
        raise RuntimeError(
            f"could not detect 4x3 axes grid (x={x_lines}, y={y_lines})"
        )
    x_lines = x_lines[:8]
    y_lines = y_lines[:6]
    boxes: list[tuple[int, int, int, int]] = []
    for row in range(3):
        for column in range(4):
            boxes.append(
                (
                    x_lines[2 * column],
                    y_lines[2 * row],
                    x_lines[2 * column + 1] + 1,
                    y_lines[2 * row + 1] + 1,
                )
            )
    return boxes


def crop_panels(source: Image.Image) -> list[Image.Image]:
    return [source.crop(box).convert("RGB") for box in detect_axes(source)]


def _panel_positions() -> tuple[float, float, float, float, float, float]:
    page_width = 960.0
    page_height = 410.0
    axis_size = 120.0
    column_gap = 28.0
    left = 82.0
    bottom_rows = (250.0, 90.0)
    return page_width, page_height, axis_size, column_gap, left, bottom_rows


def _draw_centered_text(canvas: Canvas, text: str, x: float, y: float, size: float) -> None:
    canvas.setFont("Helvetica", size)
    canvas.setFillColor(black)
    canvas.drawCentredString(x, y, text)


def build_pdf(panels: list[Image.Image], output: Path) -> None:
    page_width, page_height, axis_size, gap, left, row_bottoms = _panel_positions()
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas = Canvas(str(output), pagesize=(page_width, page_height), pageCompression=1)
    canvas.setTitle("Compact StarDial trajectory reconstruction gallery")
    canvas.setFillColor(white)
    canvas.rect(0, 0, page_width, page_height, stroke=0, fill=1)

    tick_font = 14.0
    label_font = 16.0
    title_font = 17.0
    tick_values = ("-20", "0", "20")

    for index, (label, panel) in enumerate(zip(LABELS, panels)):
        row, column = divmod(index, 6)
        x = left + column * (axis_size + gap)
        y = row_bottoms[row]
        buffer = io.BytesIO()
        panel.save(buffer, format="PNG", optimize=True)
        buffer.seek(0)
        canvas.drawImage(
            ImageReader(buffer), x, y, width=axis_size, height=axis_size, mask="auto"
        )
        _draw_centered_text(canvas, label, x + axis_size / 2, y + axis_size + 6, title_font)

        if row == 1:
            for fraction, tick in zip((0.25, 0.5, 0.75), tick_values):
                _draw_centered_text(canvas, tick, x + fraction * axis_size, y - 13, tick_font)
        if column == 0:
            canvas.setFont("Helvetica", tick_font)
            canvas.setFillColor(black)
            for fraction, tick in zip((0.25, 0.5, 0.75), tick_values):
                canvas.drawRightString(x - 7, y + fraction * axis_size - 4, tick)

    _draw_centered_text(canvas, "East (cm)", page_width / 2, 48.0, label_font)
    canvas.saveState()
    canvas.translate(16, page_height / 2 + 5)
    canvas.rotate(90)
    _draw_centered_text(canvas, "North (cm)", 0, 0, label_font)
    canvas.restoreState()

    legend_y = 18.0
    legend_size = 14.0
    legend_items = (
        ("line", GT_COLOR, "Ground truth"),
        ("line", REC_COLOR, "Recovered trajectory"),
        ("circle", (0.22, 0.22, 0.22), "Gesture Start"),
        ("square", (0.22, 0.22, 0.22), "Gesture End"),
    )
    widths = []
    for kind, _, label in legend_items:
        widths.append(24 + pdfmetrics.stringWidth(label, "Helvetica", legend_size) + 18)
    cursor = (page_width - sum(widths)) / 2
    canvas.setFont("Helvetica", legend_size)
    for item_width, (kind, color, label) in zip(widths, legend_items):
        rgb = Color(*color)
        mark_x = cursor + 8
        if kind == "line":
            canvas.setStrokeColor(rgb)
            canvas.setLineWidth(2.2)
            canvas.line(mark_x, legend_y + 4, mark_x + 18, legend_y + 4)
        elif kind == "circle":
            canvas.setFillColor(rgb)
            canvas.circle(mark_x + 9, legend_y + 4, 2.4, stroke=0, fill=1)
        else:
            canvas.setStrokeColor(rgb)
            canvas.setFillColor(white)
            canvas.rect(mark_x + 6, legend_y + 1, 6, 6, stroke=1, fill=1)
        canvas.setFillColor(black)
        canvas.drawString(mark_x + 23, legend_y, label)
        cursor += item_width

    canvas.showPage()
    canvas.save()


def _find_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = (
        Path("C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/Arial.ttf"),
    )
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def build_preview(panels: list[Image.Image], output: Path) -> None:
    """Create a 1920x700 preview with the same geometry as the PDF."""

    page_width, page_height, axis_size, gap, left, row_bottoms = _panel_positions()
    scale = 2.0
    canvas = Image.new("RGB", (int(page_width * scale), int(page_height * scale)), "white")
    draw = ImageDraw.Draw(canvas)
    tick_font = _find_font(28)
    label_font = _find_font(32)
    title_font = _find_font(34)

    def centered(text: str, x: float, y: float, font: ImageFont.ImageFont) -> None:
        box = draw.textbbox((0, 0), text, font=font)
        draw.text((x - (box[2] - box[0]) / 2, y), text, fill="black", font=font)

    for index, (label, panel) in enumerate(zip(LABELS, panels)):
        row, column = divmod(index, 6)
        x = (left + column * (axis_size + gap)) * scale
        y_pdf = row_bottoms[row]
        y = (page_height - y_pdf - axis_size) * scale
        resized = panel.resize((int(axis_size * scale), int(axis_size * scale)), Image.Resampling.LANCZOS)
        canvas.paste(resized, (int(x), int(y)))
        centered(label, x + axis_size * scale / 2, y - 39, title_font)
        if row == 1:
            for fraction, tick in zip((0.25, 0.5, 0.75), ("-20", "0", "20")):
                centered(tick, x + fraction * axis_size * scale, y + axis_size * scale + 5, tick_font)
        if column == 0:
            for fraction, tick in zip((0.25, 0.5, 0.75), ("20", "0", "-20")):
                box = draw.textbbox((0, 0), tick, font=tick_font)
                draw.text(
                    (x - 12 - (box[2] - box[0]), y + fraction * axis_size * scale - 14),
                    tick,
                    fill="black",
                    font=tick_font,
                )

    centered("East (cm)", page_width * scale / 2, page_height * scale - 124, label_font)
    y_label_layer = Image.new("RGBA", (280, 60), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_label_layer)
    y_draw.text((0, 5), "North (cm)", fill="black", font=label_font)
    rotated = y_label_layer.rotate(90, expand=True)
    canvas.paste(rotated, (5, int(page_height * scale / 2 - rotated.height / 2)), rotated)

    legend_y = int(page_height * scale - 54)
    legend = (
        (GT_COLOR, "Ground truth"),
        (REC_COLOR, "Recovered trajectory"),
        ((0.22, 0.22, 0.22), "Gesture Start"),
        ((0.22, 0.22, 0.22), "Gesture End"),
    )
    widths = [45 + draw.textlength(label, font=tick_font) + 32 for _, label in legend]
    cursor = (canvas.width - sum(widths)) / 2
    for item_index, ((color, label), width) in enumerate(zip(legend, widths)):
        rgb = tuple(round(255 * value) for value in color)
        if item_index < 2:
            draw.line((cursor, legend_y + 14, cursor + 34, legend_y + 14), fill=rgb, width=4)
        elif item_index == 2:
            draw.ellipse((cursor + 12, legend_y + 8, cursor + 24, legend_y + 20), fill=rgb)
        else:
            draw.rectangle((cursor + 12, legend_y + 8, cursor + 24, legend_y + 20), outline=rgb, width=2)
        draw.text((cursor + 42, legend_y), label, fill="black", font=tick_font)
        cursor += width

    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, format="PNG", optimize=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_png", type=Path)
    parser.add_argument("output_pdf", type=Path)
    parser.add_argument("--preview", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with Image.open(args.source_png) as source:
        panels = crop_panels(source.convert("RGB"))
    build_pdf(panels, args.output_pdf)
    if args.preview:
        build_preview(panels, args.preview)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
