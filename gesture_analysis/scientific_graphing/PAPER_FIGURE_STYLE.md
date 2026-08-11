# StarDial paper figure style

Use the following rules when generating or revising paper figures.  The target
is the size of text after a figure has been inserted into the IEEE two-column
paper, not the raw font size used on a MATLAB canvas.

## Typography

- Font family: Arial for all axes, titles, legends, annotations, and colorbars.
- Final printed size: 8 pt for axis labels and panel titles; 7.5 pt for tick
  labels and legends; never below 7 pt.
- Use sentence case.  Units belong in parentheses, for example `East (cm)` and
  `RMSE (cm)`.
- Do not enlarge one figure independently to make it readable.  Set its export
  canvas and font together so it reaches the target printed size.

For the current MATLAB exporter, the 1600 x 1200 standard canvas uses
`paper_font_sizes_local` as its baseline.  The 1920-pixel-wide trajectory
gallery uses the dedicated `gallery_font_sizes_local` values because it contains
more panels.  Exceptional 40--60 pt overrides should be removed when the
affected figures are next regenerated.

## Layout

- Single-panel plots should normally occupy one IEEE column.
- Multi-panel comparisons may occupy both columns, but repeated axis labels
  must be replaced by shared labels.
- Keep one legend per figure.  Prefer a frameless legend in unused outer space.
- Use the same physical axis height for figures that are placed side by side.
- Crop external whitespace during export; do not repair inconsistent canvases
  by scaling figures differently in LaTeX.

## Lines and marks

- Primary lines: 1.5--2.0 pt at the final printed size.
- Secondary grids: thin, light gray, and visually behind the data.
- Use both color and line/marker shape when a distinction is security-critical.
- Start/end markers in trajectory plots should remain visible but subordinate
  to the two trajectory curves.

## LaTeX insertion

- One-column figure: `width=\columnwidth`.
- Two-column figure: use `figure*` and normally `width=0.92--0.96\textwidth`.
- Avoid per-figure `scale`, `height`, or arbitrary widths used only to correct
  text size.  Regenerate the source figure with the shared style instead.
