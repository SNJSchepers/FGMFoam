# FGMFoam agent — a code-aware CLI assistant for the FGMFoam solver.
#
# Copyright (C) 2025-2026 Stijn Schepers, Eindhoven University of Technology
#
# This file is part of FGMFoam. FGMFoam is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# FGMFoam is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# FGMFoam. If not, see <https://www.gnu.org/licenses/>.

"""Field snapshots of a campaign: one image per case, with the measurement below.

How a case is shown comes from its plugin (``CaseTemplate.snapshot``), so this
module holds no case-specific code. The keys of that dict:

    field      the field to show ("T")
    unit       its unit, for the colour bar ("K")
    range      (low, high) of the colour scale; values above ``high`` share its
               colour (a fine scale turned small variations in the burnt gas into
               islands)
    mirror     "ymax": the domain is half of a symmetric one; mirror it about
               its upper edge to show the whole
    upright    True: the flow (x) runs upwards in the image
    origin     "body_end": measure x from the downstream end of the solid body
               found as the gap in the bottom row of cells (as the flame-length
               measurement does); otherwise the domain's own x
    window     (before, after) in metres around the origin along x

Cell positions come from ``postProcess -func writeCellCentres -latestTime``,
run in a case that has no ``Cx`` yet: post-processing, never the solver.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any, Optional

# Phone-readable defaults (the images go into reports and posts).
INK, INK_2, MUTED, GRID, AXIS, SURFACE = (
    "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb")
LINE = "#2a78d6"
SOLID = "#b9b8b1"


def _times(case_dir: Path) -> list[Path]:
    out = []
    for d in case_dir.iterdir():
        try:
            float(d.name)
        except ValueError:
            continue
        if d.is_dir() and d.name != "0":
            out.append(d)
    return sorted(out, key=lambda d: float(d.name))


def latest_fields(case_dir: Path, timeout: int = 600) -> Optional[Path]:
    """The latest time directory, with cell centres written if they are not."""
    times = _times(case_dir)
    if not times:
        return None
    latest = times[-1]
    if not (latest / "Cx").is_file():
        subprocess.run(["postProcess", "-func", "writeCellCentres", "-latestTime"],
                       cwd=case_dir, capture_output=True, text=True, timeout=timeout)
    return latest if (latest / "Cx").is_file() else None


def read_field(path: Path) -> list[float]:
    """The internalField of an ASCII OpenFOAM scalar field written per cell."""
    lines = path.read_text().splitlines()
    i = next(k for k, line in enumerate(lines) if line.startswith("internalField"))
    if "uniform" in lines[i] and "nonuniform" not in lines[i]:
        raise ValueError(f"{path} is uniform, not a per-cell field")
    n = int(lines[i + 1])
    return [float(v) for v in lines[i + 3:i + 3 + n]]


def _strip(ax: Any, time_dir: Path, spec: dict[str, Any]) -> Any:
    import matplotlib.tri as mtri
    import numpy as np
    x = np.array(read_field(time_dir / "Cx"))
    y = np.array(read_field(time_dir / "Cy"))
    f = np.array(read_field(time_dir / spec.get("field", "T")))
    if spec.get("origin") == "body_end":
        bottom = np.sort(x[np.isclose(y, y.min())])
        x0 = bottom[np.argmax(np.diff(bottom)) + 1]
    else:
        x0 = x.min()
    if spec.get("mirror") == "ymax":
        top = y.max() + (y.max() - np.sort(np.unique(y))[-2]) / 2
        x, y, f = np.r_[x, x], np.r_[y, 2 * top - y], np.r_[f, f]
        y = y - top
    across, along = (y, x - x0) if spec.get("upright") else (x - x0, y)
    tri = mtri.Triangulation(across * 1e3, along * 1e3)
    tx, ty = tri.x[tri.triangles], tri.y[tri.triangles]
    span = np.max(np.abs(tx - np.roll(tx, 1, axis=1)) + np.abs(ty - np.roll(ty, 1, axis=1)),
                  axis=1)
    cell = np.median(span)
    tri.set_mask(span > 8 * cell)        # no triangles across solid bodies
    low, high = spec.get("range", (float(f.min()), float(f.max())))
    image = ax.tricontourf(tri, np.clip(f, low, high), levels=np.linspace(low, high, 25),
                           cmap="inferno", extend="max")
    ax.set_facecolor(SOLID)
    before, after = spec.get("window", (None, None))
    if before is not None and spec.get("upright"):
        ax.set_ylim(-before * 1e3, after * 1e3)
    elif before is not None:
        ax.set_xlim(-before * 1e3, after * 1e3)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    return image


def render(campaign: Any, results: dict[str, Any], spec: dict[str, Any],
           metric: dict[str, str], variable_label: str, unit: str) -> tuple[Optional[Path], Optional[str]]:
    """The snapshot figure for a campaign: (path, error)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None, "matplotlib is not installed"
    from of_writer import case_path
    rows = sorted(results["rows"], key=lambda r: r[campaign.variable])
    shown = []
    for row in rows:
        latest = latest_fields(case_path(row["case_name"]))
        if latest is not None:
            shown.append((row, latest))
    if not shown:
        return None, ("no case has written fields yet (or postProcess is not on the PATH), "
                      "so there is nothing to show")
    upright = bool(spec.get("upright"))
    fig = plt.figure(figsize=(6.4, 8.0 if upright else 9.0), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    fig.text(0.04, 0.975, metric.get("title", campaign.name), color=INK, fontsize=18,
             fontweight="bold", va="top")
    fig.text(0.04, 0.943, f"{campaign.case_type} · {campaign.name} · "
             f"{spec.get('field', 'T')} at the last written time",
             color=INK_2, fontsize=12, va="top")
    image = None
    n = len(shown)
    for i, (row, latest) in enumerate(shown):
        if upright:
            width = 0.92 / n
            ax = fig.add_axes([0.04 + i * width, 0.47, width * 0.94, 0.43])
        else:
            height = 0.43 / n
            ax = fig.add_axes([0.04, 0.47 + (n - 1 - i) * height, 0.92, height * 0.92])
        image = _strip(ax, latest, spec)
        label = f"{row[campaign.variable]:g} {unit}".strip()
        if upright:
            ax.set_title(label, color=INK, fontsize=12.5, fontweight="bold", y=-0.075, va="top")
        else:
            ax.text(0.985, 0.07, label, transform=ax.transAxes, color=INK, fontsize=12,
                    fontweight="bold", ha="right", va="bottom")
    low, high = spec.get("range", (None, None))
    cax = fig.add_axes([0.16, 0.425, 0.78, 0.012])
    bar = fig.colorbar(image, cax=cax, orientation="horizontal")
    if high is not None:
        import numpy as np
        step = 10 ** np.floor(np.log10(max(high - low, 1e-12)) - 1)   # round values
        ticks = [round(t / step) * step for t in np.linspace(low, high, 4)]
        ticks[-1] = high
        bar.set_ticks(ticks)
        bar.set_ticklabels([f"{t:g}" for t in ticks[:-1]]
                           + [f"≥{high:g} {spec.get('unit', '')}".strip()])
    bar.ax.tick_params(labelsize=11, colors=MUTED, length=0)
    bar.outline.set_visible(False)
    points = [(r[campaign.variable], r[metric["key"]]) for r, _ in shown
              if r.get(metric["key"]) is not None]
    if points:
        ylabel = metric.get("label", metric["key"])
        # Short lengths read better in mm than as 0.005 m.
        if "[m]" in ylabel and max(abs(v) for _, v in points) < 0.1:
            points = [(u, v * 1e3) for u, v in points]
            ylabel = ylabel.replace("[m]", "[mm]")
        ax = fig.add_axes([0.16, 0.11, 0.78, 0.25])
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
        ax.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
        ax.plot(*zip(*points), color=LINE, linewidth=3, marker="o", markersize=10,
                markeredgecolor=SURFACE, markeredgewidth=2)
        ax.tick_params(colors=MUTED, labelsize=12, length=0, pad=5)
        ax.set_xlabel(variable_label, color=INK_2, fontsize=13)
        ax.set_ylabel(ylabel, color=INK_2, fontsize=13)
        ax.set_ylim(bottom=0)
    path = campaign.directory / f"{campaign.name}_snapshots.png"
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return path, None
