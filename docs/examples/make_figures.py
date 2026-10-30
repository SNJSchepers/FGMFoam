"""Regenerate the example figures in docs/examples/ from agent campaign output.

    agent/.venv/bin/python docs/examples/make_figures.py collect   # CSVs from campaigns
    agent/.venv/bin/python docs/examples/make_figures.py plot      # PNGs from the CSVs

`collect` (and `flames`) need the campaigns and cases (git-ignored); `plot` needs
only the CSVs committed here, so the figures can be restyled without rerunning.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "agent"))

# Validated two-slot categorical palette (blue, orange) and chart chrome.
SERIES = ("#2a78d6", "#eb6834")
SURFACE, INK, INK_2, MUTED, GRID, AXIS = (
    "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7")

SPEED_CSV = HERE / "laminar_flame_speed.csv"
LENGTH_CSV = HERE / "flame_length_vs_inlet_velocity.csv"
CHEM1D_CSV = HERE / "chem1d_flame_speed.csv"
FLAMES_PNG = HERE / "slot_burner_flames.png"
FLAME_CASES = "agent/campaigns/LvsUlong"


def collect() -> None:
    import cantera_tool
    # SLvsT20: 20 ms per flamelet, so the start-up transient has died out.
    rows = list(csv.DictReader(open(REPO / "agent/campaigns/SLvsT20/results.csv")))
    cache = json.loads(cantera_tool._FLAME_CACHE.read_text())
    chem1d = {float(r["T_unburnt_K"]): float(r["S_L_Chem1D_m_s"])
              for r in csv.DictReader(l for l in open(CHEM1D_CSV) if not l.startswith("#"))}
    with SPEED_CSV.open("w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["T_unburnt_K", "S_L_FGMFoam_m_s", "S_L_Chem1D_m_s", "S_L_Cantera_UL_m_s",
                      "vs_Chem1D_pct", "vs_Cantera_pct"])
        for r in sorted(rows, key=lambda r: float(r["T_inlet"])):
            T, s = float(r["T_inlet"]), float(r["laminar_flame_speed"])
            c = cache[cantera_tool._cache_key(0.7, T, 101325.0)]["flame_speed"]
            k = chem1d.get(T)
            out.writerow([f"{T:g}", f"{s:.4f}", f"{k:.4f}" if k else "", f"{c:.4f}",
                          f"{100 * (s - k) / k:+.1f}" if k else "", f"{100 * (s - c) / c:+.1f}"])
    print("wrote", SPEED_CSV.name)
    # LvsUlong: 35 mm outlet, so the tip stays inside the domain up to 0.6 m/s.
    source = REPO / "agent/campaigns/LvsUlong/results.csv"
    if source.is_file():
        rows = list(csv.DictReader(open(source)))
        with LENGTH_CSV.open("w", newline="") as fh:
            out = csv.writer(fh)
            out.writerow(["inlet_velocity_m_s", "flame_length_mm", "outcome",
                          "flame_length_lower_bound_mm"])
            for r in sorted(rows, key=lambda r: float(r["inlet_velocity"])):
                L = r.get("flame_length")
                # A tip beyond the outlet still bounds the length from below.
                bound = ""
                if r["status"] == "longer_than_domain":
                    bound = _lower_bound(REPO / "testcases" / r["case_name"])
                out.writerow([r["inlet_velocity"],
                              f"{float(L) * 1e3:.3f}" if L else "", r["status"], bound])
        print("wrote", LENGTH_CSV.name)


def _lower_bound(case_dir: Path) -> str:
    """Outlet minus trailing edge [mm] for a flame longer than the domain."""
    import postrun
    result = postrun.analyse(case_dir).get("flame", {})
    bound = result.get("flame_length_lower_bound")
    return f"{bound * 1e3:.3f}" if bound else ""


def _axes(title: str, subtitle: str, xlabel: str, ylabel: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Read on a phone in a social feed: near-square, large text, thick lines.
    fig, ax = plt.subplots(figsize=(6.4, 5.6), dpi=200)
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=13, length=0, pad=6)
    ax.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=INK_2, fontsize=14, labelpad=8)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=14, labelpad=8)
    fig.text(0.04, 0.965, title, color=INK, fontsize=18, fontweight="bold", va="top")
    fig.text(0.04, 0.905, subtitle, color=INK_2, fontsize=12.5, va="top")
    fig.subplots_adjust(left=0.14, right=0.95, top=0.82, bottom=0.13)
    return fig, ax, plt


def plot_speed() -> None:
    rows = list(csv.DictReader(open(SPEED_CSV)))
    T = [float(r["T_unburnt_K"]) for r in rows]
    fgm = [float(r["S_L_FGMFoam_m_s"]) for r in rows]
    cantera = [float(r["S_L_Cantera_UL_m_s"]) for r in rows]
    chem = [(t, float(r["S_L_Chem1D_m_s"])) for t, r in zip(T, rows) if r["S_L_Chem1D_m_s"]]
    fig, ax, plt = _axes(
        "Laminar flame speed",
        "Methane–air, φ = 0.7, 1 atm · 9-case campaign by the agent",
        "Unburnt temperature  [K]", "Flame speed  [m/s]")
    # Cantera stays in the CSV; the figure compares FGMFoam with its own table's source.
    ax.plot(*zip(*chem), color=SERIES[1], linewidth=3, label="Chem1D (table source)")
    ax.plot(T, fgm, color=SERIES[0], linewidth=0, marker="o", markersize=10,
            markeredgecolor=SURFACE, markeredgewidth=2, label="FGMFoam")
    # The three lines end too close together for labels at their ends.
    ax.legend(loc="upper left", frameon=False, fontsize=13, labelcolor=INK_2,
              handlelength=1.6, borderaxespad=0.2)
    ax.annotate("FGMFoam within 0.3 %\nof Chem1D", (0.97, 0.05), xycoords="axes fraction",
                ha="right", va="bottom", color=INK_2, fontsize=12.5)
    ax.set_ylim(bottom=0)
    out = HERE / "laminar_flame_speed.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig); print("wrote", out.name)


def plot_length() -> None:
    if not LENGTH_CSV.is_file():
        print("skip flame length: no", LENGTH_CSV.name, "yet"); return
    rows = list(csv.DictReader(open(LENGTH_CSV)))
    measured = [(float(r["inlet_velocity_m_s"]), float(r["flame_length_mm"]))
                for r in rows if r["flame_length_mm"]]
    bounded = [(float(r["inlet_velocity_m_s"]), float(r["flame_length_lower_bound_mm"]))
               for r in rows if not r["flame_length_mm"]
               and r.get("flame_length_lower_bound_mm")]
    other = [(float(r["inlet_velocity_m_s"]), r["outcome"]) for r in rows
             if not r["flame_length_mm"] and not r.get("flame_length_lower_bound_mm")]
    fig, ax, plt = _axes(
        "Flame length",
        "2D laminar slot burner, methane–air, φ = 0.7",
        "Inlet velocity  [m/s]", "Flame length  [mm]")
    if measured:
        ax.plot(*zip(*measured), color=SERIES[0], linewidth=3, marker="o", markersize=10,
                markeredgecolor=SURFACE, markeredgewidth=2,
                label="tip on the slot centreline (900 K)")
    if bounded:
        # Open markers with an upward arrow: the tip left the domain, so the
        # outlet only gives a lower bound.
        ax.plot(*zip(*bounded), linestyle="none", marker="o", markersize=6.5,
                markerfacecolor=SURFACE, markeredgecolor=SERIES[0], markeredgewidth=1.5,
                label="tip at or beyond the outlet: lower bound")
        for u, L in bounded:
            ax.annotate("", xy=(u, L + 4.0), xytext=(u, L + 0.6),
                        arrowprops=dict(arrowstyle="-|>", color=SERIES[0], linewidth=1.4))
        ax.annotate("longer than the domain", (bounded[0][0], bounded[0][1] + 2.3),
                    xytext=(-8, 0), textcoords="offset points", ha="right", va="center",
                    color=INK_2, fontsize=9)
    for u, outcome in other:        # an outcome, not a gap in the data
        ax.axvline(u, color=MUTED, linestyle=":", linewidth=1.2)
        ax.annotate(outcome.replace("_", " "), (u, 0.02), xycoords=("data", "axes fraction"),
                    rotation=90, color=INK_2, fontsize=9, ha="right", va="bottom")
    top = max([L for _, L in measured + bounded], default=1.0)
    ax.set_ylim(0, top + (6.0 if bounded else 1.0))   # room for the arrows
    ax.legend(loc="lower right", frameon=False, fontsize=12, labelcolor=INK_2)
    out = HERE / "flame_length_vs_inlet_velocity.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig); print("wrote", out.name)


def _read_field(path: Path) -> "list[float]":
    """The internalField of an ASCII OpenFOAM scalar field written per cell."""
    lines = path.read_text().splitlines()
    i = next(k for k, line in enumerate(lines) if line.startswith("internalField"))
    n = int(lines[i + 1])
    return [float(v) for v in lines[i + 3:i + 3 + n]]


def plot_flames() -> None:
    """The slot-burner flames standing upright, with their length against velocity.

    Needs the LvsUlong cases with their cell centres written
    (``postProcess -func writeCellCentres -latestTime`` in each case). The domain
    is half a slot, so each field is mirrored about the slot centreline to show
    the whole flame. The colour scale stops at 1800 K: the burnt gas varies by a
    few tens of kelvin, and finer bands showed that as islands above the tip.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import numpy as np
    plan = json.loads((REPO / FLAME_CASES / "campaign.json").read_text())
    rows = list(csv.DictReader(open(LENGTH_CSV)))
    cases = sorted(plan["cases"], key=lambda e: e["value"])
    fig = plt.figure(figsize=(6.4, 8.0), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    fig.text(0.04, 0.975, "Slot-burner flames", color=INK, fontsize=18,
             fontweight="bold", va="top")
    fig.text(0.04, 0.943, "Methane–air, φ = 0.7 · five cases run by the agent",
             color=INK_2, fontsize=12.5, va="top")
    below, above = 2.0, 25.0                 # mm around the plate's trailing edge
    width = 0.92 / len(cases)
    image = None
    for i, entry in enumerate(cases):
        ax = fig.add_axes([0.04 + i * width, 0.47, width * 0.94, 0.43])
        case = REPO / "testcases" / entry["case_name"]
        latest = max((d for d in case.iterdir() if d.name.replace(".", "", 1).isdigit()
                      and (d / "Cx").is_file()), key=lambda d: float(d.name))
        x, y, T = (np.array(_read_field(latest / n)) for n in ("Cx", "Cy", "T"))
        h = y.max() + (y.max() - np.sort(np.unique(y))[-2]) / 2   # the slot centreline
        # The plate is the gap in the bottom row of cells; lengths start at its end.
        bottom = np.sort(x[np.isclose(y, y.min())])
        edge_te = bottom[np.argmax(np.diff(bottom)) + 1]
        xs, ys, Ts = np.r_[x, x], np.r_[y, 2 * h - y], np.r_[T, T]
        # Standing up: across the slot horizontally, the flow upwards.
        tri = mtri.Triangulation((ys - h) * 1e3, (xs - edge_te) * 1e3)
        tx, ty = tri.x[tri.triangles], tri.y[tri.triangles]
        span = np.max(np.abs(tx - np.roll(tx, 1, axis=1)) + np.abs(ty - np.roll(ty, 1, axis=1)), axis=1)
        tri.set_mask(span > 0.5)            # no triangles across the plates
        image = ax.tricontourf(tri, np.clip(Ts, 280, 1800), levels=np.linspace(280, 1800, 25),
                               cmap="inferno", extend="max")
        ax.set_facecolor("#b9b8b1")          # the burner plates
        ax.set_xlim(-h * 1e3, h * 1e3)
        ax.set_ylim(-below - 3.5, above)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_title(f"{entry['value']:g} m/s", color=INK, fontsize=12.5,
                     fontweight="bold", y=-0.075, va="top")
    # Flame length against inlet velocity, under the flames it measures.
    measured = [(float(r["inlet_velocity_m_s"]), float(r["flame_length_mm"]))
                for r in rows if r["flame_length_mm"]]
    ax = fig.add_axes([0.16, 0.11, 0.78, 0.25])
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    ax.plot(*zip(*measured), color=SERIES[0], linewidth=3, marker="o", markersize=10,
            markeredgecolor=SURFACE, markeredgewidth=2)
    for u, L in measured:
        ax.annotate(f"{L:.1f} mm", (u, L), xytext=(0, 13), textcoords="offset points",
                    ha="center", color=INK_2, fontsize=11)
    ax.tick_params(colors=MUTED, labelsize=12, length=0, pad=5)
    ax.set_xlabel("Inlet velocity  [m/s]", color=INK_2, fontsize=13)
    ax.set_ylabel("Flame length  [mm]", color=INK_2, fontsize=13)
    ax.set_ylim(0, max(L for _, L in measured) * 1.25)
    cax = fig.add_axes([0.16, 0.425, 0.78, 0.012])
    bar = fig.colorbar(image, cax=cax, orientation="horizontal", ticks=[300, 800, 1300, 1800])
    bar.ax.set_xticklabels(["300", "800", "1300", "≥1800 K"])
    bar.ax.tick_params(labelsize=11, colors=MUTED, length=0)
    bar.outline.set_visible(False)
    fig.savefig(FLAMES_PNG, facecolor=SURFACE); plt.close(fig); print("wrote", FLAMES_PNG.name)


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "plot"
    if step in ("collect", "all"):
        collect()
    if step in ("flames", "collect", "all"):
        plot_flames()          # needs the cases themselves, like collect
    if step in ("plot", "all"):
        plot_speed(); plot_length()
