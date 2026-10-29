"""Regenerate the example figures in docs/examples/ from agent campaign output.

    agent/.venv/bin/python docs/examples/make_figures.py collect   # CSVs from campaigns
    agent/.venv/bin/python docs/examples/make_figures.py plot      # PNGs from the CSVs

`collect` needs the campaigns under agent/campaigns/ (git-ignored); `plot` needs
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


def collect() -> None:
    import cantera_tool
    rows = list(csv.DictReader(open(REPO / "agent/campaigns/SLvsT/results.csv")))
    cache = json.loads(cantera_tool._FLAME_CACHE.read_text())
    with SPEED_CSV.open("w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["T_unburnt_K", "S_L_FGMFoam_m_s", "S_L_Cantera_UL_m_s", "difference_pct"])
        for r in sorted(rows, key=lambda r: float(r["T_inlet"])):
            T, s = float(r["T_inlet"]), float(r["laminar_flame_speed"])
            c = cache[cantera_tool._cache_key(0.7, T, 101325.0)]["flame_speed"]
            out.writerow([f"{T:g}", f"{s:.4f}", f"{c:.4f}", f"{100 * (s - c) / c:+.1f}"])
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
    fig, ax = plt.subplots(figsize=(7.0, 4.4), dpi=200)
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9, length=0, pad=6)
    ax.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=INK_2, fontsize=10, labelpad=8)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=10, labelpad=8)
    fig.text(0.075, 0.955, title, color=INK, fontsize=13, fontweight="bold", va="top")
    fig.text(0.075, 0.895, subtitle, color=INK_2, fontsize=9.5, va="top")
    fig.subplots_adjust(left=0.1, right=0.83, top=0.80, bottom=0.14)
    return fig, ax, plt


def plot_speed() -> None:
    rows = list(csv.DictReader(open(SPEED_CSV)))
    T = [float(r["T_unburnt_K"]) for r in rows]
    fgm = [float(r["S_L_FGMFoam_m_s"]) for r in rows]
    ref = [float(r["S_L_Cantera_UL_m_s"]) for r in rows]
    fig, ax, plt = _axes(
        "Laminar flame speed against unburnt temperature",
        "Methane–air, φ = 0.7, 1 atm · FGMFoam measured by the agent from a 9-case campaign",
        "Unburnt temperature  [K]", "Laminar flame speed  [m/s]")
    ax.plot(T, ref, color=SERIES[1], linewidth=2, linestyle=(0, (5, 3)),
            label="Cantera 1D flame, unity Lewis, no Soret (GRI-3.0)")
    ax.plot(T, fgm, color=SERIES[0], linewidth=2, marker="o", markersize=6.5,
            markeredgecolor=SURFACE, markeredgewidth=1.5,
            label="FGMFoam (tabulated, measured by the agent)")
    # Direct labels at the line ends, in ink — the line beside them carries colour.
    # Nudge them apart vertically when the line ends are close, so they never collide.
    upper, lower = ("Cantera", ref[-1]), ("FGMFoam", fgm[-1])
    if lower[1] > upper[1]:
        upper, lower = lower, upper
    ax.annotate(upper[0], (T[-1], upper[1]), xytext=(8, 7), textcoords="offset points",
                va="center", color=INK, fontsize=9.5)
    ax.annotate(lower[0], (T[-1], lower[1]), xytext=(8, -7), textcoords="offset points",
                va="center", color=INK, fontsize=9.5)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=INK_2)
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
        "Flame length against inlet velocity",
        "2D laminar slot burner, methane–air, φ = 0.7, 1 atm · measured by the FGMFoam agent",
        "Inlet velocity  [m/s]", "Flame length  [mm]")
    if measured:
        ax.plot(*zip(*measured), color=SERIES[0], linewidth=2, marker="o", markersize=6.5,
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                label="flame length (tip on the slot centreline, 900 K)")
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
    ax.legend(loc="lower right", frameon=False, fontsize=9, labelcolor=INK_2)
    out = HERE / "flame_length_vs_inlet_velocity.png"
    fig.savefig(out, facecolor=SURFACE); plt.close(fig); print("wrote", out.name)


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "plot"
    if step in ("collect", "all"):
        collect()
    if step in ("plot", "all"):
        plot_speed(); plot_length()
