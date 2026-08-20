#!/usr/bin/env python
"""Phase 5 / FEEMS step 6: publication-ready figures. STYLING ONLY.

Reads the finished artefacts and re-draws them. It refits NOTHING and recomputes NO statistic:
edge weights, bootstrap CIs and the supported/not-supported call are read verbatim from
`<tag>_<scheme>_bootstrap.npz`, so the numbers in results/phase5/feems_report.md cannot move.

NON-NEGOTIABLE STRUCTURE
  Every figure is TWO panels: (a) the point-estimate surface, (b) only those edges whose
  bootstrap CI excludes zero. Panel (b) is near-empty (1/329 within-Mf, 0/329 all-cluster) and
  that emptiness IS the result. The point estimate is never emitted on its own, because on its
  own it implies a barrier the data do not support.

STYLE
  Basemap and layout follow the project's existing Sabah maps (figures/sabah_admix_map.*,
  tess_map.*, clonal_cluster_plot.*), which are ggplot2/svglite at 720x504 pt, Arial,
  8.8 pt text / 11 pt titles, black coastline at 1.07 pt, plain long/lat axes, no graticule,
  legend outside the panel on the right.
  The diverging orange<->blue ramp and the w/wbar colourbar are the FEEMS convention and are
  deliberately NOT restyled to the project palette -- the 13 stops are copied from
  feems/viz.py so the colourway is bit-identical to the standard.

Outputs (svg + pdf + 300-dpi png):
  figures/phase5/feems_Mf_only.*        primary, within-Mf
  figures/phase5/feems_all_cluster.*    secondary, all clusters

Run: micromamba run -n feems_e python scripts/2x_feems_figures.py
"""
from __future__ import annotations

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as clr  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# --- house style, read off figures/sabah_admix_map.svg et al. ------------------
FONT = "Arial"
TXT, TITLE = 8.8, 11.0          # svglite: axis/legend text 8.80 pt, titles 11.00 pt
LW_FRAME = 1.07                 # panel border / coastline stroke-width in the reference SVGs
LW_THIN = 0.43
COAST = "#000000"

# FEEMS convention -- do not restyle. Copied verbatim from feems/viz.py.
EEMS_COLORS = ["#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC",
               "#FBFBFB", "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99"]
EDGE_CMAP = clr.LinearSegmentedColormap.from_list("eems_colors", EEMS_COLORS, N=256)

UNSUPPORTED = "#d9d9d9"         # lattice drawn but greyed: present, not distinguishable
DEME_FACE = "#bdbdbd"

FIGS = [
    {"out": "feems_Mf_only", "tag": "Mf_only", "scheme": "individual",
     "cohort": "within-Mf", "panel_a": "(a) fitted effective-migration surface"},
    {"out": "feems_all_cluster", "tag": "cleaned_501", "scheme": "individual",
     "cohort": "all clusters (Mf + Mn)", "panel_a": "(a) fitted effective-migration surface"},
]


def set_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT, "Helvetica", "DejaVu Sans"],
        "font.size": TXT,
        "axes.labelsize": TITLE,
        "axes.titlesize": TITLE,
        "xtick.labelsize": TXT,
        "ytick.labelsize": TXT,
        "legend.fontsize": TXT,
        "legend.title_fontsize": TITLE,
        "axes.linewidth": LW_FRAME,
        "xtick.major.width": LW_THIN,
        "ytick.major.width": LW_THIN,
        "xtick.major.size": 2.75,
        "ytick.major.size": 2.75,
        "axes.grid": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        # mathtext (lambda, w/wbar, exponents) must not fall back to DejaVu in the vector output
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT,
        "mathtext.it": f"{FONT}:italic",
        "mathtext.bf": f"{FONT}:bold",
        "svg.fonttype": "none",      # keep text as text in the vector output
        "pdf.fonttype": 42,
    })


def load_coastline(shp):
    import fiona
    from shapely.geometry import shape

    rings = []
    for rec in fiona.open(shp):
        g = shape(rec["geometry"])
        for poly in (g.geoms if g.geom_type == "MultiPolygon" else [g]):
            rings.append(np.asarray(poly.exterior.coords))
            rings.extend(np.asarray(i.coords) for i in poly.interiors)
    return rings


def deme_sizes(n, smin=6.0, smax=90.0, ncap=50):
    """Marker area ~ sqrt(n), capped so the 50-sample deme stays legible without dominating."""
    n = np.asarray(n, dtype=float)
    f = np.sqrt(np.clip(n, 1, ncap)) / np.sqrt(ncap)
    return smin + (smax - smin) * f


def draw_panel(ax, rings, ends, r, mask, demes, abs_max, title, edge_lw=0.9,
               emphasise=False):
    for ring in rings:
        ax.plot(ring[:, 0], ring[:, 1], color=COAST, lw=LW_FRAME,
                solid_capstyle="round", solid_joinstyle="round", zorder=6)

    norm = clr.Normalize(vmin=-abs_max, vmax=abs_max)
    from matplotlib.collections import LineCollection

    def segs_of(idx):
        return [[(ends[i, 0], ends[i, 1]), (ends[i, 2], ends[i, 3])] for i in idx]

    # unsupported underneath, greyed: the lattice is there, the signal on it is not
    off = np.where(~mask)[0]
    if off.size:
        ax.add_collection(LineCollection(segs_of(off), colors=[UNSUPPORTED] * off.size,
                                         linewidths=edge_lw * 0.8, zorder=3, capstyle="round"))
    on = np.where(mask)[0]
    if on.size:
        # In the supported-only panel a lone surviving edge can be almost invisible if its
        # weight is near the surface mean. Give supported edges a thin dark casing so they are
        # FINDABLE; colour still encodes the magnitude, which is not exaggerated.
        if emphasise and on.size < len(mask):
            ax.add_collection(LineCollection(segs_of(on), colors=["#3f3f3f"] * on.size,
                                             linewidths=edge_lw * 3.0, zorder=4,
                                             capstyle="round"))
        ax.add_collection(LineCollection(segs_of(on), colors=[EDGE_CMAP(norm(r[i])) for i in on],
                                         linewidths=edge_lw * (1.9 if emphasise else 1.0),
                                         zorder=5, capstyle="round"))

    ax.scatter(demes[:, 0], demes[:, 1], s=deme_sizes(demes[:, 2]),
               facecolor=DEME_FACE, edgecolor="black", linewidths=LW_THIN,
               zorder=7, clip_on=True)

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title, fontsize=TITLE, pad=6)
    return norm


def common_abs_max(cfgs, in_dir):
    """One colour scale for every figure, so the two are directly comparable."""
    m = 0.0
    for c in cfgs:
        z = np.load(os.path.join(in_dir, f"{c['tag']}_{c['scheme']}_bootstrap.npz"))
        m = max(m, float(np.abs(z["r_hat"]).max()))
    return float(np.clip(np.round(m + 0.05, 1), 0.2, 2.0))


def build(cfg, args, abs_max):
    stem = os.path.join(args.in_dir, cfg["tag"])
    boot_stem = os.path.join(args.in_dir, f"{cfg['tag']}_{cfg['scheme']}")

    with open(stem + "_fit.json") as fh:
        fit = json.load(fh)
    with open(os.path.join(args.res_dir, f"{cfg['tag']}_{cfg['scheme']}_bootstrap.json")) as fh:
        bj = json.load(fh)
    bz = np.load(boot_stem + "_bootstrap.npz")

    r = bz["r_hat"]
    mask = bz["supported"].astype(bool)
    ends = bz["edge_endpoints"]
    demes = np.loadtxt(stem + "_demes.tsv", delimiter="\t", skiprows=1)
    rings = load_coastline(args.coastline)

    n = fit["inputs"]["n_samples"]
    n_snp = fit["inputs"]["n_polymorphic_sites"]
    lamb = fit["lamb_chosen"]
    n_sup, n_edges, n_boot = bj["n_supported"], bj["n_edges"], bj["n_boot_usable"]
    n_demes_obs = fit["graph"]["n_observed_demes"]
    folds = fit["cv"]["n_folds"] if fit.get("cv") else None
    set_style()
    # Panel geometry is computed in inches from the data aspect so the map fills its frame
    # and the figure has no dead space -- the reference maps are 10 in wide, so we keep that.
    ends_x, ends_y = ends[:, [0, 2]], ends[:, [1, 3]]
    x0, x1, y0, y1 = ends_x.min(), ends_x.max(), ends_y.min(), ends_y.max()
    padx, pady = 0.10 * (x1 - x0), 0.06 * (y1 - y0)
    x0, x1, y0, y1 = x0 - padx, x1 + padx, y0 - pady, y1 + pady
    asp = 1.0 / np.cos(np.radians(0.5 * (y0 + y1)))          # quickmap, as in the reference maps

    FW = 10.0
    L, GAP, LEGW, RM = 0.62, 0.34, 1.58, 0.12
    PW = (FW - L - GAP - LEGW - RM) / 2.0
    PH = PW * ((y1 - y0) * asp) / (x1 - x0)
    BOT, TOPGAP = 1.62, 0.76
    FH = BOT + PH + TOPGAP

    fig = plt.figure(figsize=(FW, FH))
    def axes_at(x_in, w_in):
        return fig.add_axes([x_in / FW, BOT / FH, w_in / FW, PH / FH])
    ax1 = axes_at(L, PW)
    ax2 = axes_at(L + PW + GAP, PW)

    panel_b = (f"(b) edges supported by bootstrap — {n_sup}/{n_edges} at 95% CI"
               if n_sup else f"(b) edges supported by bootstrap — none, 0/{n_edges} at 95% CI")
    norm = draw_panel(ax1, rings, ends, r, np.ones(len(r), bool), demes, abs_max, cfg["panel_a"])
    draw_panel(ax2, rings, ends, r, mask, demes, abs_max, panel_b, emphasise=True)

    for ax in (ax1, ax2):
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect(asp)
        ax.set_xticks(np.arange(np.ceil(x0), np.floor(x1) + 1, 1))
        ax.set_yticks(np.arange(np.ceil(y0), np.floor(y1) + 1, 1))
        ax.tick_params(direction="out", top=False, right=False)
    ax2.set_ylabel("")
    plt.setp(ax2.get_yticklabels(), visible=False)

    # ---- legend column: FEEMS colourbar + deme size key ----------------------
    LX = L + 2 * PW + GAP + 0.30
    cax = fig.add_axes([LX / FW, (BOT + PH * 0.52) / FH, 0.17 / FW, (PH * 0.33) / FH])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=EDGE_CMAP), cax=cax,
                      orientation="vertical")
    cb.set_ticks([-abs_max, 0, abs_max])
    cb.set_ticklabels([f"$10^{{-{abs_max:g}}}$", "$10^{0}$", f"$10^{{{abs_max:g}}}$"])
    cb.outline.set_linewidth(LW_THIN)
    cb.ax.tick_params(width=LW_THIN, length=2.5, labelsize=TXT, pad=2)
    fig.text(LX / FW, (BOT + PH * 1.00) / FH, "$w\\,/\\,\\bar{w}$",
             fontsize=TITLE, ha="left", va="top")
    fig.text(LX / FW, (BOT + PH * 0.945) / FH, "effective migration",
             fontsize=TXT, ha="left", va="top", color="#444444")

    key_n = [1, 5, 10, 25, 50]
    handles = [Line2D([], [], marker="o", linestyle="none",
                      markerfacecolor=DEME_FACE, markeredgecolor="black",
                      markeredgewidth=LW_THIN,
                      markersize=np.sqrt(deme_sizes([k])[0]), label=str(k)) for k in key_n]
    lax = fig.add_axes([LX / FW, BOT / FH, LEGW / FW, (PH * 0.44) / FH])
    lax.set_axis_off()
    leg = lax.legend(handles=handles, title="Samples per deme", loc="upper left",
                     bbox_to_anchor=(0.0, 1.0), frameon=False, labelspacing=0.75,
                     handletextpad=0.8, borderpad=0.0, borderaxespad=0.0)
    leg._legend_box.align = "left"

    # ---- title and the honest caption block ----------------------------------
    import textwrap

    fig.text(L / FW, 1.0 - 0.20 / FH,
             f"FEEMS effective migration across Sabah — {cfg['cohort']} "
             f"(n = {n}, {n_snp:,} polymorphic sites)",
             fontsize=TITLE, ha="left", va="top")

    cv_txt = f"leave-one-deme-out CV ({folds} folds)" if folds else "CV"
    foot = (
        f"$\\lambda$ = {lamb:g} chosen by {cv_txt}; {n_demes_obs} observed demes on a 0.30° "
        f"triangular lattice ({n_edges} edges). Panel (b): bootstrap over individuals, "
        f"{n_boot} replicates; edges supported at 95% CI: {n_sup}/{n_edges}. "
        f"The surface is near-uniform and no barrier is statistically supported — grey edges in "
        f"(b) are present in the lattice but indistinguishable from uniform migration. "
        f"Power-limited at this n; not evidence of absence. Colour scale is shared across figures."
    )
    wrapped = "\n".join(textwrap.wrap(foot, width=150, break_long_words=False))
    fig.text(L / FW, 0.16 / FH, wrapped, fontsize=TXT, ha="left", va="bottom", linespacing=1.6)

    os.makedirs(args.fig_dir, exist_ok=True)
    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(os.path.join(args.fig_dir, f"{cfg['out']}.{ext}"), **kw)
    plt.close(fig)

    print(f"  {cfg['out']}: n={n} lambda={lamb:g} supported={n_sup}/{n_edges} "
          f"(from {n_boot} {cfg['scheme']} replicates) -> svg/pdf/png")
    return {"out": cfg["out"], "tag": cfg["tag"], "cohort": cfg["cohort"], "n_samples": n,
            "n_polymorphic_sites": n_snp, "lamb": lamb, "cv_folds": folds,
            "n_observed_demes": n_demes_obs, "n_edges": n_edges, "n_supported": n_sup,
            "n_boot": n_boot, "abs_max": abs_max}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", default="data/processed/feems")
    p.add_argument("--res-dir", default="results/phase5")
    p.add_argument("--fig-dir", default="figures/phase5")
    p.add_argument("--coastline", default="data/raw/MSC/spatial/malaysia.shp")
    args = p.parse_args()

    print("=== publication figures (styling only; no refit, no recomputation) ===")
    abs_max = common_abs_max(FIGS, args.in_dir)
    print(f"  shared colour scale: +/-{abs_max:g} log10 units")
    summary = [build(cfg, args, abs_max) for cfg in FIGS]
    with open(os.path.join(args.res_dir, "feems_figures_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"=== wrote {args.res_dir}/feems_figures_summary.json ===")


if __name__ == "__main__":
    main()
