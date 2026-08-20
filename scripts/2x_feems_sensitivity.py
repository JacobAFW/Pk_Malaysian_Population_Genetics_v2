#!/usr/bin/env python
"""Phase 5 / FEEMS step 5: is the surface an artefact of the grid or of one over-sampled deme?

Refits the same data under a set of configurations and asks whether a named spatial feature
survives all of them. Each configuration re-runs its own CV, because the right lambda depends
on the deme resolution -- reusing the primary lambda across grids would beg the question.

Configurations:
  - grid spacing 0.50 / 0.30 / 0.20 degrees   (coarser and finer than the primary)
  - max 10 and max 5 samples per deme at 0.30  (the largest deme holds 50 of 276/335)

The feature under test defaults to the northern-Sabah low-migration band, whose box was read
off the primary fit's eight lowest-weight edges (midpoints 116.9-117.3 E, 6.5-6.9 N) and then
padded. It is scored as a CONTRAST:

    band_contrast = mean(r over edges with midpoint inside the box) - mean(r over all edges)

with r = log10(w) - mean(log10 w). A real barrier gives a strongly negative contrast in every
configuration. The contrast is used rather than raw weights because grids differ between
configurations, so individual edges are not comparable but a regional mean is.

Outputs:
  results/phase5/<tag>_sensitivity.json
  figures/phase5/<tag>_feems_sensitivity.png    one surface panel per configuration

Run: micromamba run -n feems_e python scripts/2x_feems_sensitivity.py --tag Mf_only
"""
from __future__ import annotations

import argparse
import io
import json
import os
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

CONFIGS = [  # (label, grid spacing in degrees, max samples per deme or None)
    ("grid 0.50 deg", 0.50, None),
    ("grid 0.30 deg (primary)", 0.30, None),
    ("grid 0.20 deg", 0.20, None),
    ("0.30 deg, max 10/deme", 0.30, 10),
    ("0.30 deg, max 5/deme", 0.30, 5),
]


def load_modules():
    import importlib.util as ilu

    here = os.path.dirname(os.path.abspath(__file__))
    mods = {}
    for name, fname in (("fit", "2x_feems_fit.py"), ("boot", "2x_feems_bootstrap.py")):
        spec = ilu.spec_from_file_location(f"{name}_mod", os.path.join(here, fname))
        m = ilu.module_from_spec(spec)
        spec.loader.exec_module(m)
        mods[name] = m
    mods["fit"].feems_compat.patch_cholmod(verbose=False)
    return mods["fit"], mods["boot"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="Mf_only")
    p.add_argument("--in-dir", default="data/processed/feems")
    p.add_argument("--fig-dir", default="figures/phase5")
    p.add_argument("--out-dir", default="results/phase5")
    p.add_argument("--band-box", default="116.8,117.4,6.4,7.0",
                   help="long_min,long_max,lat_min,lat_max of the feature under test")
    p.add_argument("--band-name", default="northern Sabah low-migration band")
    p.add_argument("--lamb-grid", default="1e-6,1e6,17")
    p.add_argument("--cv-folds", type=int, default=10)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.fig_dir, exist_ok=True)
    fit, boot = load_modules()
    from feems.cross_validation import run_cv

    stem_in = os.path.join(args.in_dir, args.tag)
    with open(stem_in + "_fit.json") as fh:
        meta = json.load(fh)
    scale_snps = bool(meta["inputs"]["scale_snps"])

    x0, x1, y0, y1 = (float(v) for v in args.band_box.split(","))
    print(f"=== sensitivity for {args.tag}; feature = {args.band_name} "
          f"[{x0}-{x1} E, {y0}-{y1} N] ===")

    geno0, ids0 = fit.load_genotypes(stem_in + "_feems")
    geno0, _, _ = fit.mean_impute(geno0)
    geno0 = fit.assert_polymorphic(geno0, label="sensitivity")
    coord_by_id = {}
    with open(stem_in + "_coords.tsv") as fh:
        for line in fh:
            iid, lo, la = line.rstrip("\n").split("\t")
            coord_by_id[iid] = (float(lo), float(la))
    coord0 = np.array([coord_by_id[i] for i in ids0])
    outer = np.loadtxt(stem_in + "_outer.csv", delimiter=",", skiprows=1)

    lo_, hi_, n_ = args.lamb_grid.split(",")
    lamb_grid = np.geomspace(float(lo_), float(hi_), int(n_))[::-1]

    rows, panels = [], []
    for label, spacing, cap in CONFIGS:
        print(f"--- {label} ---")
        grid, edges = fit.make_tri_grid(outer, spacing)
        geno, coord = geno0, coord0
        n_dropped = 0
        if cap:
            keep, assign = fit.thin_demes(coord0, grid, cap)
            n_dropped = int((~keep).sum())
            geno, coord = geno0[keep], coord0[keep]
            geno = fit.assert_polymorphic(geno, label=f"  thin{cap}")
        sp_graph = boot.quiet_graph(fit, geno, coord, grid, edges, scale_snps)

        saved = sys.stdout
        sys.stdout = io.StringIO()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cv_err = run_cv(sp_graph, lamb_grid, n_folds=args.cv_folds or None,
                                factr=1e10, outer_verbose=False, inner_verbose=False)
        finally:
            sys.stdout = saved
        mean_cv = np.mean(cv_err, axis=(0, 2))
        lamb = float(lamb_grid[np.argmin(mean_cv)])
        cv_uniform = float(mean_cv[0])          # lamb_grid is descending: [0] is the largest lamb
        gain = (cv_uniform - mean_cv.min()) / cv_uniform

        sp_graph.fit(lamb=lamb, optimize_q=None, verbose=False)
        w = boot.w_by_edge(sp_graph, edges, grid.shape[0])
        r = boot.rel_log_weights(w)
        e = np.asarray(edges)
        mid = (grid[e[:, 0] - 1] + grid[e[:, 1] - 1]) / 2
        in_band = (mid[:, 0] >= x0) & (mid[:, 0] <= x1) & (mid[:, 1] >= y0) & (mid[:, 1] <= y1)
        contrast = float(r[in_band].mean() - r.mean()) if in_band.any() else float("nan")
        obs, r_full = fit.deme_distances(sp_graph)
        from scipy.stats import pearsonr
        fit_r = float(pearsonr(obs, r_full)[0])

        print(f"  n={geno.shape[0]} (dropped {n_dropped})  demes={grid.shape[0]} "
              f"observed={sp_graph.n_observed_nodes}  lambda={lamb:g}")
        print(f"  CV gain over uniform = {gain * 100:.3f}%   fit r = {fit_r:.4f}   "
              f"sd(r) = {r.std():.3f}")
        print(f"  edges in band = {int(in_band.sum())}   BAND CONTRAST = {contrast:+.3f} log10 units")

        rows.append({
            "label": label, "grid_spacing_deg": spacing, "max_per_deme": cap,
            "n_samples": int(geno.shape[0]), "n_samples_dropped": n_dropped,
            "n_demes": int(grid.shape[0]), "n_edges": int(e.shape[0]),
            "n_observed_demes": int(sp_graph.n_observed_nodes),
            "lamb_chosen": lamb,
            "cv_err_min": float(mean_cv.min()), "cv_err_uniform": cv_uniform,
            "cv_gain_over_uniform_frac": float(gain),
            "fit_r": fit_r, "sd_r": float(r.std()),
            "r_range": [float(r.min()), float(r.max())],
            "n_edges_in_band": int(in_band.sum()),
            "band_contrast": contrast,
        })
        panels.append((label, sp_graph, w, grid.shape[0], edges, lamb, contrast))

    # ---- figure --------------------------------------------------------------
    abs_max = 1.0
    ncol = 3
    nrow = int(np.ceil(len(panels) / ncol))
    fig = plt.figure(figsize=(5.0 * ncol, 4.0 * nrow), dpi=200)
    for k, (label, sg, w, n_nodes, edges, lamb, contrast) in enumerate(panels):
        idx = boot.edge_to_w_index(sg, edges, n_nodes)
        w_draw = np.empty(len(w))
        w_draw[idx] = w
        boot.draw_surface_on_ax(fig, (nrow, ncol, k + 1), sg, w_draw,
                                f"{label}\n$\\lambda$={lamb:g}, band contrast = {contrast:+.3f}",
                                abs_max)
    fig.suptitle(f"{args.tag}: grid-resolution and sampling-balance sensitivity\n"
                 f"(colour scale fixed at +/-{abs_max} log10 units for comparability)", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(args.fig_dir, f"{args.tag}_feems_sensitivity.{ext}"),
                    bbox_inches="tight")
    plt.close(fig)

    contrasts = [r["band_contrast"] for r in rows]
    verdict = ("survives: negative in every configuration"
               if all(c < -0.1 for c in contrasts) else
               "DOES NOT survive: sign or magnitude changes across configurations")
    print(f"=== band contrast across configurations: "
          f"{', '.join(f'{c:+.3f}' for c in contrasts)} -> {verdict} ===")

    with open(os.path.join(args.out_dir, f"{args.tag}_sensitivity.json"), "w") as fh:
        json.dump({"tag": args.tag, "band_name": args.band_name,
                   "band_box": [x0, x1, y0, y1], "cv_folds": args.cv_folds,
                   "configurations": rows, "verdict": verdict}, fh, indent=2)
    print(f"=== wrote {args.out_dir}/{args.tag}_sensitivity.json and "
          f"{args.fig_dir}/{args.tag}_feems_sensitivity.png ===")


if __name__ == "__main__":
    main()
