#!/usr/bin/env python
"""Phase 5 / FEEMS step 2: build the spatial graph, sweep lambda, fit the effective-migration surface.

Inputs (from scripts/2x_feems_prep.sh):
  data/processed/feems/<tag>_feems.{bed,bim,fam}   monomorphic-free genotypes
  data/processed/feems/<tag>_coords.tsv            IID -> long/lat  (joined by ID, not position)
  data/processed/feems/<tag>_outer.csv             outer region ring

Outputs:
  data/processed/feems/<tag>_fit.json              counts, sweep, CV, chosen lambda
  data/processed/feems/<tag>_weights.npz           fitted edge weights per lambda
  figures/phase5/<tag>_feems_surface_lamb*.png/svg
  figures/phase5/<tag>_feems_cv.png
  figures/phase5/<tag>_feems_fit_vs_dist.png

THE #30 CONVERGENCE FIX
  feems standardises each SNP by 1/sqrt(mu*(1-mu)) when scale_snps=True. A monomorphic
  site gives mu == 0 -> division by zero -> nan/inf in S -> "did not converge" from
  L-BFGS. Prep drops them with `plink --maf`; `assert_polymorphic()` below re-checks the
  actual matrix handed to SpatialGraph and aborts loudly if any zero-variance column
  survives (e.g. one created by mean-imputation of an all-missing column).

Run: micromamba run -n feems_e python scripts/2x_feems_fit.py --tag cleaned_501
"""
from __future__ import annotations

import argparse
import json
import os
import time
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from pandas_plink import read_plink1_bin  # noqa: E402

import sys  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Must patch before feems.spatial_graph is used: scikit-sparse 0.5 dropped the
# callable-Factor API feems 1.0.0 relies on. See 2x_feems_compat.py.
import importlib.util as _ilu  # noqa: E402

_spec = _ilu.spec_from_file_location(
    "feems_compat", os.path.join(os.path.dirname(os.path.abspath(__file__)), "2x_feems_compat.py"))
feems_compat = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(feems_compat)

from feems import SpatialGraph, Viz  # noqa: E402
from feems.cross_validation import run_cv  # noqa: E402
from feems.utils import prepare_graph_inputs  # noqa: E402

FEEMS_GRIDS = {  # bundled discrete global grids, coarse -> fine
    "500": "grid_500.shp",
    "250": "grid_250.shp",
    "100": "grid_100.shp",
}


def feems_grid_path(res):
    import feems

    return os.path.join(os.path.dirname(feems.__file__), "data", FEEMS_GRIDS[res])


def load_genotypes(bed_base):
    """Read PLINK1 -> (genotypes float64 (n, p) with NaN, sample IDs in .fam order)."""
    G = read_plink1_bin(bed_base + ".bed", bed_base + ".bim", bed_base + ".fam", verbose=False)
    geno = np.asarray(G.values, dtype=np.float64)
    ids = [str(s) for s in G.sample.values]
    return geno, ids


def mean_impute(geno):
    """Replace NaN with the per-SNP mean. feems' allele-frequency step uses np.mean and
    cannot tolerate NaN. Returns (imputed, per-SNP missing rate, n all-missing columns)."""
    miss = np.isnan(geno)
    rate = miss.mean(axis=0)
    all_missing = int((rate == 1.0).sum())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        col_mean = np.nanmean(geno, axis=0)
    col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
    out = np.where(miss, col_mean[None, :], geno)
    return out, rate, all_missing


def assert_polymorphic(geno, label=""):
    """THE #30 FIX, re-checked on the real matrix. Drop and report zero-variance columns."""
    var = geno.var(axis=0)
    freq = geno.mean(axis=0) / 2.0          # allele frequency in [0, 1]
    bad = (var <= 0) | (freq <= 0) | (freq >= 1)
    n_bad = int(bad.sum())
    if n_bad:
        print(f"  {label}: dropping {n_bad} non-polymorphic column(s) missed upstream")
        geno = geno[:, ~bad]
    var = geno.var(axis=0)
    assert geno.shape[1] > 0, "no polymorphic sites left"
    assert var.min() > 0, f"{label}: {int((var <= 0).sum())} zero-variance columns remain"
    # feems computes mu = mean(freq_per_deme)/2, then divides by sqrt(mu*(1-mu));
    # report the sample-level analogue so a near-zero denominator is visible.
    f = geno.mean(axis=0) / 2.0
    print(f"  {label}: n_polymorphic_sites={geno.shape[1]}  min_var={var.min():.3e}  "
          f"freq range=[{f.min():.5f}, {f.max():.5f}]")
    return geno


def make_tri_grid(outer, spacing):
    """Triangular (hexagonal-connectivity) deme lattice clipped to the outer ring.

    feems ships grid_100/250/500.shp (ISEA3H global grids); grid_100 is the only complete
    one in this install and yields ~10 demes over Sabah, far too coarse to resolve a
    migration surface. So build the lattice directly at an explicit, documentable spacing:
    rows are offset by half a step and separated by spacing*sqrt(3)/2, giving each interior
    node the 6 equidistant neighbours EEMS/FEEMS assume.

    Returns (grid (n,2) long/lat, edges (m,2) 1-indexed).
    """
    from shapely.geometry import Point, Polygon
    from shapely.prepared import prep

    poly = prep(Polygon(outer))
    x0, y0, x1, y1 = Polygon(outer).bounds
    dy = spacing * np.sqrt(3.0) / 2.0
    idx, pts = {}, []
    for j in range(int(np.floor((y1 - y0) / dy)) + 2):
        y = y0 + j * dy
        xoff = (j % 2) * spacing / 2.0
        for i in range(int(np.floor((x1 - x0) / spacing)) + 2):
            x = x0 + xoff + i * spacing
            if poly.intersects(Point(x, y)):
                idx[(i, j)] = len(pts)
                pts.append((x, y))
    grid = np.array(pts)
    assert grid.shape[0] > 0, "empty grid; reduce --grid-spacing"

    edges = set()
    for (i, j), a in idx.items():
        # right neighbour, plus the two neighbours in the row above (offset by parity)
        up = (0, 1) if j % 2 == 0 else (1, 1)
        for di, dj in ((1, 0), (up[0] - 1, 1), (up[0], 1)):
            b = idx.get((i + di, j + dj))
            if b is not None:
                edges.add((min(a, b) + 1, max(a, b) + 1))
    return grid, np.array(sorted(edges))


def thin_demes(coord, grid, max_per_deme, seed=20260817):
    """Cap how many samples any one deme may contribute; returns a boolean keep mask.

    Sampling here is very unbalanced (one deme holds 50 of 335). A deme's allele frequency is
    estimated from its own samples, so an over-sampled deme is estimated far more precisely
    than its neighbours and can anchor the surface around it. Thinning tests whether a feature
    is a property of the data or of that one deme.
    """
    from scipy.spatial import cKDTree

    assign = cKDTree(grid).query(coord, k=1)[1]
    rng = np.random.default_rng(seed)
    keep = np.ones(coord.shape[0], dtype=bool)
    for d in np.unique(assign):
        idx = np.where(assign == d)[0]
        if len(idx) > max_per_deme:
            keep[rng.permutation(idx)[max_per_deme:]] = False
    return keep, assign


def build_graph(geno, coord, outer, grid_res, spacing=None, scale_snps=True):
    if spacing:
        grid, edges = make_tri_grid(outer, spacing)
    else:
        _, edges, grid, _ = prepare_graph_inputs(
            coord=coord, ggrid=feems_grid_path(grid_res), translated=False, outer=outer
        )
    import networkx as nx

    g = nx.Graph()
    g.add_nodes_from(range(1, grid.shape[0] + 1))
    g.add_edges_from(map(tuple, edges))
    assert nx.is_connected(g), (
        f"deme graph is disconnected ({nx.number_connected_components(g)} components); "
        "increase --grid-spacing or the region buffer"
    )
    sp_graph = SpatialGraph(geno, coord, grid, edges, scale_snps=scale_snps)
    assert np.isfinite(sp_graph.S).all(), (
        "sample covariance S contains nan/inf -- non-polymorphic sites survived (feems #30)"
    )
    return sp_graph, grid, edges


def plot_surface(sp_graph, lamb, out_stem, title):
    fig = plt.figure(figsize=(8, 6), dpi=200)
    import cartopy.crs as ccrs

    proj = ccrs.EquidistantConic(central_longitude=float(np.mean(sp_graph.node_pos[:, 0])),
                                 central_latitude=float(np.mean(sp_graph.node_pos[:, 1])))
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    # feems defaults the colour scale to +/-2 log10 units. A well-regularised surface can
    # span far less than that and then renders almost blank, so scale to the actual spread.
    lw = np.log10(np.array(sp_graph.w, dtype=float))
    abs_max = float(np.clip(np.round(np.abs(lw - lw.mean()).max() + 0.05, 1), 0.2, 2.0))
    v = Viz(ax, sp_graph, projection=proj, edge_width=0.5, edge_alpha=1.0, edge_zorder=100,
            sample_pt_size=8, obs_node_size=5, sample_pt_color="black", cbar_font_size=8,
            cbar_ticklabelsize=8, abs_max=abs_max)
    v.draw_map()
    v.draw_edges(use_weights=True)
    v.draw_obs_nodes(use_ids=False)
    v.draw_edge_colorbar()
    ax.set_title(title, fontsize=9)
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", bbox_inches="tight")
    plt.close(fig)


def deme_distances(sp_graph):
    """Return (observed, fitted) between-deme genetic distances, lower triangle.

    feems' own `comp_mats` guards its centring step with `hasattr(sp_graph.q, 'mu')`;
    `q` is the residual-variance ndarray, which never has `.mu`, so emp_cov comes back
    UNCENTRED and is not comparable to fit_cov. Centre the deme allele frequencies here
    (the standard EEMS/FEEMS empirical covariance) and compare on the distance scale,
    which is what the feems goodness-of-fit plot uses.
    """
    from feems.objective import Objective, comp_mats
    from feems.utils import cov_to_dist

    obj = Objective(sp_graph)
    fit_cov, _, _ = comp_mats(obj)
    F = sp_graph.frequencies
    Fc = F - F.mean(axis=0)
    emp_cov = Fc @ Fc.T / sp_graph.n_snps
    tril = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
    return cov_to_dist(emp_cov)[tril], cov_to_dist(fit_cov)[tril]


def plot_fit_vs_dist(sp_graph, out_path, title):
    """Fitted vs observed genetic distance -- the standard feems goodness-of-fit check."""
    from scipy.stats import pearsonr

    x, y = deme_distances(sp_graph)
    r, _ = pearsonr(x, y)
    fig, ax = plt.subplots(figsize=(4.5, 4.5), dpi=200)
    ax.scatter(x, y, s=6, alpha=0.5, edgecolors="none")
    lim = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xlabel("observed genetic distance (deme pairs)")
    ax.set_ylabel("fitted genetic distance")
    ax.set_title(f"{title}\nPearson r = {r:.3f} (R² = {r ** 2:.3f})", fontsize=8)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return float(r)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", required=True)
    p.add_argument("--out-tag", default=None,
                   help="prefix for outputs (default: --tag); use for sensitivity runs off the same input")
    p.add_argument("--in-dir", default="data/processed/feems")
    p.add_argument("--fig-dir", default="figures/phase5")
    p.add_argument("--grid-res", default="100", choices=sorted(FEEMS_GRIDS),
                   help="bundled feems global grid; used only when --grid-spacing is 0")
    p.add_argument("--grid-spacing", type=float, default=0.30,
                   help="triangular lattice spacing in degrees (0 = use the bundled feems grid)")
    p.add_argument("--lamb-grid", default="1e-6,1e2,10",
                   help="geomspace lo,hi,n for the lambda sweep")
    p.add_argument("--max-per-deme", type=int, default=0,
                   help="cap samples per deme (0 = no cap); sampling-balance sensitivity")
    p.add_argument("--cv-folds", type=int, default=10, help="0 = leave-one-deme-out")
    p.add_argument("--no-cv", action="store_true")
    p.add_argument("--no-scale-snps", action="store_true")
    args = p.parse_args()

    os.makedirs(args.fig_dir, exist_ok=True)
    out_tag = args.out_tag or args.tag
    stem_in = os.path.join(args.in_dir, args.tag)
    stem_out = os.path.join(args.in_dir, out_tag)
    stem_fig = os.path.join(args.fig_dir, out_tag)
    t0 = time.time()

    print("=== environment ===")
    compat_msg = feems_compat.patch_cholmod()
    feems_compat.self_test()

    # ---- inputs -------------------------------------------------------------
    print("=== inputs ===")
    geno_raw, ids = load_genotypes(stem_in + "_feems")
    print(f"  genotypes: {geno_raw.shape[0]} samples x {geno_raw.shape[1]} variants")

    coord_by_id = {}
    with open(stem_in + "_coords.tsv") as fh:
        for line in fh:
            iid, lo, la = line.rstrip("\n").split("\t")
            coord_by_id[iid] = (float(lo), float(la))
    missing = [i for i in ids if i not in coord_by_id]
    assert not missing, f"{len(missing)} .fam samples have no coordinate"
    coord = np.array([coord_by_id[i] for i in ids])
    assert coord.shape == (geno_raw.shape[0], 2)
    assert np.isfinite(coord).all()
    # coord row i corresponds to .fam row i by ID lookup -- this is the positional
    # contract SpatialGraph relies on, established by join rather than assumed.
    print(f"  coords joined BY SAMPLE ID to .fam order: {coord.shape[0]} rows, "
          f"long [{coord[:, 0].min():.3f}, {coord[:, 0].max():.3f}] "
          f"lat [{coord[:, 1].min():.3f}, {coord[:, 1].max():.3f}]")

    outer = np.loadtxt(stem_in + "_outer.csv", delimiter=",", skiprows=1)
    print(f"  outer ring: {outer.shape[0]} vertices")

    geno, miss_rate, n_all_missing = mean_impute(geno_raw)
    print(f"  missingness: mean per-SNP {miss_rate.mean():.4f}, max {miss_rate.max():.4f}, "
          f"all-missing columns {n_all_missing} (mean-imputed)")
    geno = assert_polymorphic(geno, label="post-impute")
    n_polymorphic_sites = int(geno.shape[1])

    n_thinned = 0
    if args.max_per_deme:
        tmp_grid, _ = make_tri_grid(outer, args.grid_spacing) if args.grid_spacing else (None, None)
        if tmp_grid is None:
            raise SystemExit("--max-per-deme requires --grid-spacing > 0")
        keep, assign = thin_demes(coord, tmp_grid, args.max_per_deme)
        n_thinned = int((~keep).sum())
        before = np.bincount(assign).max()
        geno, coord = geno[keep], coord[keep]
        ids = [i for i, k in zip(ids, keep) if k]
        print(f"  deme thinning: cap {args.max_per_deme}/deme -> dropped {n_thinned} samples, "
              f"{geno.shape[0]} remain (largest deme was {before})")
        geno = assert_polymorphic(geno, label="post-thin")

    # ---- graph --------------------------------------------------------------
    print("=== spatial graph ===")
    sp_graph, grid, edges = build_graph(geno, coord, outer, args.grid_res,
                                        spacing=args.grid_spacing,
                                        scale_snps=not args.no_scale_snps)
    grid_desc = (f"triangular lattice, spacing {args.grid_spacing} deg"
                 if args.grid_spacing else f"feems {FEEMS_GRIDS[args.grid_res]}")
    print(f"  {grid_desc}: {grid.shape[0]} demes, {edges.shape[0]} edges, "
          f"{sp_graph.n_observed_nodes} observed demes")

    lo, hi, n = args.lamb_grid.split(",")
    lamb_grid = np.geomspace(float(lo), float(hi), int(n))[::-1]  # descending: warm-start order

    # ---- CV -----------------------------------------------------------------
    cv_summary = None
    lamb_cv = None
    if not args.no_cv:
        print(f"=== cross-validation ({args.cv_folds or sp_graph.n_observed_nodes} folds) ===")
        cv_err = run_cv(sp_graph, lamb_grid, n_folds=args.cv_folds or None,
                        factr=1e10, outer_verbose=False, inner_verbose=False)
        mean_cv = np.mean(cv_err, axis=(0, 2))
        lamb_cv = float(lamb_grid[np.argmin(mean_cv)])
        cv_summary = {"lamb": lamb_grid.tolist(), "mean_cv_err": mean_cv.tolist(),
                      "n_folds": int(cv_err.shape[0])}
        for l_, e_ in zip(lamb_grid, mean_cv):
            print(f"  lamb={l_:<12.6g} cv_err={e_:.6f}{'   <-- min' if l_ == lamb_cv else ''}")
        fig, ax = plt.subplots(figsize=(5, 3.5), dpi=200)
        ax.plot(lamb_grid, mean_cv, "o-", ms=4)
        ax.axvline(lamb_cv, color="crimson", ls="--", lw=1, label=f"chosen $\\lambda$={lamb_cv:.4g}")
        ax.set_xscale("log")
        ax.set_xlabel("$\\lambda$")
        ax.set_ylabel("mean CV error")
        ax.set_title(f"{args.tag}: FEEMS {cv_summary['n_folds']}-fold CV", fontsize=9)
        ax.legend(fontsize=7)
        fig.savefig(f"{stem_fig}_feems_cv.png", bbox_inches="tight")
        plt.close(fig)

    lamb_star = lamb_cv if lamb_cv is not None else float(np.sqrt(float(lo) * float(hi)))

    # ---- sweep of full fits -------------------------------------------------
    print("=== fitting full model across the lambda sweep ===")
    weights, sweep = {}, []
    for lamb in lamb_grid:
        lamb = float(lamb)
        try:
            sp_graph.fit(lamb=lamb, optimize_q=None, verbose=False)
            ok, err = True, None
        except AssertionError as e:  # feems raises AssertionError("did not converge")
            ok, err = False, str(e)
            print(f"  lamb={lamb:<12.6g} FAILED: {err}")
            sweep.append({"lamb": lamb, "converged": False, "error": err})
            continue
        w = np.array(sp_graph.w, dtype=float)
        weights[f"{lamb:.6g}"] = w
        tag = "chosen" if lamb == lamb_star else ""
        stem = f"{stem_fig}_feems_surface_lamb{lamb:.6g}" + ("_CHOSEN" if tag else "")
        plot_surface(sp_graph, lamb, stem,
                     f"{args.tag}: FEEMS effective migration, $\\lambda$={lamb:.4g}"
                     + (" (chosen)" if tag else ""))
        r = plot_fit_vs_dist(sp_graph, f"{stem_fig}_feems_fit_vs_dist_lamb{lamb:.6g}.png",
                             f"{args.tag}, $\\lambda$={lamb:.4g}")
        rec = {"lamb": lamb, "converged": ok, "fit_r": r,
               "log10_w": {"min": float(np.log10(w).min()), "max": float(np.log10(w).max()),
                           "sd": float(np.log10(w).std())},
               "s2": float(sp_graph.s2), "train_loss": float(sp_graph.train_loss)}
        sweep.append(rec)
        print(f"  lamb={lamb:<12.6g} r={r:.4f}  sd(log10 w)={rec['log10_w']['sd']:.4f}  "
              f"range log10 w=[{rec['log10_w']['min']:.2f}, {rec['log10_w']['max']:.2f}]")

    # ---- refit at chosen lambda and persist ---------------------------------
    sp_graph.fit(lamb=lamb_star, optimize_q=None, verbose=False)
    np.savez_compressed(stem_out + "_weights.npz",
                        node_pos=sp_graph.node_pos, edges=np.array(edges),
                        grid=grid, coord=coord, outer=outer,
                        w_chosen=np.array(sp_graph.w, dtype=float),
                        **{f"w_lamb_{k}": v for k, v in weights.items()})

    # deme-level table for the IBD/Fst consistency check
    obs_pos = sp_graph.node_pos[sp_graph.perm_idx[: sp_graph.n_observed_nodes]]
    n_per = np.asarray(sp_graph.n_samples_per_obs_node_permuted)

    out = {
        "tag": out_tag,
        "input_tag": args.tag,
        "feems_version": __import__("feems").__version__,
        "sksparse_compat": compat_msg,
        "inputs": {
            "bfile": stem_in + "_feems",
            "coords": stem_in + "_coords.tsv",
            "coords_joined_by": "sample ID (not row position)",
            "outer_ring_vertices": int(outer.shape[0]),
            "n_samples": int(geno.shape[0]),
            "max_per_deme": args.max_per_deme or None,
            "n_samples_thinned_out": n_thinned,
            "n_variants_in_bed": int(geno_raw.shape[1]),
            "n_polymorphic_sites": n_polymorphic_sites,
            "invariant_sites_dropped": True,
            "n_zero_variance_dropped_here": int(geno_raw.shape[1] - n_polymorphic_sites),
            "mean_snp_missingness": float(miss_rate.mean()),
            "max_snp_missingness": float(miss_rate.max()),
            "imputation": "per-SNP mean",
            "scale_snps": not args.no_scale_snps,
        },
        "graph": {
            "grid": grid_desc,
            "grid_spacing_deg": args.grid_spacing or None,
            "n_demes": int(grid.shape[0]),
            "n_edges": int(np.array(edges).shape[0]),
            "n_observed_demes": int(sp_graph.n_observed_nodes),
            "samples_per_observed_deme": {
                "min": int(n_per.min()), "median": float(np.median(n_per)), "max": int(n_per.max())
            },
        },
        "lamb_sweep": sweep,
        "cv": cv_summary,
        "lamb_chosen": lamb_star,
        "lamb_chosen_rule": "min mean CV error" if lamb_cv is not None else "geometric midpoint (CV skipped)",
        "runtime_sec": round(time.time() - t0, 1),
    }
    with open(stem_out + "_fit.json", "w") as fh:
        json.dump(out, fh, indent=2)

    np.savetxt(stem_out + "_demes.tsv",
               np.column_stack([obs_pos, n_per]), delimiter="\t",
               header="long\tlat\tn_samples", comments="", fmt=["%.6f", "%.6f", "%d"])

    print(f"=== done in {out['runtime_sec']}s; chosen lambda = {lamb_star:.6g} ===")
    print(f"  wrote {stem_out}_fit.json, {stem_out}_weights.npz, {stem_out}_demes.tsv")
    print(f"  figures in {args.fig_dir}/")


if __name__ == "__main__":
    main()
