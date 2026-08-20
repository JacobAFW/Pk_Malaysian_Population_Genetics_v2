#!/usr/bin/env python
"""Phase 5 / FEEMS step 4: block-bootstrap the migration surface.

Refits the surface at the CHOSEN lambda on `--n-boot` resampled SNP sets, so every edge
gets a sampling distribution instead of a point estimate. Nothing in the surface should be
called a feature unless it survives this.

Resampling is a BLOCK bootstrap over contiguous runs of SNPs within each chromosome
(`--block-size` SNPs per block, blocks drawn with replacement to the original SNP count).
Resampling individual SNPs would treat linked sites as independent and give false precision.

Per-edge statistic is the same quantity feems plots: r_e = log10(w_e) - mean(log10(w)),
i.e. log10 of the edge weight relative to the surface mean. An edge is called SUPPORTED when
its `--ci` percentile interval across bootstraps excludes zero, which is the same thing as
its sign agreeing with the point estimate in >= (1 - (1-ci)/2) of bootstraps.

Outputs:
  data/processed/feems/<tag>_bootstrap.npz    r_boot (n_boot x n_edges), r_hat, edge endpoints
  results/phase5/<tag>_bootstrap.json         support counts, CI widths, top stable edges
  figures/phase5/<tag>_feems_bootstrap.png    point estimate | supported-only  (2 panels)
  figures/phase5/<tag>_feems_bootstrap_support.png   sign-consistency + CI-width distributions

Run: micromamba run -n feems_e python scripts/2x_feems_bootstrap.py --tag Mf_only --n-boot 500
"""
from __future__ import annotations

import argparse
import io
import json
import os
import sys
import time
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def load_fit_module():
    import importlib.util as ilu

    here = os.path.dirname(os.path.abspath(__file__))
    spec = ilu.spec_from_file_location("fit_mod", os.path.join(here, "2x_feems_fit.py"))
    mod = ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.feems_compat.patch_cholmod(verbose=False)
    return mod


def read_bim_chrom(bim_path):
    """Chromosome code per variant, in .bim order."""
    chrom = []
    with open(bim_path) as fh:
        for line in fh:
            f = line.split()
            if f:
                chrom.append(f[0])
    return np.array(chrom)


def make_blocks(chrom, block_size):
    """Contiguous blocks of `block_size` SNPs, never spanning a chromosome boundary."""
    blocks = []
    for c in np.unique(chrom):
        idx = np.where(chrom == c)[0]
        for s in range(0, len(idx), block_size):
            blocks.append(idx[s:s + block_size])
    return blocks


def quiet_graph(fit, geno, coord, grid, edges, scale_snps):
    """SpatialGraph construction, with feems' progress chatter suppressed."""
    from feems import SpatialGraph

    saved = sys.stdout
    sys.stdout = io.StringIO()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return SpatialGraph(geno, coord, grid, edges, scale_snps=scale_snps)
    finally:
        sys.stdout = saved


def rel_log_weights(w):
    """r_e = log10(w_e) - mean(log10 w): the quantity feems colours the map by."""
    lw = np.log10(np.asarray(w, dtype=float))
    return lw - lw.mean()


def node_perm_positions(sp_graph, n_nodes):
    """grid node index (0-based) -> its slot in feems' permuted node order."""
    perm = np.asarray(sp_graph.perm_idx)
    pos = np.full(n_nodes, -1, dtype=int)
    pos[perm] = np.arange(len(perm))
    assert (pos >= 0).all()
    return pos


def w_by_edge(sp_graph, edges, n_nodes):
    """Fitted weights re-indexed onto the CANONICAL `edges` order.

    `sp_graph.w` is ordered over the permuted upper-triangular adjacency, and that permutation
    depends on WHICH demes are observed. Under an individual bootstrap the observed set
    changes between replicates, so raw `w` vectors are not comparable. Reading the weights off
    the permuted adjacency and back into `edges` order fixes them in a common frame.
    """
    pos = node_perm_positions(sp_graph, n_nodes)
    W = sp_graph.inv_triu(np.asarray(sp_graph.w, dtype=float), perm=True).toarray()  # already symmetric
    a = pos[np.asarray(edges)[:, 0] - 1]
    b = pos[np.asarray(edges)[:, 1] - 1]
    w = W[a, b]
    assert (w > 0).all(), "edge missing from the fitted adjacency"
    return w


def edge_to_w_index(sp_graph, edges, n_nodes):
    """canonical `edges` row -> index into sp_graph.w (needed to draw a canonical-order vector)."""
    pos = node_perm_positions(sp_graph, n_nodes)
    n_edges = np.asarray(edges).shape[0]
    M = sp_graph.inv_triu(np.arange(1, n_edges + 1, dtype=float), perm=True).toarray()  # symmetric
    a = pos[np.asarray(edges)[:, 0] - 1]
    b = pos[np.asarray(edges)[:, 1] - 1]
    idx = M[a, b].astype(int) - 1
    assert (idx >= 0).all() and len(np.unique(idx)) == n_edges, "edge <-> w index map is not a bijection"
    return idx


def drop_invariant(geno):
    """Resampling individuals can make a column invariant; feems 1.0.0 dies on that
    (its own guard raises TypeError before removing anything). Strip them here."""
    f = geno.mean(axis=0) / 2.0
    keep = (geno.var(axis=0) > 0) & (f > 0) & (f < 1)
    return geno[:, keep], int((~keep).sum())


def draw_surface_on_ax(fig, gs_pos, sp_graph, w_draw, title, abs_max):
    import cartopy.crs as ccrs
    from feems import Viz

    proj = ccrs.EquidistantConic(central_longitude=float(np.mean(sp_graph.node_pos[:, 0])),
                                 central_latitude=float(np.mean(sp_graph.node_pos[:, 1])))
    ax = fig.add_subplot(*gs_pos, projection=proj)
    saved_w = np.array(sp_graph.w, dtype=float)
    sp_graph.w = w_draw
    try:
        v = Viz(ax, sp_graph, projection=proj, edge_width=0.6, edge_alpha=1.0, edge_zorder=100,
                sample_pt_size=6, obs_node_size=4, sample_pt_color="black",
                cbar_font_size=7, cbar_ticklabelsize=7, abs_max=abs_max)
        v.draw_map()
        v.draw_edges(use_weights=True)
        v.draw_obs_nodes(use_ids=False)
        v.draw_edge_colorbar()
    finally:
        sp_graph.w = saved_w
    ax.set_title(title, fontsize=8)
    return ax


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", required=True)
    p.add_argument("--in-dir", default="data/processed/feems")
    p.add_argument("--fig-dir", default="figures/phase5")
    p.add_argument("--out-dir", default="results/phase5")
    p.add_argument("--n-boot", type=int, default=500)
    p.add_argument("--resample", choices=["snp", "individual"], default="snp",
                   help="snp: block bootstrap over loci (uncertainty in WHICH SITES we typed). "
                        "individual: bootstrap over samples (uncertainty in WHICH INFECTIONS "
                        "we sampled) -- normally the binding constraint.")
    p.add_argument("--out-tag", default=None, help="output prefix (default: <tag>_<resample>)")
    p.add_argument("--block-size", type=int, default=500, help="SNPs per bootstrap block")
    p.add_argument("--ci", type=float, default=0.95)
    p.add_argument("--seed", type=int, default=20260817)
    p.add_argument("--lamb", type=float, default=None, help="default: chosen lambda from <tag>_fit.json")
    p.add_argument("--grid-spacing", type=float, default=None, help="default: from <tag>_fit.json")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.fig_dir, exist_ok=True)
    stem_in = os.path.join(args.in_dir, args.tag)
    out_tag = args.out_tag or f"{args.tag}_{args.resample}"
    fit = load_fit_module()

    with open(stem_in + "_fit.json") as fh:
        meta = json.load(fh)
    lamb = args.lamb if args.lamb is not None else float(meta["lamb_chosen"])
    spacing = args.grid_spacing if args.grid_spacing is not None else float(meta["graph"]["grid_spacing_deg"])
    scale_snps = bool(meta["inputs"]["scale_snps"])
    print(f"=== bootstrap {args.tag} [{args.resample}]: lambda={lamb:g}, spacing={spacing}, "
          f"scale_snps={scale_snps}, n_boot={args.n_boot} -> tag {out_tag} ===")

    geno_raw, ids = fit.load_genotypes(stem_in + "_feems")
    geno, _, _ = fit.mean_impute(geno_raw)
    geno = fit.assert_polymorphic(geno, label="bootstrap")
    coord_by_id = {}
    with open(stem_in + "_coords.tsv") as fh:
        for line in fh:
            iid, lo, la = line.rstrip("\n").split("\t")
            coord_by_id[iid] = (float(lo), float(la))
    coord = np.array([coord_by_id[i] for i in ids])
    outer = np.loadtxt(stem_in + "_outer.csv", delimiter=",", skiprows=1)
    grid, edges = fit.make_tri_grid(outer, spacing)

    chrom = read_bim_chrom(stem_in + "_feems.bim")
    assert len(chrom) == geno_raw.shape[1], "bim/genotype length mismatch"
    if geno.shape[1] != geno_raw.shape[1]:
        raise SystemExit("ERROR: columns were dropped after the bim was read; re-run prep")
    blocks = make_blocks(chrom, args.block_size)
    print(f"  {geno.shape[1]} SNPs on {len(np.unique(chrom))} chromosomes -> {len(blocks)} blocks "
          f"of <= {args.block_size} SNPs")

    # ---- point estimate ------------------------------------------------------
    n_nodes = grid.shape[0]
    n_edges = np.asarray(edges).shape[0]
    sp_graph = quiet_graph(fit, geno, coord, grid, edges, scale_snps)
    sp_graph.fit(lamb=lamb, optimize_q=None, verbose=False)
    w_hat = w_by_edge(sp_graph, edges, n_nodes)
    r_hat = rel_log_weights(w_hat)
    print(f"  point estimate: {n_edges} edges, {sp_graph.n_observed_nodes} observed demes, "
          f"r_hat in [{r_hat.min():+.3f}, {r_hat.max():+.3f}]")

    # ---- bootstrap -----------------------------------------------------------
    rng = np.random.default_rng(args.seed)
    n_snps = geno.shape[1]
    n_samples = geno.shape[0]
    r_boot = np.full((args.n_boot, n_edges), np.nan)
    n_obs_boot = np.zeros(args.n_boot, dtype=int)
    n_dropped = []
    n_failed = 0
    t0 = time.time()
    for b in range(args.n_boot):
        if args.resample == "snp":
            pick = rng.integers(0, len(blocks), len(blocks))
            cols = np.concatenate([blocks[k] for k in pick])[:n_snps]
            g_b, c_b = geno[:, cols], coord
        else:
            rows = rng.integers(0, n_samples, n_samples)
            g_b, dropped = drop_invariant(geno[rows])
            c_b = coord[rows]
            n_dropped.append(dropped)
        try:
            sg = quiet_graph(fit, g_b, c_b, grid, edges, scale_snps)
            if not np.isfinite(sg.S).all():
                raise ValueError("non-finite covariance")
            sg.fit(lamb=lamb, optimize_q=None, verbose=False)
            r_boot[b] = rel_log_weights(w_by_edge(sg, edges, n_nodes))
            n_obs_boot[b] = sg.n_observed_nodes
        except Exception as e:  # noqa: BLE001 - a failed replicate is data, not a crash
            n_failed += 1
            if n_failed <= 3:
                print(f"  replicate {b} failed: {type(e).__name__}: {e}")
        if (b + 1) % 100 == 0:
            print(f"  {b + 1}/{args.n_boot} replicates ({time.time() - t0:.0f}s)")
    ok = np.isfinite(r_boot).all(axis=1)
    r_boot = r_boot[ok]
    n_obs_boot = n_obs_boot[ok]
    print(f"  {r_boot.shape[0]}/{args.n_boot} replicates usable ({n_failed} failed)")
    if args.resample == "individual":
        print(f"  observed demes per replicate: min {n_obs_boot.min()}, "
              f"median {np.median(n_obs_boot):.0f}, max {n_obs_boot.max()} "
              f"(point estimate {sp_graph.n_observed_nodes})")
        print(f"  invariant columns dropped per replicate: median {int(np.median(n_dropped))}")
    assert r_boot.shape[0] >= 20, "too few usable bootstrap replicates to summarise"

    # ---- per-edge support ----------------------------------------------------
    a = (1 - args.ci) / 2
    lo = np.percentile(r_boot, 100 * a, axis=0)
    hi = np.percentile(r_boot, 100 * (1 - a), axis=0)
    supported = (lo > 0) | (hi < 0)
    sign_consistency = np.where(
        r_hat >= 0, (r_boot > 0).mean(axis=0), (r_boot < 0).mean(axis=0))
    ci_width = hi - lo

    print(f"=== support at {args.ci:.0%} CI ===")
    print(f"  edges with CI excluding 0 : {supported.sum()} / {n_edges} "
          f"({supported.mean():.1%})")
    print(f"    of which reduced migration (r<0): {int((supported & (r_hat < 0)).sum())}, "
          f"elevated (r>0): {int((supported & (r_hat > 0)).sum())}")
    print(f"  median CI width           : {np.median(ci_width):.3f} log10 units")
    print(f"  median |r_hat|            : {np.median(np.abs(r_hat)):.3f} log10 units")
    print(f"  CI width / |r_hat| ratio  : {np.median(ci_width) / max(np.median(np.abs(r_hat)), 1e-12):.1f}x")

    # edge endpoint coordinates, for naming where the stable features are
    e_arr = np.asarray(edges)
    ends = np.column_stack([grid[e_arr[:, 0] - 1], grid[e_arr[:, 1] - 1]])
    w_index = edge_to_w_index(sp_graph, edges, n_nodes)

    order = np.argsort(r_hat)
    def describe(k, label):
        rows = []
        for e in k:
            rows.append({
                "r_hat": float(r_hat[e]), "ci_lo": float(lo[e]), "ci_hi": float(hi[e]),
                "sign_consistency": float(sign_consistency[e]), "supported": bool(supported[e]),
                "midpoint_long": float((ends[e, 0] + ends[e, 2]) / 2),
                "midpoint_lat": float((ends[e, 1] + ends[e, 3]) / 2),
            })
        print(f"  {label}:")
        for r in rows:
            print(f"    r={r['r_hat']:+.3f} CI[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}] "
                  f"sign={r['sign_consistency']:.2f} {'SUPPORTED' if r['supported'] else 'not supported'} "
                  f"@ {r['midpoint_long']:.2f}E {r['midpoint_lat']:.2f}N")
        return rows

    print("=== the extremes of the point-estimate surface ===")
    lowest = describe(order[:8], "8 lowest-weight edges (candidate barriers)")
    highest = describe(order[::-1][:8], "8 highest-weight edges (candidate corridors)")

    # ---- figures -------------------------------------------------------------
    abs_max = float(np.clip(np.round(np.abs(r_hat).max() + 0.05, 1), 0.2, 2.0))
    # back into sp_graph.w order for drawing
    w_draw = np.empty(n_edges)
    w_draw[w_index] = w_hat
    w_mask = np.empty(n_edges)
    w_mask[w_index] = np.where(supported, w_hat, 10 ** np.log10(w_hat).mean())
    scheme = (f"block bootstrap, {args.block_size}-SNP blocks" if args.resample == "snp"
              else "bootstrap over individuals")
    fig = plt.figure(figsize=(13, 5), dpi=200)
    draw_surface_on_ax(fig, (1, 2, 1), sp_graph, w_draw,
                       f"{args.tag}: point estimate, $\\lambda$={lamb:g}", abs_max)
    draw_surface_on_ax(fig, (1, 2, 2), sp_graph, w_mask,
                       f"bootstrap-supported edges only ({supported.sum()}/{n_edges} at "
                       f"{args.ci:.0%} CI, {r_boot.shape[0]} replicates)", abs_max)
    fig.suptitle(f"{args.tag}: FEEMS effective migration — {scheme}", fontsize=10)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(args.fig_dir, f"{out_tag}_feems_bootstrap.{ext}"),
                    bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), dpi=200)
    axes[0].hist(sign_consistency, bins=np.linspace(0.5, 1.0, 26), color="tab:blue")
    axes[0].axvline(1 - a, color="crimson", ls="--", lw=1, label=f"{args.ci:.0%} CI threshold")
    axes[0].set_xlabel("bootstrap sign consistency")
    axes[0].set_ylabel("edges")
    axes[0].legend(fontsize=7)
    axes[0].set_title("how often each edge keeps its sign", fontsize=9)

    axes[1].scatter(np.abs(r_hat), ci_width, s=8, alpha=0.6, edgecolors="none",
                    c=np.where(supported, "tab:green", "tab:grey"))
    lim = [0, max(np.abs(r_hat).max(), ci_width.max()) * 1.05]
    axes[1].plot(lim, lim, "k--", lw=0.8, label="CI width = |effect|")
    axes[1].set_xlabel("|r| point estimate (log10 units)")
    axes[1].set_ylabel(f"{args.ci:.0%} CI width")
    axes[1].legend(fontsize=7)
    axes[1].set_title("effect size vs uncertainty\n(green = CI excludes 0)", fontsize=9)

    # vlines, not errorbar: under a strong bootstrap the point estimate can sit outside the
    # percentile interval (bootstrap bias), which errorbar rejects as a negative arm.
    axes[2].vlines(np.arange(n_edges), lo[order], hi[order], color="lightgrey", lw=0.6)
    axes[2].scatter(np.arange(n_edges), r_hat[order], s=5,
                    c=np.where(supported[order], "tab:green", "tab:grey"))
    axes[2].axhline(0, color="k", lw=0.8)
    axes[2].set_xlabel("edges, ranked by fitted weight")
    axes[2].set_ylabel("r = log10(w / mean w)")
    axes[2].set_title("per-edge estimate with bootstrap CI", fontsize=9)
    fig.tight_layout()
    fig.suptitle(f"{args.tag} — {scheme}", fontsize=10)
    fig.savefig(os.path.join(args.fig_dir, f"{out_tag}_feems_bootstrap_support.png"),
                bbox_inches="tight")
    plt.close(fig)

    # ---- persist -------------------------------------------------------------
    np.savez_compressed(os.path.join(args.in_dir, out_tag) + "_bootstrap.npz", r_boot=r_boot, r_hat=r_hat, w_hat=w_hat,
                        ci_lo=lo, ci_hi=hi, supported=supported,
                        sign_consistency=sign_consistency, edge_endpoints=ends)
    out = {
        "tag": out_tag, "input_tag": args.tag, "resample": args.resample, "lamb": lamb, "grid_spacing_deg": spacing, "scale_snps": scale_snps,
        "n_boot_requested": args.n_boot, "n_boot_usable": int(r_boot.shape[0]),
        "block_size_snps": args.block_size, "n_blocks": len(blocks), "n_snps": int(n_snps),
        "seed": args.seed, "ci": args.ci,
        "n_edges": int(n_edges),
        "n_supported": int(supported.sum()),
        "frac_supported": float(supported.mean()),
        "n_supported_reduced": int((supported & (r_hat < 0)).sum()),
        "n_supported_elevated": int((supported & (r_hat > 0)).sum()),
        "median_ci_width": float(np.median(ci_width)),
        "median_abs_r_hat": float(np.median(np.abs(r_hat))),
        "r_hat_range": [float(r_hat.min()), float(r_hat.max())],
        "lowest_weight_edges": lowest,
        "highest_weight_edges": highest,
    }
    if args.resample == "individual":
        out["observed_demes_per_replicate"] = {
            "min": int(n_obs_boot.min()), "median": float(np.median(n_obs_boot)),
            "max": int(n_obs_boot.max()), "point_estimate": int(sp_graph.n_observed_nodes)}
        out["median_invariant_cols_dropped"] = int(np.median(n_dropped))
    with open(os.path.join(args.out_dir, f"{out_tag}_bootstrap.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"=== wrote {args.out_dir}/{out_tag}_bootstrap.json, "
          f"{args.in_dir}/{out_tag}_bootstrap.npz, and 2 figures in {args.fig_dir}/ ===")


if __name__ == "__main__":
    main()
