#!/usr/bin/env python
"""Phase 5 / FEEMS step 3: does the fitted migration surface agree with our IBD and Fst?

The surface is only credible if it reproduces structure we already measured by other means:
  high effective migration  <->  high IBD sharing / low differentiation
  barriers (low w)          <->  low IBD sharing / high Fst

Three checks, all on the SAME 335 samples the surface was fitted to:

  1. IBD vs FEEMS resistance. hmmIBD `fract_sites_IBD` for every sample pair, against the
     FEEMS fitted between-deme genetic distance and against great-circle distance.
     Expect IBD to fall with both; the surface adds value only if it beats raw geography.
  2. Barrier vs cluster boundary. Global Mf-vs-Mn Fst is high, so edges spanning the
     Mf/Mn boundary should carry LOW weight. Per edge, compare fitted log10(w) against
     the difference in Mn-ancestry fraction between its two endpoint demes.
  3. Deme-level IBD vs deme-level fitted distance, which removes within-deme pairs and
     is the resolution the model actually works at.

Outputs: results/phase5/<tag>_consistency.json, figures/phase5/<tag>_feems_consistency.png

Run: micromamba run -n feems_e python scripts/2x_feems_consistency.py --tag cleaned_501
"""
from __future__ import annotations

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

EARTH_R_KM = 6371.0088


def haversine(a, b):
    """Great-circle km between two (long, lat) arrays, broadcast over rows."""
    lo1, la1 = np.radians(a[:, 0]), np.radians(a[:, 1])
    lo2, la2 = np.radians(b[:, 0]), np.radians(b[:, 1])
    d = (np.sin((la2 - la1) / 2) ** 2
         + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2)
    return 2 * EARTH_R_KM * np.arcsin(np.sqrt(np.clip(d, 0, 1)))


def rebuild_graph(tag, in_dir, spacing, grid_res, scale_snps=True):
    """Reconstruct the fitted SpatialGraph so we can read fitted distances and weights."""
    import importlib.util as ilu

    here = os.path.dirname(os.path.abspath(__file__))
    spec = ilu.spec_from_file_location("fit_mod", os.path.join(here, "2x_feems_fit.py"))
    fit = ilu.module_from_spec(spec)
    spec.loader.exec_module(fit)
    fit.feems_compat.patch_cholmod(verbose=False)

    stem = os.path.join(in_dir, tag)
    geno, ids = fit.load_genotypes(stem + "_feems")
    geno, _, _ = fit.mean_impute(geno)
    geno = fit.assert_polymorphic(geno, label="consistency")
    coord_by_id = {}
    with open(stem + "_coords.tsv") as fh:
        for line in fh:
            iid, lo, la = line.rstrip("\n").split("\t")
            coord_by_id[iid] = (float(lo), float(la))
    coord = np.array([coord_by_id[i] for i in ids])
    outer = np.loadtxt(stem + "_outer.csv", delimiter=",", skiprows=1)
    sp_graph, grid, edges = fit.build_graph(geno, coord, outer, grid_res,
                                            spacing=spacing, scale_snps=scale_snps)
    return fit, sp_graph, grid, edges, coord, ids


def sample_to_deme(sp_graph, n_samples):
    """Map sample index -> observed-deme index (0..n_observed-1), matching feems' permutation."""
    import networkx as nx

    from feems.spatial_graph import query_node_attributes

    sample_idx = nx.get_node_attributes(sp_graph, "sample_idx")
    permuted_idx = query_node_attributes(sp_graph, "permuted_idx")
    out = np.full(n_samples, -1, dtype=int)
    for i, node_id in enumerate(permuted_idx[: sp_graph.n_observed_nodes]):
        for s in sample_idx[node_id]:
            out[s] = i
    assert (out >= 0).all(), "some samples were not assigned to an observed deme"
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", default="cleaned_501")
    p.add_argument("--in-dir", default="data/processed/feems")
    p.add_argument("--fig-dir", default="figures/phase5")
    p.add_argument("--out-dir", default="results/phase5")
    p.add_argument("--ibd", default="data/gadi/validation/hmmIBD/Pk.hmm_fract.txt",
                   help="hmmIBD .hmm_fract table (sample1 sample2 ... fract_sites_IBD)")
    p.add_argument("--clusters", default="data/processed/cluster_maf/cluster_labels.tsv")
    p.add_argument("--fst-mean", default="data/gadi/validation/fst/Fst_mean_mf_vs_mn.txt")
    p.add_argument("--grid-spacing", type=float, default=0.30)
    p.add_argument("--grid-res", default="100")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.fig_dir, exist_ok=True)
    stem_in = os.path.join(args.in_dir, args.tag)

    with open(stem_in + "_fit.json") as fh:
        fitmeta = json.load(fh)
    lamb = fitmeta["lamb_chosen"]
    print(f"=== rebuilding the fitted surface at the chosen lambda = {lamb:g} ===")

    fit, sp_graph, grid, edges, coord, ids = rebuild_graph(
        args.tag, args.in_dir, args.grid_spacing, args.grid_res,
        scale_snps=fitmeta["inputs"]["scale_snps"])
    sp_graph.fit(lamb=lamb, optimize_q=None, verbose=False)
    w = np.array(sp_graph.w, dtype=float)

    emp_d, fit_d = fit.deme_distances(sp_graph)
    n_obs = sp_graph.n_observed_nodes
    tril = np.tril_indices(n_obs, k=-1)
    fit_D = np.zeros((n_obs, n_obs))
    fit_D[tril] = fit_d
    fit_D += fit_D.T
    emp_D = np.zeros((n_obs, n_obs))
    emp_D[tril] = emp_d
    emp_D += emp_D.T

    s2d = sample_to_deme(sp_graph, len(ids))
    idx_of = {iid: i for i, iid in enumerate(ids)}
    deme_pos = sp_graph.node_pos[sp_graph.perm_idx[:n_obs]]
    geo_D = np.array([[haversine(deme_pos[i:i + 1], deme_pos[j:j + 1])[0] for j in range(n_obs)]
                      for i in range(n_obs)])

    results = {"tag": args.tag, "lamb": lamb, "n_samples": len(ids), "n_observed_demes": int(n_obs)}

    # ---- check 1: sample-pair IBD vs FEEMS fitted distance / geography ------
    print("=== check 1: hmmIBD vs the fitted surface (sample pairs) ===")
    pairs = []
    with open(args.ibd) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c1, c2 = hdr.index("sample1"), hdr.index("sample2")
        cf = hdr.index("fract_sites_IBD")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b = idx_of.get(f[c1]), idx_of.get(f[c2])
            if a is None or b is None:
                continue
            pairs.append((a, b, float(f[cf])))
    print(f"  hmmIBD pairs with both samples in the FEEMS set: {len(pairs)}")
    if not pairs:
        raise SystemExit(f"ERROR: no usable hmmIBD pairs in {args.ibd}")

    pa = np.array([p[0] for p in pairs])
    pb = np.array([p[1] for p in pairs])
    ibd = np.array([p[2] for p in pairs])
    da, db = s2d[pa], s2d[pb]
    between = da != db
    fd = fit_D[da, db]
    gd = geo_D[da, db]

    def corr(x, y, label):
        r, pr = pearsonr(x, y)
        rho, ps = spearmanr(x, y)
        print(f"  {label:<46s} pearson r={r:+.4f} (p={pr:.2e})  spearman rho={rho:+.4f}")
        return {"pearson_r": float(r), "pearson_p": float(pr),
                "spearman_rho": float(rho), "spearman_p": float(ps), "n": int(len(x))}

    results["check1_ibd_vs_surface"] = {
        "n_pairs_total": int(len(pairs)),
        "n_pairs_between_demes": int(between.sum()),
        "ibd_vs_feems_fitted_distance": corr(fd[between], ibd[between],
                                             "IBD ~ FEEMS fitted distance (between-deme)"),
        "ibd_vs_geographic_distance": corr(gd[between], ibd[between],
                                           "IBD ~ great-circle distance (between-deme)"),
        "mean_ibd_within_deme": float(ibd[~between].mean()) if (~between).any() else None,
        "mean_ibd_between_deme": float(ibd[between].mean()),
    }
    wi = results["check1_ibd_vs_surface"]["mean_ibd_within_deme"]
    print(f"  mean IBD within-deme  = {wi:.4f}" if wi is not None else "  no within-deme pairs")
    print(f"  mean IBD between-deme = {results['check1_ibd_vs_surface']['mean_ibd_between_deme']:.4f}")

    # ---- check 2: deme-level IBD vs fitted / empirical distance -------------
    # Deme-mean IBD from a handful of pairs is very noisy, and singleton demes ALSO get an
    # inflated fitted distance (feems' fitted covariance carries a per-deme residual-variance
    # term that grows as deme sample size falls). Those two artefacts correlate positively and
    # can flip the sign of an unfiltered Pearson r. So: require a minimum number of observed
    # sample pairs per deme pair, lead with Spearman, and report a clonality-based summary
    # (fraction of pairs above an IBD threshold) alongside the mean.
    MIN_PAIRS = 5
    IBD_HIGH = 0.25
    print(f"=== check 2: deme-level IBD vs fitted distance (>={MIN_PAIRS} sample pairs per deme pair) ===")
    acc = np.zeros((n_obs, n_obs))
    hi = np.zeros((n_obs, n_obs))
    cnt = np.zeros((n_obs, n_obs))
    np.add.at(acc, (da, db), ibd)
    np.add.at(hi, (da, db), (ibd > IBD_HIGH).astype(float))
    np.add.at(cnt, (da, db), 1)
    for m in (acc, hi, cnt):
        m += m.T
    off = np.arange(n_obs)[:, None] != np.arange(n_obs)[None, :]
    ii, jj = np.where(np.tril((cnt >= MIN_PAIRS) & off, k=-1))
    deme_ibd = acc[ii, jj] / cnt[ii, jj]
    deme_hi = hi[ii, jj] / cnt[ii, jj]
    n_any = int(np.tril((cnt > 0) & off, k=-1).sum())
    print(f"  deme pairs kept: {len(ii)} (of {n_any} with any data, "
          f"{n_obs * (n_obs - 1) // 2} possible)")

    def partial_spearman(x, y, z):
        """Spearman of x vs y after regressing both on z (all rank-transformed)."""
        from scipy.stats import rankdata

        rx, ry, rz = (rankdata(v) for v in (x, y, z))
        Z = np.column_stack([np.ones_like(rz), rz])
        ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
        ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        return float(pearsonr(ex, ey)[0])

    pr = partial_spearman(fit_D[ii, jj], deme_ibd, geo_D[ii, jj])
    print(f"  {'partial spearman IBD ~ fitted dist | geography':<46s} rho={pr:+.4f}")
    results["check2_deme_level"] = {
        "min_sample_pairs_per_deme_pair": MIN_PAIRS,
        "ibd_high_threshold": IBD_HIGH,
        "n_deme_pairs_kept": int(len(ii)),
        "n_deme_pairs_with_any_data": n_any,
        "demeIBD_vs_feems_fitted_distance": corr(fit_D[ii, jj], deme_ibd,
                                                 "deme-mean IBD ~ FEEMS fitted distance"),
        "demeIBD_vs_geographic_distance": corr(geo_D[ii, jj], deme_ibd,
                                               "deme-mean IBD ~ great-circle distance"),
        "frac_highIBD_vs_feems_fitted_distance": corr(fit_D[ii, jj], deme_hi,
                                                      f"frac pairs IBD>{IBD_HIGH} ~ fitted distance"),
        "frac_highIBD_vs_geographic_distance": corr(geo_D[ii, jj], deme_hi,
                                                    f"frac pairs IBD>{IBD_HIGH} ~ geography"),
        "partial_spearman_ibd_fitted_given_geo": pr,
        "empirical_vs_fitted_distance": corr(emp_D[ii, jj], fit_D[ii, jj],
                                             "observed ~ fitted genetic distance"),
    }

    # ---- check 3: barriers vs the Mf/Mn cluster boundary --------------------
    print("=== check 3: barrier edges vs the Mf/Mn boundary (Fst = "
          f"{open(args.fst_mean).read().strip().splitlines()[-1]}) ===")
    lab = {}
    with open(args.clusters) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                lab[f[0]] = f[1]
    have_lab = [i for i, s in enumerate(ids) if s in lab]
    print(f"  cluster labels available for {len(have_lab)}/{len(ids)} samples: "
          + str({c: sum(1 for i in have_lab if lab[ids[i]] == c) for c in sorted(set(lab.values()))}))

    is_mn = np.array([1.0 if lab.get(ids[i]) == "Mn" else 0.0 for i in have_lab])
    demes_l = s2d[np.array(have_lab)]
    frac_mn = np.full(n_obs, np.nan)
    for d in range(n_obs):
        m = demes_l == d
        if m.sum():
            frac_mn[d] = is_mn[m].mean()

    # edges are 1-indexed over the FULL grid; map to observed-deme indices via perm_idx
    perm = np.asarray(sp_graph.perm_idx)
    full_to_obs = np.full(len(perm), -1)
    for obs_i, node in enumerate(perm[:n_obs]):
        full_to_obs[node] = obs_i
    ea = full_to_obs[np.asarray(edges)[:, 0] - 1]
    eb = full_to_obs[np.asarray(edges)[:, 1] - 1]
    # sp_graph.w is ordered over the permuted upper-triangular nonzeros, not `edges`;
    # read the weights straight off the permuted adjacency instead.
    W = sp_graph.inv_triu(w, perm=True).toarray()
    ok = (ea >= 0) & (eb >= 0) & np.isfinite(frac_mn[np.clip(ea, 0, None)]) \
         & np.isfinite(frac_mn[np.clip(eb, 0, None)])
    ea_o, eb_o = ea[ok], eb[ok]
    w_edge = W[ea_o, eb_o]
    keep = w_edge > 0
    ea_o, eb_o, w_edge = ea_o[keep], eb_o[keep], w_edge[keep]
    dmn = np.abs(frac_mn[ea_o] - frac_mn[eb_o])
    print(f"  edges joining two observed, cluster-labelled demes: {len(w_edge)}")
    if not np.any(dmn > 0):
        # e.g. a single-cluster run: there is no Mf/Mn boundary to test against.
        print("  SKIPPED: no Mf/Mn variation among the sampled demes in this run")
        results["check3_barrier_vs_cluster"] = {
            "skipped": "no Mf/Mn variation among sampled demes (single-cluster run)",
            "n_edges": int(len(w_edge)),
            "clusters_present": sorted({lab[ids[i]] for i in have_lab}),
        }
        dmn_plot, w_plot = None, None
    else:
        dmn_plot, w_plot = dmn, w_edge
        results["check3_barrier_vs_cluster"] = {
            "global_fst_mf_vs_mn": float(open(args.fst_mean).read().strip().splitlines()[-1]),
            "n_edges": int(len(w_edge)),
            "log10w_vs_delta_mn_fraction": corr(dmn, np.log10(w_edge),
                                                "log10(edge weight) ~ |delta Mn fraction|"),
            "median_log10w_boundary_edges":
                float(np.median(np.log10(w_edge[dmn > 0.5]))) if (dmn > 0.5).any() else None,
            "median_log10w_interior_edges":
                float(np.median(np.log10(w_edge[dmn <= 0.5]))) if (dmn <= 0.5).any() else None,
        }
        b = results["check3_barrier_vs_cluster"]
        if b["median_log10w_boundary_edges"] is not None:
            print(f"  median log10(w): boundary edges (|dMn|>0.5) = {b['median_log10w_boundary_edges']:+.3f}  "
                  f"vs interior = {b['median_log10w_interior_edges']:+.3f}")

    # ---- figure -------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), dpi=200)
    axes[0].scatter(fit_D[ii, jj], deme_ibd, s=8, alpha=0.5, edgecolors="none")
    axes[0].set_xlabel("FEEMS fitted genetic distance (deme pair)")
    axes[0].set_ylabel("mean hmmIBD fract_sites_IBD")
    c0 = results["check2_deme_level"]["demeIBD_vs_feems_fitted_distance"]
    axes[0].set_title(f"IBD vs fitted distance\nspearman $\\rho$ = {c0['spearman_rho']:+.3f} "
                      f"(pearson r = {c0['pearson_r']:+.3f})", fontsize=9)

    axes[1].scatter(geo_D[ii, jj], deme_ibd, s=8, alpha=0.5, edgecolors="none", color="tab:orange")
    axes[1].set_xlabel("great-circle distance (km)")
    axes[1].set_ylabel("mean hmmIBD fract_sites_IBD")
    c1 = results["check2_deme_level"]["demeIBD_vs_geographic_distance"]
    axes[1].set_title(f"IBD vs geography\nspearman $\\rho$ = {c1['spearman_rho']:+.3f} "
                      f"(pearson r = {c1['pearson_r']:+.3f})", fontsize=9)

    if dmn_plot is None:
        axes[2].text(0.5, 0.5, "check 3 not applicable\n(single-cluster run:\nno Mf/Mn boundary)",
                     ha="center", va="center", fontsize=9, transform=axes[2].transAxes)
        axes[2].set_axis_off()
    else:
        axes[2].scatter(dmn_plot, np.log10(w_plot), s=10, alpha=0.6, edgecolors="none",
                        color="tab:green")
        axes[2].set_xlabel("|difference in Mn fraction| across edge")
        axes[2].set_ylabel("log10 fitted edge weight $w$")
        c2 = results["check3_barrier_vs_cluster"]["log10w_vs_delta_mn_fraction"]
        axes[2].set_title(f"barriers vs Mf/Mn boundary\nspearman $\\rho$ = {c2['spearman_rho']:+.3f} "
                          f"(pearson r = {c2['pearson_r']:+.3f})", fontsize=9)
    fig.suptitle(f"{args.tag}: FEEMS surface consistency ($\\lambda$={lamb:g})", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(args.fig_dir, f"{args.tag}_feems_consistency.png"), bbox_inches="tight")
    plt.close(fig)

    out = os.path.join(args.out_dir, f"{args.tag}_consistency.json")
    with open(out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"=== wrote {out} and {args.fig_dir}/{args.tag}_feems_consistency.png ===")


if __name__ == "__main__":
    main()
