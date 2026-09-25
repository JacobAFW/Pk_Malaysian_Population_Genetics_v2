#!/usr/bin/env python3
"""score_rules.py — build the rule scoreboard from the phase-6 aggregations.

Three measure families, per rule and per floor policy (see score_rules.sh for
the policy definitions):

  1. coordinate window reproduction vs the poster's 217 truth windows
     (CHROM + window START; never WINDOW id, never BIN, integer keys only)
  2. cluster-size dependence — median windows per sample per cluster,
     zero-filled over every cluster member; reported as Spearman rho AND
     max/min ratio, because with K=3 rho takes only four values
  3. perturbation stability — window-set Jaccard of drop10 / drop20 against
     that rule's own full-scenario set, under the same floor policy

Plus the filter-4 decomposition the phase-5 review requires: of the windows
whose status differs from `absolute`, how many moved because the multi-cluster
mask moved rather than because detection changed.

Deliberately written in Python: the coordinate keying bug in phase 5 came from
R's scientific-notation stringification of doubles, and every load-bearing
number in this run is now computed in a second language.
"""
import csv, os, re, sys, statistics as st
from itertools import combinations

RES = sys.argv[1] if len(sys.argv) > 1 else "results/2026-09-25_introgression-module"
ARM = sys.argv[2] if len(sys.argv) > 2 else "full"
RULES = ["absolute", "relative", "distance", "distance_adaptive"]
SCEN = ["full", "drop10", "drop20"]
OUT = os.path.join(RES, "rules", ARM)


def read_truth():
    t = {}
    with open(os.path.join(RES, "inputs/truth/truth_windows_coords.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t[(int(r["CHROM"]), int(float(r["start"])))] = r["truth_cluster"]
    return t


def read_filtered(path):
    """-> (set of windows, {cluster: set of windows}, {sample: set of windows}, {sample: cluster})"""
    allw, bycl, bysamp, samp_cl = set(), {}, {}, {}
    if not os.path.exists(path):
        return allw, bycl, bysamp, samp_cl
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (int(r["CHROM"]), int(float(r["START"])))
            allw.add(k)
            bycl.setdefault(r["Cluster"], set()).add(k)
            bysamp.setdefault(r["SAMPLE"], set()).add(k)
            samp_cl[r["SAMPLE"]] = r["Cluster"]
    return allw, bycl, bysamp, samp_cl


def read_mask(path):
    m = set()
    if not os.path.exists(path):
        return m
    for line in open(path).read().splitlines()[1:]:
        g = re.match(r"w(\d+)_(\d+)$", line.strip())
        if g:
            m.add((int(g.group(1)), int(g.group(2)) - 5000))  # BIN -> START
    return m


def cluster_members():
    cl = {}
    f = ("admix_clusters.tsv" if ARM == "full" else "admix_clusters_declonal.tsv")
    with open(os.path.join(RES, "inputs", f)) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cl.setdefault(r["Cluster"], set()).add(r["Sample"])
    return cl


def jaccard(a, b):
    u = len(a | b)
    return (len(a & b) / u) if u else float("nan")


def spearman(xs, ys):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(s):
            r[i] = pos + 1.0
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


truth = read_truth()
T = set(truth)
members = cluster_members()
rows = []

# "floor6" = every rule at the published analysis's floor of 6, the
# configuration the Indo 558-basis scoreboard used, so our ranking is
# directly comparable to it. Stability columns are absent there by design
# (no drop scenarios were aggregated at this floor).
for policy in ("common", "own", "floor6"):
    base_full = None
    for rule in RULES:
        p = os.path.join(OUT, rule, f"agg_{policy}", "full",
                         "introgressed_windows_filtered.tsv")
        allw, bycl, bysamp, _ = read_filtered(p)
        if not allw:
            print(f"  [skip] {rule}/{policy}: no filtered output", file=sys.stderr)
            continue
        if rule == "absolute":
            base_full = (allw, read_mask(os.path.join(
                OUT, rule, f"agg_{policy}", "full", "hypervariable_masked_windows.tsv")))

        # 1. reproduction
        repro = len(allw & T)
        # 2. cluster-size dependence: median windows/sample, zero-filled
        med, sizes = {}, {}
        for cl, mem in members.items():
            counts = [len(bysamp.get(s, ())) for s in mem]
            med[cl] = st.median(counts) if counts else 0.0
            sizes[cl] = len(mem)
        # A cluster with zero surviving windows makes a raw max/min ratio
        # infinite and uninformative, so report the ratio over clusters that
        # retain any signal and carry the empty ones as an explicit count.
        order = sorted(members)
        nz = {c: v for c, v in med.items() if v > 0}
        n_empty = len(med) - len(nz)
        ratio = (max(nz.values()) / min(nz.values())) if len(nz) >= 2 else float("nan")
        rho = spearman([sizes[c] for c in order], [med[c] for c in order])
        # 3. stability
        stab = {}
        for s in ("drop10", "drop20"):
            w2, _, _, _ = read_filtered(os.path.join(
                OUT, rule, f"agg_{policy}", s, "introgressed_windows_filtered.tsv"))
            stab[s] = jaccard(allw, w2) if w2 else float("nan")  # NaN where the scenario was not aggregated
        # filter-4 decomposition vs absolute
        mask = read_mask(os.path.join(OUT, rule, f"agg_{policy}", "full",
                                      "hypervariable_masked_windows.tsv"))
        if base_full and rule != "absolute":
            b_all, b_mask = base_full
            diff = allw ^ b_all
            via_mask = len([w for w in diff if (w in mask) != (w in b_mask)])
        else:
            diff, via_mask = set(), 0

        rows.append(dict(
            arm=ARM, rule=rule, floor_policy=policy,
            called_windows=len(allw),
            reproduced_of_217=repro,
            pct_reproduced=round(100 * repro / len(T), 2),
            jaccard_vs_truth=round(jaccard(allw, T), 4),
            median_win_per_sample=";".join(f"{c}={med[c]:g}" for c in order),
            size_bias_ratio_nonzero=(round(ratio, 2) if ratio == ratio else "NA"),
            clusters_with_zero_windows=n_empty,
            size_bias_rho=round(rho, 3),
            stability_drop10=round(stab["drop10"], 4),
            stability_drop20=round(stab["drop20"], 4),
            hypervariable_masked=len(mask),
            windows_differing_vs_absolute=len(diff),
            of_which_via_mask_movement=via_mask,
        ))

dest = os.path.join(RES, "benchmark", "rule_scoreboard.tsv")
os.makedirs(os.path.dirname(dest), exist_ok=True)
with open(dest, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

hdr = ["rule", "floor_policy", "called_windows", "reproduced_of_217", "pct_reproduced",
       "jaccard_vs_truth", "median_win_per_sample", "size_bias_ratio_nonzero",
       "clusters_with_zero_windows", "size_bias_rho",
       "stability_drop10", "stability_drop20", "hypervariable_masked",
       "windows_differing_vs_absolute", "of_which_via_mask_movement"]
print("\t".join(hdr))
for r in rows:
    print("\t".join(str(r[h]) for h in hdr))
print(f"\nwrote {dest}", file=sys.stderr)
