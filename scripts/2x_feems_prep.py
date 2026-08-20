#!/usr/bin/env python
"""Phase 5 / FEEMS step 1: assemble the sample keep-list, coordinates and region polygon.

Emits, for a given PLINK base and EEMS keep-list:
  <out>/<tag>_keep.txt        FID IID, for `plink --keep`
  <out>/<tag>_coords.tsv      IID<TAB>long<TAB>lat for every kept sample (ID-keyed, NOT positional)
  <out>/<tag>_outer.csv       long,lat ring of the outer region polygon
  <out>/<tag>_prep.json       provenance + counts

Sample set = (EEMS keep-list) INTERSECT (PLINK cohort) INTERSECT (inside outer polygon).
The spatial filter is applied *here* so that the monomorphic-site drop downstream
(`plink --maf`) is computed on the final sample set, not a superset.

Coordinates are keyed by sample ID and re-joined against the .fam order in the fit
script; nothing downstream relies on the row order of the input .coords file.

Run: micromamba run -n feems_e python scripts/2x_feems_prep.py --help
"""
from __future__ import annotations

import argparse
import json
import os

import fiona
import numpy as np
from shapely.geometry import Point, shape
from shapely.ops import unary_union


def read_fam_ids(path):
    """Return [(FID, IID), ...] in .fam order."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if f:
                out.append((f[0], f[1]))
    return out


def read_coord_table(path):
    """Read a .coords / .fam-style table: FID IID 0 0 0 -9 LONG LAT -> {IID: (long, lat)}."""
    coords = {}
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            f = line.split()
            if not f:
                continue
            if len(f) < 8:
                raise ValueError(f"{path}:{ln}: expected >=8 columns (FID IID .. LONG LAT), got {len(f)}")
            coords[f[1]] = (float(f[6]), float(f[7]))
    return coords


def load_outer(shp, buffer_deg):
    """Union the region shapefile, buffer it, and return a single exterior ring as (q, 2)."""
    geoms = [shape(rec["geometry"]) for rec in fiona.open(shp)]
    poly = unary_union(geoms)
    if buffer_deg:
        poly = poly.buffer(buffer_deg)
    if poly.geom_type == "MultiPolygon":
        # feems `prepare_graph_inputs` builds a single shapely Polygon from `outer`,
        # so a disjoint region cannot be represented. Take the largest part and say so.
        parts = sorted(poly.geoms, key=lambda g: g.area, reverse=True)
        print(
            f"  WARNING: outer is a MultiPolygon after buffer={buffer_deg} "
            f"({len(parts)} parts); keeping the largest (area={parts[0].area:.3f}, "
            f"dropped area={sum(p.area for p in parts[1:]):.3f})"
        )
        poly = parts[0]
    ring = np.asarray(poly.exterior.coords)
    return poly, ring


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bfile", required=True, help="PLINK base path (expects <bfile>.fam)")
    p.add_argument("--keeplist", required=True, help="EEMS keep-list (.fam-style, IID in col 2)")
    p.add_argument("--coords", required=True, help="coordinate table (.fam-style, LONG col 7, LAT col 8)")
    p.add_argument("--outer-shp", required=True, help="region polygon shapefile (EPSG:4326)")
    p.add_argument("--buffer", type=float, default=0.25,
                   help="degrees to buffer the region polygon (default 0.25 ~ 27 km)")
    p.add_argument("--cluster-labels", default=None,
                   help="TSV with header, col1=sample col2=cluster; used with --cluster")
    p.add_argument("--cluster", default=None,
                   help="restrict to this cluster label (e.g. Mf). FEEMS models ONE spatially "
                        "structured population; sympatric divergent lineages violate that.")
    p.add_argument("--out-dir", default="data/processed/feems")
    p.add_argument("--tag", required=True, help="output file prefix")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    stem = os.path.join(args.out_dir, args.tag)

    fam = read_fam_ids(args.bfile + ".fam")
    cohort = {iid: fid for fid, iid in fam}
    keep_ids = [iid for _, iid in read_fam_ids(args.keeplist)]
    coords = read_coord_table(args.coords)
    poly, ring = load_outer(args.outer_shp, args.buffer)

    print(f"  PLINK cohort           : {len(fam)}")
    print(f"  EEMS keep-list         : {len(keep_ids)} ({len(set(keep_ids))} unique)")
    print(f"  coordinate table       : {len(coords)}")
    print(f"  outer polygon bounds   : {tuple(round(v, 4) for v in poly.bounds)}")

    in_cohort = [i for i in keep_ids if i in cohort]
    print(f"  keep-list INTERSECT cohort : {len(in_cohort)}  (dropped {len(keep_ids) - len(in_cohort)})")

    n_before_cluster = len(in_cohort)
    if args.cluster:
        if not args.cluster_labels:
            raise SystemExit("ERROR: --cluster requires --cluster-labels")
        labels = {}
        with open(args.cluster_labels) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    labels[f[0]] = f[1]
        in_cohort = [i for i in in_cohort if labels.get(i) == args.cluster]
        print(f"  restricted to cluster '{args.cluster}': {len(in_cohort)} "
              f"(dropped {n_before_cluster - len(in_cohort)})")

    missing_coord = [i for i in in_cohort if i not in coords]
    with_coord = [i for i in in_cohort if i in coords]
    if missing_coord:
        print(f"  WARNING: {len(missing_coord)} kept samples have no coordinate -> dropped")

    inside, outside = [], []
    for iid in with_coord:
        (inside if poly.intersects(Point(*coords[iid])) else outside).append(iid)
    print(f"  inside outer polygon   : {len(inside)}  (outside {len(outside)})")
    if outside:
        pts = sorted({(round(coords[i][0], 1), round(coords[i][1], 1)) for i in outside})
        print(f"  dropped (outside region) at approx long/lat: {pts}")

    if not inside:
        raise SystemExit("ERROR: no samples fall inside the outer polygon")

    with open(stem + "_keep.txt", "w") as fh:
        for iid in inside:
            fh.write(f"{cohort[iid]}\t{iid}\n")
    with open(stem + "_coords.tsv", "w") as fh:
        for iid in inside:
            lo, la = coords[iid]
            fh.write(f"{iid}\t{lo!r}\t{la!r}\n")
    np.savetxt(stem + "_outer.csv", ring, delimiter=",", header="long,lat", comments="")

    xy = np.array([coords[i] for i in inside])
    meta = {
        "tag": args.tag,
        "bfile": args.bfile,
        "keeplist": args.keeplist,
        "coords_file": args.coords,
        "outer_shp": args.outer_shp,
        "outer_buffer_deg": args.buffer,
        "outer_bounds": [float(v) for v in poly.bounds],
        "outer_ring_vertices": int(ring.shape[0]),
        "n_plink_cohort": len(fam),
        "n_keeplist": len(keep_ids),
        "cluster_restriction": args.cluster,
        "n_keeplist_in_cohort": n_before_cluster,
        "n_after_cluster_restriction": len(in_cohort) + len(missing_coord) + len(outside),
        "n_missing_coord": len(missing_coord),
        "n_outside_region": len(outside),
        "n_final": len(inside),
        "long_range": [float(xy[:, 0].min()), float(xy[:, 0].max())],
        "lat_range": [float(xy[:, 1].min()), float(xy[:, 1].max())],
        "n_unique_localities": int(np.unique(xy, axis=0).shape[0]),
    }
    with open(stem + "_prep.json", "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"  FINAL n_samples        : {len(inside)}")
    print(f"  unique localities      : {meta['n_unique_localities']}")
    print(f"  wrote {stem}_{{keep.txt,coords.tsv,outer.csv,prep.json}}")


if __name__ == "__main__":
    main()
