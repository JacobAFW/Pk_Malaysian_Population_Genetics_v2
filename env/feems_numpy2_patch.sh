#!/usr/bin/env bash
# feems_numpy2_patch.sh — apply the 3-line NumPy 2 compat patch to an
# installed bioconda feems-2.0.1. Idempotent. Run after any reinstall of
# feems into the feems_e env. See env/ENV-PLAN.md + MEMORY.md for context.
#
# Why this patch exists:
# - bioconda feems-2.0.1 declares numpy<2, but the conda-forge tskit 1.0.3
#   and msprime 1.4.2 binaries that get pulled in are compiled against the
#   NumPy 2 C-API. So with numpy<2 the import fails with a C-API ABI error.
# - With numpy>=2, the tskit/msprime ABI is satisfied but feems source still
#   uses np.Inf (removed in NumPy 2.0) and np.in1d (removed in NumPy 2.0).
# - Upstream NovembreLab/feems and VivaswatS/feems HEAD have NOT been
#   patched as of Phase 3 install. The 3 lines below are the minimal fix.

set -euo pipefail

ENV_PREFIX="${ENV_PREFIX:-$HOME/mamba/envs/feems_e}"
FEEMS_DIR="$ENV_PREFIX/lib/python3.12/site-packages/feems"

if [[ ! -d "$FEEMS_DIR" ]]; then
  echo "ERROR: $FEEMS_DIR not found — install feems first." >&2
  exit 1
fi

# Patch 1+2: spatial_graph.py:759-760  np.Inf  -> np.inf
/usr/bin/sed -i.bak -e 's/lb=-np\.Inf,/lb=-np.inf,/' \
                    -e 's/ub=np\.Inf,/ub=np.inf,/' \
                    "$FEEMS_DIR/spatial_graph.py"

# Patch 3: helper_funcs.py:67  np.in1d  -> np.isin
/usr/bin/sed -i.bak -e 's/np\.in1d(/np.isin(/' \
                    "$FEEMS_DIR/helper_funcs.py"

# Sanity check: no remaining np.Inf or np.in1d
if /usr/bin/grep -qE 'np\.Inf|np\.in1d' "$FEEMS_DIR"/*.py; then
  echo "ERROR: residual np.Inf or np.in1d in $FEEMS_DIR" >&2
  /usr/bin/grep -nE 'np\.Inf|np\.in1d' "$FEEMS_DIR"/*.py >&2
  exit 1
fi
echo "==> feems NumPy-2 patch applied (3 sites)."
