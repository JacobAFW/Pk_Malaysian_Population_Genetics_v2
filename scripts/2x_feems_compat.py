#!/usr/bin/env python
"""scikit-sparse API shim for feems 1.0.0.

feems calls the scikit-sparse 0.4.x CHOLMOD interface:

    factor = cholmod.cholesky(A)      # -> callable Factor object
    x      = factor(b)                # solve A x = b
    factor = factor.cholesky(A2)      # refactorise, same sparsity pattern

scikit-sparse 0.5.0 removed that interface. `cholesky(A)` now returns a plain
`(R, P)` tuple, so `feems/objective.py::_rank_one_solver` fails with
`TypeError: 'tuple' object is not callable` inside `fit_null_model`. The
factor-object API moved to `cho_factor` -> `CholeskyFactor.solve/.factorize`.

This module installs a 0.4-shaped wrapper into `feems.spatial_graph.cholmod`
when (and only when) the installed scikit-sparse has the 0.5 signature. It is a
pure adapter -- CHOLMOD does the identical numerical work either way.

This is NOT the feems #30 convergence bug (monomorphic sites); it is an
independent dependency-version break. Both have to be fixed for a fit to run.

Usage (before constructing a SpatialGraph):
    from feems_compat import patch_cholmod; patch_cholmod()
"""
from __future__ import annotations

import numpy as np


class _Factor:
    """0.4.x-style callable CHOLMOD factor backed by 0.5.x CholeskyFactor."""

    __slots__ = ("_f",)

    def __init__(self, f):
        self._f = f

    def __call__(self, b):
        b = np.asarray(b) if not hasattr(b, "toarray") else b
        one_d = getattr(b, "ndim", 2) == 1
        x = self._f.solve(b.reshape(-1, 1) if one_d else b)
        x = np.asarray(x.toarray() if hasattr(x, "toarray") else x)
        return x.ravel() if one_d else x

    def cholesky(self, A):
        """Refactorise in place and return self (0.4.x returned a new Factor)."""
        self._f.factorize(A)
        return self

    def logdet(self):
        return self._f.logdet()

    def __getattr__(self, name):
        return getattr(self._f, name)


class _CholmodShim:
    def __init__(self, mod):
        self._mod = mod

    def cholesky(self, A):
        return _Factor(self._mod.cho_factor(A))

    def __getattr__(self, name):
        return getattr(self._mod, name)


def needs_patch():
    """True if the installed scikit-sparse returns a tuple from cholesky()."""
    import scipy.sparse as sp
    import sksparse.cholmod as cholmod

    probe = sp.eye(3, format="csc") * 2.0
    return isinstance(cholmod.cholesky(probe), tuple)


def patch_cholmod(verbose=True):
    """Install the shim into feems.spatial_graph if the 0.5 API is present.

    Returns a short description of what was done, for the provenance record.
    """
    import importlib.metadata as md

    import feems.spatial_graph as fsg

    try:
        ver = md.version("scikit-sparse")
    except md.PackageNotFoundError:
        ver = "unknown"

    if not needs_patch():
        msg = f"scikit-sparse {ver}: 0.4-style API, no shim needed"
    else:
        import sksparse.cholmod as cholmod

        fsg.cholmod = _CholmodShim(cholmod)
        msg = (f"scikit-sparse {ver}: 0.5-style API detected; installed callable-Factor "
               f"shim into feems.spatial_graph.cholmod (cho_factor/solve/factorize)")
    if verbose:
        print(f"  compat: {msg}")
    return msg


def self_test():
    """Verify the shim solves the same systems as a dense reference."""
    import scipy.sparse as sp
    import sksparse.cholmod as cholmod

    rng = np.random.default_rng(0)
    n = 40
    A = (sp.eye(n) * 4 - sp.eye(n, k=1) - sp.eye(n, k=-1)).tocsc()
    f = _CholmodShim(cholmod).cholesky(A)
    B = rng.normal(size=(n, 5))
    ones = np.ones(n)
    assert np.allclose(f(B), np.linalg.solve(A.toarray(), B)), "matrix solve mismatch"
    v = f(ones)
    assert v.ndim == 1 and np.allclose(v, np.linalg.solve(A.toarray(), ones)), "vector solve mismatch"
    A2 = (sp.eye(n) * 6 - sp.eye(n, k=1) - sp.eye(n, k=-1)).tocsc()
    f2 = f.cholesky(A2)
    assert np.allclose(f2(B), np.linalg.solve(A2.toarray(), B)), "refactorisation mismatch"
    print("  compat self-test: OK (solve, 1-D solve, refactorise all match dense reference)")


if __name__ == "__main__":
    print(patch_cholmod())
    self_test()
