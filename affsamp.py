"""Exact computations for affine-flat (coset) sampling designs on F_p^m, p prime.

Conventions
-----------
* Units are the vectors x of G = F_p^m, encoded as integers in C order:
  code(x) = sum_j x_j p^(m-1-j).  Arrays of length N = p^m are indexed by code.
* A factor (stratification) is given by a subspace W of G; its levels are the
  cosets x + W.  We store the factor through a basis of its annihilator
  W^perp (a subspace of the dual, identified with G by the dot product).
* A coset design is given by a probability vector nu over "dual directions"
  D (subspaces of dimension m - h of the dual); the sampled flat has direction
  V = D^perp (dimension h).  r distinct cosets of V are drawn by simple random
  sampling, so the sample size is k = r p^h.
* Fourier coefficients: yhat(u) = N^-1 sum_x y(x) exp(-2 pi i u.x / p), which is
  numpy.fft.fftn(Y) / N for Y = y.reshape((p,)*m).
"""
from __future__ import annotations

import itertools
import math
import numpy as np


# --------------------------------------------------------------------------
# Linear algebra over GF(p)
# --------------------------------------------------------------------------
def rref_mod(A, p):
    """Reduced row echelon form of A over GF(p). Returns (R, pivots)."""
    R = np.array(A, dtype=np.int64) % p
    if R.ndim == 1:
        R = R.reshape(1, -1)
    rows, cols = R.shape
    pivots = []
    r = 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(R[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + nz[0]
        R[[r, i]] = R[[i, r]]
        R[r] = (R[r] * pow(int(R[r, c]), -1, p)) % p
        for i2 in range(rows):
            if i2 != r and R[i2, c]:
                R[i2] = (R[i2] - R[i2, c] * R[r]) % p
        pivots.append(c)
        r += 1
    return R[:r], pivots


def rank_mod(A, p):
    return len(rref_mod(A, p)[1]) if np.size(A) else 0


def nullspace_mod(A, p, m=None):
    """Basis (rows) of {x in F_p^m : A x = 0}."""
    A = np.array(A, dtype=np.int64)
    if A.size == 0:
        if m is None:
            raise ValueError("m is required when A is empty")
        return np.eye(int(m), dtype=np.int64)
    m = A.shape[-1]
    R, piv = rref_mod(A, p)
    free = [c for c in range(m) if c not in piv]
    basis = []
    for f in free:
        x = np.zeros(m, dtype=np.int64)
        x[f] = 1
        for i, pc in enumerate(piv):
            x[pc] = (-R[i, f]) % p
        basis.append(x)
    return np.array(basis, dtype=np.int64).reshape(len(basis), m)


def subspaces(m, k, p):
    """Yield RREF bases (k x m) of all k-dimensional subspaces of F_p^m."""
    if k == 0:
        yield np.zeros((0, m), dtype=np.int64)
        return
    for piv in itertools.combinations(range(m), k):
        free_pos = [(i, j) for i in range(k) for j in range(piv[i] + 1, m) if j not in piv]
        for vals in itertools.product(range(p), repeat=len(free_pos)):
            B = np.zeros((k, m), dtype=np.int64)
            for i, c in enumerate(piv):
                B[i, c] = 1
            for (i, j), v in zip(free_pos, vals):
                B[i, j] = v
            yield B


def gaussian_binomial(m, k, q):
    if k < 0 or k > m:
        return 0
    num = den = 1
    for i in range(k):
        num *= q ** (m - i) - 1
        den *= q ** (i + 1) - 1
    return num // den


# --------------------------------------------------------------------------
# The ambient space
# --------------------------------------------------------------------------
class Space:
    def __init__(self, p, m):
        self.p, self.m, self.N = p, m, p ** m
        self.vecs = np.array(list(itertools.product(range(p), repeat=m)), dtype=np.int64)
        self.place = p ** np.arange(m - 1, -1, -1)

    def code(self, v):
        return (np.asarray(v, dtype=np.int64) % self.p) @ self.place

    def span(self, basis):
        """Sorted codes of the subspace spanned by the rows of basis."""
        basis = np.asarray(basis, dtype=np.int64).reshape(-1, self.m)
        if basis.shape[0] == 0:
            return np.array([0], dtype=np.int64)
        coeffs = np.array(list(itertools.product(range(self.p), repeat=basis.shape[0])), dtype=np.int64)
        return np.unique(self.code((coeffs @ basis) % self.p))

    def dot(self, u_codes, x_codes):
        return (self.vecs[u_codes] @ self.vecs[x_codes].T) % self.p

    def fourier_power(self, y):
        """|yhat(u)|^2 for all u (array of length N, indexed by code)."""
        Y = np.asarray(y, dtype=float).reshape((self.p,) * self.m)
        yhat = np.fft.fftn(Y) / self.N
        return (np.abs(yhat) ** 2).reshape(-1)


# --------------------------------------------------------------------------
# Factors, main-effect set, admissible directions
# --------------------------------------------------------------------------
def main_effect_set(space, factors):
    """M_S: nonzero vectors of the union of the annihilators W^perp (boolean mask)."""
    mask = np.zeros(space.N, dtype=bool)
    for Wp in factors:
        mask[space.span(Wp)] = True
    mask[0] = False
    return mask


def admissible_directions(space, factors, h):
    """All dual directions D (dim m-h) with D cap M_S = {0}, i.e. V = D^perp admissible."""
    M = main_effect_set(space, factors)
    out = []
    for D in subspaces(space.m, space.m - h, space.p):
        if not M[space.span(D)].any():
            out.append(D)
    return out


def coset_labels(space, D):
    """Label of the coset of V = D^perp containing each unit: x -> D x (as an integer)."""
    k = D.shape[0]
    lab = (space.vecs @ D.T) % space.p
    return lab @ (space.p ** np.arange(k - 1, -1, -1)) if k else np.zeros(space.N, dtype=np.int64)


# --------------------------------------------------------------------------
# Designs
# --------------------------------------------------------------------------
class CosetDesign:
    """V ~ nu over admissible directions, then r of the Q = p^(m-h) cosets of V by SRS."""

    def __init__(self, space, directions, probs=None, r=1):
        self.space = space
        self.D = list(directions)
        n = len(self.D)
        self.nu = np.full(n, 1.0 / n) if probs is None else np.asarray(probs, dtype=float)
        self.nu = self.nu / self.nu.sum()
        self.r = r
        self.mh = self.D[0].shape[0]  # m - h
        self.h = space.m - self.mh
        self.Q = space.p ** self.mh
        self.k = r * space.p ** self.h
        self._cover = None

    def coverage(self):
        """c(u) = Pr_nu(u in D) for every u (array of length N)."""
        if self._cover is None:
            c = np.zeros(self.space.N)
            for D, pr in zip(self.D, self.nu):
                if pr > 0:
                    c[self.space.span(D)] += pr
            self._cover = c
        return self._cover

    def weights(self):
        """Spectral weights w(u) = E|eta_s(u)|^2, u != 0; w(0) is set to 0."""
        Q, r = self.Q, self.r
        w = self.coverage() * (Q - r) / (r * (Q - 1))
        w = w.copy()
        w[0] = 0.0
        return w

    def variance(self, y):
        """Design variance of the sample mean (spectral formula)."""
        return float(np.dot(self.weights(), self.space.fourier_power(y)))

    def pair_probability(self):
        """pi(z) = P(x in s, x+z in s) for z != 0 (closed form)."""
        Q, r = self.Q, self.r
        # Pr(z in V) for V = D^perp: z in V iff D z = 0
        pz = np.zeros(self.space.N)
        for D, pr in zip(self.D, self.nu):
            if pr > 0:
                inV = ((self.space.vecs @ D.T) % self.space.p == 0).all(axis=1)
                pz[inV] += pr
        pi1 = self.k / self.space.N
        return pi1 * (pz + (1 - pz) * (r - 1) / (Q - 1))

    def samples(self):
        """Yield (probability, sample codes) over the whole support (small cases only)."""
        for D, pr in zip(self.D, self.nu):
            if pr == 0:
                continue
            lab = coset_labels(self.space, D)
            ncomb = math.comb(self.Q, self.r)
            for sel in itertools.combinations(range(self.Q), self.r):
                s = np.nonzero(np.isin(lab, sel))[0]
                yield pr / ncomb, s


def brute_force_moments(design, y):
    """Exact mean and variance of the sample mean by enumerating the support."""
    y = np.asarray(y, dtype=float)
    m1 = m2 = 0.0
    for pr, s in design.samples():
        v = y[s].mean()
        m1 += pr * v
        m2 += pr * v * v
    return m1, m2 - m1 * m1


def brute_force_pi2(design):
    N = design.space.N
    P = np.zeros((N, N))
    for pr, s in design.samples():
        ind = np.zeros(N)
        ind[s] = 1
        P += pr * np.outer(ind, ind)
    return P


def is_balanced(space, s, factors):
    """True if sample s has the same count in every level of every factor."""
    for Wp in factors:
        Wp = np.asarray(Wp).reshape(-1, space.m)
        lab = coset_labels(space, Wp)  # levels of W are the values of x -> W^perp x
        counts = np.bincount(lab[s], minlength=space.p ** Wp.shape[0])
        if counts.min() != counts.max():
            return False
    return True


# --------------------------------------------------------------------------
# Population summaries
# --------------------------------------------------------------------------
def decomposition(space, factors, y):
    """(E_M, sigma_R^2, |M_S|, |R_S|): main-effect and residual variances (divisor N)."""
    pw = space.fourier_power(y)
    M = main_effect_set(space, factors)
    R = ~M
    R[0] = False
    return float(pw[M].sum()), float(pw[R].sum()), int(M.sum()), int(R.sum())


def srs_variance(N, k, y):
    y = np.asarray(y, dtype=float)
    S2 = y.var(ddof=1)
    return (1 - k / N) * S2 / k


def adjusted_r2(space, factors, y):
    EM, sR, _, dR = decomposition(space, factors, y)
    N = space.N
    return 1 - (sR / dR) / ((EM + sR) / (N - 1))
