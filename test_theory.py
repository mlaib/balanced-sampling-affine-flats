"""Brute-force checks of every exact statement used in the manuscript.

Run:  python3 -I test_theory.py   (from this directory: python3 test_theory.py)
"""
import itertools
import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import (Space, CosetDesign, admissible_directions, brute_force_moments,
                     brute_force_pi2, coset_labels, is_balanced,
                     main_effect_set, nullspace_mod, rank_mod, subspaces, srs_variance,
                     adjusted_r2)

rng = np.random.default_rng(20261015)
TOL = 1e-10
fails = 0
counts = {}


def check(cond, msg):
    global fails
    key = msg.split(" p=")[0]
    counts[key] = counts.get(key, 0) + 1
    if not cond:
        fails += 1
        print("FAIL:", msg)


def random_factors(sp, nfac, codims):
    """Random factors: annihilator bases of random dimension in codims."""
    fac = []
    for _ in range(nfac):
        n = rng.choice(codims)
        while True:
            B = rng.integers(0, sp.p, size=(n, sp.m))
            if rank_mod(B, sp.p) == n:
                break
        fac.append(B)
    return fac


# ---------------------------------------------------------------- Lemma 1 (balance <=> V + W = G)
for (p, m) in [(2, 4), (3, 3), (2, 5)]:
    sp = Space(p, m)
    for trial in range(4):
        fac = random_factors(sp, 2, [1, 2] if m > 3 else [1])
        nmax = max(f.shape[0] for f in fac)
        for h in range(nmax, m):
            for D in subspaces(m, m - h, p):
                V = nullspace_mod(D, p)
                lab = coset_labels(sp, D)
                s = np.nonzero(lab == lab[rng.integers(sp.N)])[0]  # one coset of V
                bal = is_balanced(sp, s, fac)
                Dspan = sp.span(D)
                dual_ok = all(np.intersect1d(Dspan, sp.span(f)).size == 1 for f in fac)
                primal_ok = all(rank_mod(np.vstack([V, nullspace_mod(f, p)]), p) == m for f in fac)
                check(bal == dual_ok == primal_ok, f"Lemma 1 p={p} m={m} h={h}")
print("Lemma 1 (balance iff V+W=G iff D cap W^perp = 0): done")

# ---------------------------------------------------------------- Theorem 1 + Prop (pi, w) on random designs
for (p, m, h, r) in [(2, 4, 2, 1), (2, 4, 2, 2), (2, 5, 3, 1), (2, 5, 2, 3), (3, 3, 1, 1), (3, 3, 1, 2), (3, 3, 2, 2), (5, 2, 1, 2)]:
    sp = Space(p, m)
    fac = random_factors(sp, 2, [1])
    A = admissible_directions(sp, fac, h)
    if not A:
        continue
    nu = rng.random(len(A)) ** 3
    des = CosetDesign(sp, A, nu, r=r)
    w = des.weights()
    # conservation law
    check(abs(w.sum() - (sp.N / des.k - 1)) < 1e-9, f"conservation p={p} m={m} h={h} r={r}")
    # zero weight on the main-effect set
    M = main_effect_set(sp, fac)
    check(np.all(w[M] < TOL), "w vanishes on M_S")
    # every sample balanced, fixed size
    for pr, s in des.samples():
        check(len(s) == des.k and is_balanced(sp, s, fac), "support balanced")
    # spectral variance vs brute force, several y
    for _ in range(3):
        y = rng.normal(size=sp.N) + 3 * rng.random()
        mean, var = brute_force_moments(des, y)
        check(abs(mean - y.mean()) < 1e-9, "unbiased")
        check(abs(var - des.variance(y)) < 1e-9 * max(1, var), f"spectral variance p={p} m={m} h={h} r={r}")
    # second-order inclusion probabilities
    P = brute_force_pi2(des)
    pz = des.pair_probability()
    x0 = 0
    for x in range(sp.N):
        z = sp.code(sp.vecs[x] - sp.vecs[x0])
        if x != x0:
            check(abs(P[x0, x] - pz[z]) < 1e-12, "pi_xy closed form")
    check(np.allclose(np.diag(P), des.k / sp.N), "pi_x = k/N")
    # Fourier link pi(z) = pi^2 (1 + sum_u w(u) chi_u(z))
    pi1 = des.k / sp.N
    phase = np.exp(2j * np.pi * sp.dot(np.arange(sp.N), np.arange(sp.N)) / p)  # [u, z]
    pz_from_w = pi1 ** 2 * (1 + (w[:, None] * phase).sum(axis=0)).real
    check(np.allclose(pz_from_w[1:], pz[1:]), "pi(z) from spectral weights")
    check(abs(pz_from_w[0] - pi1) < 1e-12, "pi(0) from spectral weights")
print("Spectral variance, conservation, pi_xy, balance of support: done")

# ---------------------------------------------------------------- single factor: uniform coset design == stratified SRS
for (p, m, n, h) in [(2, 4, 1, 2), (2, 4, 2, 2), (2, 5, 2, 3), (3, 3, 1, 1), (3, 3, 1, 2), (3, 4, 2, 3)]:
    sp = Space(p, m)
    Wp = random_factors(sp, 1, [n])
    A = admissible_directions(sp, Wp, h)
    des = CosetDesign(sp, A)
    w = des.weights()
    M = main_effect_set(sp, Wp)
    R = ~M
    R[0] = False
    flat = (p ** (m - h) - 1) / (sp.N - 1 - M.sum())
    check(np.allclose(w[R], flat), f"single factor flat p={p} m={m} n={n} h={h}")
    # stratified SRS with proportional allocation: strata = cosets of W, k_h = k / p^n each
    lab = coset_labels(sp, Wp[0])
    k = des.k
    Nh, kh = sp.N // p ** n, k // p ** n
    for _ in range(3):
        y = rng.normal(size=sp.N)
        v_str = sum((1 - kh / Nh) * y[lab == g].var(ddof=1) / kh * (Nh / sp.N) ** 2 for g in range(p ** n))
        check(abs(v_str - des.variance(y)) < 1e-10, "coset design == stratified SRS variance")
print("Single factor: uniform admissible coset design has the second-order behaviour of stratified SRS: done")

# ---------------------------------------------------------------- single factor with r >= 2 cosets == stratified SRS; deff independent of r
for (p, m, n, h, r) in [(2, 4, 1, 2, 2), (2, 5, 2, 3, 3), (3, 3, 1, 2, 2), (3, 4, 1, 2, 4)]:
    sp = Space(p, m)
    Wp = random_factors(sp, 1, [n])
    A = admissible_directions(sp, Wp, h)
    des = CosetDesign(sp, A, r=r)
    lab = coset_labels(sp, Wp[0])
    k = des.k
    Nh, kh = sp.N // p ** n, k // p ** n
    for _ in range(3):
        y = rng.normal(size=sp.N)
        v_str = sum((1 - kh / Nh) * y[lab == g].var(ddof=1) / kh * (Nh / sp.N) ** 2 for g in range(p ** n))
        check(abs(v_str - des.variance(y)) < 1e-10, f"r-coset design == stratified SRS p={p} m={m} r={r}")
    # brute-force pi_xy against stratified SRS second-order probabilities
    P = brute_force_pi2(des)
    same = lab[:, None] == lab[None, :]
    pi_str = np.where(same, kh * (kh - 1) / (Nh * (Nh - 1)), (kh / Nh) ** 2)
    np.fill_diagonal(pi_str, kh / Nh)
    check(np.allclose(P, pi_str), f"r-coset pi_xy == stratified SRS pi_xy p={p} m={m} r={r}")
for (p, m, h) in [(2, 5, 3), (3, 3, 2)]:
    sp = Space(p, m)
    fac = random_factors(sp, 2, [1])
    A = admissible_directions(sp, fac, h)
    nu = rng.random(len(A))
    y = rng.normal(size=sp.N)
    Q = p ** (m - h)
    deffs = []
    for r in range(1, Q):
        d = CosetDesign(sp, A, nu, r=r)
        deffs.append(d.variance(y) / srs_variance(sp.N, d.k, y))
    check(np.allclose(deffs, deffs[0]), f"deff independent of r p={p} m={m}")
print("Single factor r>=2 == stratified SRS; design effect independent of r: done")

# ---------------------------------------------------------------- flat designs: m-h = 1, m = 2 (any S), spreads
for (p, m) in [(2, 4), (3, 3), (2, 5)]:
    sp = Space(p, m)
    fac = random_factors(sp, 3, [1, 2] if m > 3 else [1])
    A = admissible_directions(sp, fac, m - 1)
    if A:
        des = CosetDesign(sp, A)
        w = des.weights()
        M = main_effect_set(sp, fac)
        R = ~M
        R[0] = False
        check(np.allclose(w[R], w[R][0]), f"m-h=1 flat p={p} m={m}")
for p in [3, 5, 7]:
    sp = Space(p, 2)
    lines = [B for B in subspaces(2, 1, p)]
    for nS in range(0, p + 1):
        S = lines[:nS]
        A = admissible_directions(sp, S, 1) if S else lines
        des = CosetDesign(sp, A)
        w = des.weights()
        M = main_effect_set(sp, S) if S else np.zeros(sp.N, bool)
        R = ~M
        R[0] = False
        check(np.allclose(w[R], 1 / (p + 1 - nS)), f"m=2 flat p={p} |S|={nS}")
print("Flat designs (m-h=1; m=2 with any set of directions): done")

# ---------------------------------------------------------------- adjusted R^2 identity for flat designs
for p in [5, 7]:
    sp = Space(p, 2)
    lines = [B for B in subspaces(2, 1, p)]
    for nS in [1, 2, 3]:
        S = lines[:nS]
        des = CosetDesign(sp, admissible_directions(sp, S, 1))
        for _ in range(3):
            y = rng.normal(size=sp.N) + np.repeat(rng.normal(size=p), p) * rng.random()
            ratio = des.variance(y) / srs_variance(sp.N, des.k, y)
            check(abs(ratio - (1 - adjusted_r2(sp, S, y))) < 1e-10, "deff = 1 - adjR2")
print("Design effect of flat designs = 1 - adjusted R^2: done")

# ---------------------------------------------------------------- uniform Latin (permutation) sampling has the same weights
for p in [5, 7]:
    sp = Space(p, 2)
    w = np.zeros(sp.N, dtype=float)
    perms = list(itertools.permutations(range(p)))
    for sig in perms:
        s = [sp.code([i, sig[i]]) for i in range(p)]
        ind = np.zeros(sp.N)
        ind[s] = 1
        eta = np.fft.fftn(ind.reshape(p, p)).reshape(-1) / p
        w += np.abs(eta) ** 2 / len(perms)
    w[0] = 0
    rows_cols = [np.array([[1, 0]]), np.array([[0, 1]])]
    des = CosetDesign(sp, admissible_directions(sp, rows_cols, 1))
    check(np.allclose(w, des.weights()), f"Latin == flat coset p={p}")
print("Uniform one-per-row-and-column sampling == flat coset design (rows, columns): done")

# ---------------------------------------------------------------- r-coset variance estimator is unbiased
for (p, m, h, r) in [(2, 4, 1, 2), (2, 4, 1, 3), (3, 3, 1, 2), (2, 5, 2, 2), (5, 2, 1, 3)]:
    sp = Space(p, m)
    fac = random_factors(sp, 2, [1])
    A = admissible_directions(sp, fac, h)
    if not A:
        continue
    des = CosetDesign(sp, A, rng.random(len(A)), r=r)
    y = rng.normal(size=sp.N)
    Ev = 0.0
    for D, pr in zip(des.D, des.nu):
        lab = coset_labels(sp, D)
        cm = np.array([y[lab == g].mean() for g in range(des.Q)])
        for sel in itertools.combinations(range(des.Q), r):
            c = cm[list(sel)]
            v = (1 - r / des.Q) * c.var(ddof=1) / r
            Ev += pr / math.comb(des.Q, r) * v
    check(abs(Ev - des.variance(y)) < 1e-10, f"variance estimator unbiased p={p} m={m} h={h} r={r}")
print("Unbiasedness of the between-coset variance estimator (r >= 2): done")

# ---------------------------------------------------------------- feasibility (Bose-Burton) and positivity
nfeas = ninf = 0
for (p, m) in [(2, 4), (2, 5), (3, 3), (3, 4)]:
    sp = Space(p, m)
    for trial in range(40):
        nfac = rng.integers(1, 6)
        fac = random_factors(sp, nfac, [1, 2] if m >= 4 else [1])
        M = main_effect_set(sp, fac)
        npts = M.sum() // (p - 1)  # points of P_S in PG(m-1, p)
        for h in range(max(f.shape[0] for f in fac), m):
            theta = (p ** (h + 1) - 1) // (p - 1)
            A = admissible_directions(sp, fac, h)
            if npts < theta:
                check(len(A) > 0, "Bose-Burton sufficient condition")
            if len(A) > 0:
                nfeas += 1
            else:
                ninf += 1
            if npts == theta and not A:
                # extremal case: P_S must be a projective h-space (vector dim h+1)
                pts = np.nonzero(M)[0]
                check(rank_mod(sp.vecs[pts], p) == h + 1, "Bose-Burton extremal case")
            if A:
                des = CosetDesign(sp, A)
                pz = des.pair_probability()
                if npts < (p ** h - 1) // (p - 1):
                    check(np.all(pz[1:] > 0), "positivity sufficient condition")
                nmax = max(f.shape[0] for f in fac)
                if h == nmax:
                    for f in fac:
                        if f.shape[0] == h:
                            Wcodes = sp.span(nullspace_mod(f, p))
                            check(np.all(pz[Wcodes[Wcodes > 0]] == 0), "h = n_W gives pi = 0 within levels")
print(f"Feasibility and positivity: done ({nfeas} feasible, {ninf} infeasible cases examined)")

# extremal infeasible example: union of W^perp equal to an (h+1)-space
sp = Space(2, 4)
T = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])  # 3-dim
fac = [np.array([[1, 0, 0, 0], [0, 1, 0, 0]]), np.array([[0, 0, 1, 0], [1, 1, 0, 0]]),
       np.array([[1, 0, 1, 0], [0, 1, 1, 0]])]
M = main_effect_set(sp, fac)
check(set(np.nonzero(M)[0]) == set(sp.span(T)) - {0}, "constructed union is a 3-space")
check(len(admissible_directions(sp, fac, 2)) == 0, "union of W^perp = (h+1)-space blocks")
print("Extremal infeasible example: done")

# deterministic tests of Theorem 1(c): |P_S| = theta_h, blocking iff M_S + {0} is an (h+1)-space
for (p, m, h, pts, is_subspace) in [
        (2, 4, 1, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]], False),          # 3 non-collinear points
        (2, 4, 1, [[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]], True),           # a projective line
        (3, 3, 1, [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], False),        # 4 points, no 3 collinear
        (3, 3, 1, [[1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 2, 0]], True),         # a projective line of PG(2,3)
        (5, 2, 1, [[1, 0], [0, 1], [1, 1], [1, 2], [1, 3], [1, 4]], True)]:    # all points of PG(1,5)
    sp = Space(p, m)
    fac = [np.array([v]) for v in pts]
    M = main_effect_set(sp, fac)
    theta = (p ** (h + 1) - 1) // (p - 1)
    check(M.sum() // (p - 1) == theta, "Theorem 1(c) setup")
    A = admissible_directions(sp, fac, h)
    check((len(A) == 0) == is_subspace, f"Bose-Burton extremal case p={p} m={m}")

# Owen-randomised OA(9,3,3,2) (lattice sampling) == equal-weight coset design (random planes), equivalence in Wynn's sense
p_, m_ = 3, 3
sp = Space(p_, m_)
base = np.array([[i, j, (i + j) % 3] for i in range(3) for j in range(3)])
perms = list(itertools.permutations(range(3)))
P = np.zeros((sp.N, sp.N)); cnt = 0
for p1 in perms:
    for p2 in perms:
        for p3 in perms:
            A = np.stack([np.array(p1)[base[:, 0]], np.array(p2)[base[:, 1]], np.array(p3)[base[:, 2]]], axis=1)
            ind = np.zeros(sp.N); ind[sp.code(A)] = 1
            P += np.outer(ind, ind); cnt += 1
P /= cnt
fac = [np.eye(3, dtype=int)[[i]] for i in range(3)] + [np.eye(3, dtype=int)[[i, j]] for i, j in itertools.combinations(range(3), 2)]
d = CosetDesign(sp, admissible_directions(sp, fac, 2))
pz = d.pair_probability()
Qm = np.array([[pz[sp.code(sp.vecs[y] - sp.vecs[x])] if x != y else d.k / sp.N for y in range(sp.N)] for x in range(sp.N)])
check(np.allclose(P, Qm), "lattice sampling OA(9,3,3,2) == coset design")
print("Lattice sampling with a randomised OA(9,3,3,2) == equal-weight coset design: done")

# eigenvalues of Delta = (pi_xy - pi_x pi_y) are (k^2/N) w(u); trace = N pi (1 - pi)  (link to the trace identity / eigenvalue bound)
for (p_, m_, h_, r_) in [(2, 4, 2, 1), (3, 3, 1, 2), (5, 2, 1, 2)]:
    sp = Space(p_, m_)
    fac = [np.array([[1] + [0] * (m_ - 1)]), np.array([[0, 1] + [0] * (m_ - 2)])]
    A = admissible_directions(sp, fac, h_)
    d = CosetDesign(sp, A, rng.random(len(A)), r=r_)
    Pm = brute_force_pi2(d); pi1 = d.k / sp.N
    Delta = Pm - pi1 ** 2; np.fill_diagonal(Delta, pi1 * (1 - pi1))
    ev = np.sort(np.linalg.eigvalsh(Delta))
    pred = np.sort(np.r_[0.0, (d.k ** 2 / sp.N) * d.weights()[1:]])
    check(np.allclose(ev, pred), f"eigenvalues of Delta p={p_} m={m_}")
print("Eigenvalues of the design variance matrix equal (k^2/N) w(u): done")

# Proposition 3 (main effects bounded by C sigma^2): for translation-invariant designs the worst case over
# Y_C is sigma^2 (C max_M w + max_R w), attained by y* = c1 Re chi_u1 + c2 Re chi_u2; for every design some
# translate of y* has mean squared error at least sigma^2 {C t_M/d + (N/k - 1 - t_M)/|R|}.
def mse_of(design, y):
    ybar = y.mean()
    return sum(pr * (y[s].mean() - ybar) ** 2 for pr, s in design)


def worst_response(sp, M, R, w, C):
    uM = np.nonzero(M)[0][np.argmax(w[M])]
    uR = np.nonzero(R)[0][np.argmax(w[R])]
    y = np.zeros(sp.N)
    for u, target in ((uM, C), (uR, 1.0)):
        part = np.cos(2 * np.pi * (sp.vecs @ sp.vecs[u] % sp.p) / sp.p)
        pw = sp.fourier_power(part)
        y += part * np.sqrt(target / pw[1:].sum())
    return y


for trial in range(60):
    p_, m_ = [(2, 3), (2, 4), (3, 2), (3, 3), (5, 2)][trial % 5]
    sp = Space(p_, m_)
    fac = [np.array([[1] + [0] * (m_ - 1)]), np.array([[0] * (m_ - 1) + [1]])][: 1 + trial % 2]
    M = main_effect_set(sp, fac)
    R = ~M
    R[0] = False
    d_, nR = int(M.sum()), int(R.sum())
    k = p_ ** (1 + trial % (m_ - 1)) if m_ > 2 else p_
    T = sp.N / k - 1
    # a design that is not translation invariant and not balanced: a few random samples of size k
    supp = [rng.choice(sp.N, size=k, replace=False) for _ in range(1 + trial % 4)]
    pr_ = rng.random(len(supp))
    design = list(zip(pr_ / pr_.sum(), supp))
    w = np.zeros(sp.N)
    for pr, s in design:
        eta = np.abs(np.exp(2j * np.pi * (sp.vecs @ sp.vecs[s].T % p_) / p_).sum(axis=1) / k) ** 2
        w += pr * eta
    tM = w[M].sum()
    check(abs(w[1:].sum() - T) < 1e-9, f"trace identity for an arbitrary design p={p_} m={m_}")
    translated = [(pr, (s_vec + a) % p_) for pr, s_vec in [(pr, sp.vecs[s]) for pr, s in design]
                  for a in sp.vecs]
    tinv = [(pr / sp.N, np.array([sp.code(x) for x in pts])) for pr, pts in translated]
    for C in (0.25, 1.0, 4.0):
        y = worst_response(sp, M, R, w, C)
        check(abs(sum(sp.fourier_power(y)[M]) - C) < 1e-9 and abs(sum(sp.fourier_power(y)[R]) - 1) < 1e-9,
              f"worst response lies in Y_C p={p_} m={m_}")
        sup_ti = C * w[M].max() + w[R].max()
        check(abs(mse_of(tinv, y) - sup_ti) < 1e-9, f"Proposition 3(b): worst case of the translated design p={p_} m={m_}")
        for _ in range(20):          # random responses in Y_C never exceed it
            Z = np.fft.fftn(rng.normal(size=(p_,) * m_)).reshape(-1)
            zM = np.real(np.fft.ifftn(np.where(M, Z, 0).reshape((p_,) * m_))).reshape(-1)   # main-effect part
            zR = np.real(np.fft.ifftn(np.where(R, Z, 0).reshape((p_,) * m_))).reshape(-1)   # residual part
            yz = zM * np.sqrt(C * rng.random() / sp.fourier_power(zM).sum()) \
                + zR * np.sqrt(rng.random() / sp.fourier_power(zR).sum())
            check(mse_of(tinv, yz) <= sup_ti + 1e-9, f"Proposition 3(b): no response in Y_C exceeds it p={p_} m={m_}")
        lower = C * tM / d_ + (T - tM) / nR
        best = max(mse_of(design, np.array([y[sp.code((x + a) % p_)] for x in sp.vecs])) for a in sp.vecs)
        check(best >= lower - 1e-9, f"Proposition 3(a): lower bound for an arbitrary design p={p_} m={m_}")
print("Proposition 3 (main effects bounded by C sigma^2): done")

print("\nchecks executed per statement:")
for k_, v_ in counts.items():
    print(f"  {v_:6d}  {k_}")
print("TOTAL FAILURES:", fails)
