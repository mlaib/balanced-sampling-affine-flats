"""Audit of the cube method balanced on rows and columns of the 61 x 61 grid (Section 7.2, Appendix A).

A Python mirror of `cube` in BalancedSampling 2.1.1 (src/CubeClass.cc): units in random order; the flight
phase moves p+1 fractional units at a time along an RREF null vector and stops as soon as at most p units
are fractional; the landing phase then uses, for the n remaining units, only the first n-1 balancing
variables.  In the flight phase the null vectors are computed from the bipartite-graph structure of the
balancing matrix (graph_null_vector), which gives the RREF vectors exactly and is about ten times faster;
check_graph_null_vector compares the two step by step. For each draw the flight is run once, and from its end point two landings are applied:

  bs        the landing of BalancedSampling (suppression of the last balancing variables);
  continued the flight continued while a nonzero u with A_J' u = 0 exists (all p equations kept),
            followed by the same landing for any units still fractional.

With the balancing variables of rows and columns the constraint polytope is integral, so the continued
flight ends at an exactly balanced sample. Writes results/cube_audit_S2.csv (summary) and
results/cube_audit_S2_draws.csv (one line per draw and landing, with the selected units), and
results/cube_audit_S2_decomposition.csv, which splits each mean squared error into the part due to the
row and column main effects (zero for balanced samples), the residual part and the cross term.
Usage: python3 audit_cube.py [B] [processes]      python3 audit_cube.py --decompose (from stored draws)
"""
import csv, os, sys, time
from multiprocessing import Pool
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
q = 61
EPS = 1e-12      # integrality tolerance of BalancedSampling (snaps units during flight and landing)
TOL = 1e-9       # units within TOL of 0 or 1 count as integral in the audit and in the continuation


def rref(M):
    """In-place reduced row echelon form with exact zero tests, as ReducedRowEchelonForm.cc.  This copies the
    package's defect too: when the elimination is not exact in floating point (three or more line stratifications),
    rounding errors are taken as pivots and the direction can violate the balancing equations, the sample size
    included (check_rref_defect.py; cube_grid_fixed.py gives numerically sound directions)."""
    nrows, ncols = M.shape
    lead = 0
    for r in range(nrows):
        if ncols <= lead:
            return
        while True:
            nz = np.flatnonzero(M[r:, lead] != 0.0)
            if nz.size:
                i = r + nz[0]
                break
            lead += 1
            if lead == ncols:
                return
        if i != r:
            M[[i, r]] = M[[r, i]]
        br = M[r]
        br[lead + 1:] /= br[lead]
        br[lead] = 1.0
        t = M[:, lead].copy()
        t[r] = 0.0
        M[:, lead] = 0.0
        M[r, lead] = 1.0
        rows = np.flatnonzero(t)          # rows with t = 0 are unchanged (x - 0*b = x)
        if rows.size:
            M[rows, lead + 1:] -= np.outer(t[rows], br[lead + 1:])
        lead += 1


def null_vector(M):
    """Null vector of an RREF matrix with ncols = nrows + 1, as CubeVectorInNullSpace.cc."""
    nrows, ncols = M.shape
    if M[nrows - 1, nrows - 1] == 1.0:
        u = np.empty(ncols)
        u[-1] = 1.0
        u[:nrows] = -M[:, ncols - 1]
        return u
    u = np.where(np.arange(ncols) % 2 == 0, -1.0, 1.0)
    for i in range(nrows):
        ones = np.flatnonzero(M[i] == 1.0)
        if not ones.size:
            continue
        lead = ones[0]
        u[lead] = -np.dot(u[lead + 1:], M[i, lead + 1:])
    return u


def graph_null_vector(units):
    """The null vector that null_vector(rref(B)) returns in exact arithmetic, for B the balancing matrix
    [1, rows 1-60, columns 1-60]/pi restricted to `units` (in candidate order).  The columns of B have the
    matroid of the bipartite graph of rows and columns (unit (i, j) is an edge), so the pivot columns of
    the RREF form the spanning forest built greedily in candidate order; the free columns get -1/+1 by
    the parity of their position, and the forest edges then solve the zero-sum equation at every row and
    column vertex.  Used by the flight phase for speed; check_graph_null_vector compares the two."""
    i, j = np.divmod(np.asarray(units), q)
    n = len(units)
    parent = list(range(2 * q))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree = [False] * n
    for e in range(n):
        a, b = find(int(i[e])), find(q + int(j[e]))
        if a != b:
            parent[a] = b
            tree[e] = True
    u = np.where(np.arange(n) % 2 == 0, -1.0, 1.0)
    d = np.zeros(2 * q)
    adj = [[] for _ in range(2 * q)]
    for e in range(n):
        if tree[e]:
            adj[i[e]].append(e)
            adj[q + j[e]].append(e)
        else:
            d[i[e]] += u[e]
            d[q + j[e]] += u[e]
    deg = [len(a_) for a_ in adj]
    done = [False] * n
    stack = [v for v in range(2 * q) if deg[v] == 1]
    while stack:
        v = stack.pop()
        if deg[v] != 1:
            continue
        e = next(e_ for e_ in adj[v] if not done[e_])
        u[e] = -d[v]
        done[e] = True
        w = q + int(j[e]) if v == i[e] else int(i[e])
        d[w] += u[e]
        d[v] = 0.0
        deg[v] -= 1
        deg[w] -= 1
        if deg[w] == 1:
            stack.append(w)
    assert np.abs(d).max() < 1e-9 and all(done[e] for e in range(n) if tree[e])
    return u


class UniformStream:
    """Supplies given uniform numbers in order, standing in for R's unif_rand (validate_transcription.py)."""
    def __init__(self, U):
        self.U, self.i = U, 0

    def __call__(self):
        u = self.U[self.i]
        self.i += 1
        return u


class Cube:
    def __init__(self, prob, x, rng, stream=None):
        self.N, self.p = x.shape
        self.A = x / prob[:, None]
        self.pi = prob.astype(float).copy()
        self.rng, self.stream = rng, stream
        if stream is None:
            self.lst = list(rng.permutation(self.N))
        else:                                   # IndexList::Shuffle: k = i + floor(u (len - i))
            self.lst = list(range(self.N))
            for i in range(self.N - 1):
                k = i + int(stream() * (self.N - i))
                self.lst[i], self.lst[k] = self.lst[k], self.lst[i]
        self.pos = {u: i for i, u in enumerate(self.lst)}
        self.sample = []

    def unif(self, scale=1.0):
        return self.stream() * scale if self.stream is not None else self.rng.uniform(0.0, scale)

    def erase(self, unit):
        k = self.pos.pop(unit)
        last = self.lst.pop()
        if k < len(self.lst):
            self.lst[k] = last
            self.pos[last] = k

    def update(self, cand, u):
        pi = self.pi[cand]
        with np.errstate(divide="ignore", invalid="ignore"):
            l_down, l_up = np.abs(pi / u), np.abs((1.0 - pi) / u)
        pos = u >= 0.0
        lam1 = min(np.min(l_up[pos], initial=np.inf), np.min(l_down[~pos], initial=np.inf))
        lam2 = min(np.min(l_down[pos], initial=np.inf), np.min(l_up[~pos], initial=np.inf))
        lam = lam1 if self.unif(lam1 + lam2) < lam2 else -lam2
        self.pi[cand] = pi + lam * u
        for unit in cand:
            v = self.pi[unit]
            if v <= EPS or v >= 1.0 - EPS:
                self.erase(unit)
                if v >= 1.0 - EPS:
                    self.sample.append(unit)

    def flight(self, graph=True):
        """Fast flight phase of BalancedSampling: p+1 fractional units at a time, until at most p remain.
        graph=True computes the null vectors with graph_null_vector, graph=False with the RREF."""
        m = self.p + 1
        while len(self.lst) >= m:
            cand = self.lst[:m]
            if graph:
                u = graph_null_vector(cand)
            else:
                B = self.A[cand].T.copy()
                rref(B)
                u = null_vector(B)
            self.update(cand, u)

    def fractional(self):
        return [u for u in self.lst if TOL < self.pi[u] < 1.0 - TOL]

    def nullity(self):
        F = self.fractional()
        if not F:
            return 0
        s = np.linalg.svd(self.A[F].T, compute_uv=False)
        tol = s.max() * max(self.p, len(F)) * np.finfo(float).eps * 10
        return len(F) - int((s > tol).sum())

    def continue_flight(self):
        """Move the fractional units along null vectors of A_F' (all p equations) until none exists."""
        steps = 0
        while True:
            F = self.fractional()
            if not F:
                break
            BF = self.A[F].T
            _, s, Vt = np.linalg.svd(BF)
            tol = s.max() * max(BF.shape) * np.finfo(float).eps * 10
            if int((s > tol).sum()) == len(F):
                break
            self.update(F, Vt[-1])
            steps += 1
        return steps

    def landing(self):
        while len(self.lst) > 1:
            m = len(self.lst)
            cand = list(self.lst)
            B = self.A[cand][:, :m - 1].T.copy()
            rref(B)
            self.update(cand, null_vector(B))
        if len(self.lst) == 1:
            unit = self.lst[0]
            if self.unif() < self.pi[unit]:
                self.sample.append(unit)
            self.erase(unit)
        return np.array(sorted(self.sample))

    def copy(self):
        c = Cube.__new__(Cube)
        c.N, c.p, c.A, c.rng, c.stream = self.N, self.p, self.A, self.rng, self.stream
        c.pi, c.lst, c.pos, c.sample = self.pi.copy(), list(self.lst), dict(self.pos), list(self.sample)
        return c


def counts(s):
    """Stratum counts on rows, columns, anti-diagonals i+j and diagonals i-j (mod 61)."""
    i, j = np.divmod(s, q)
    return [np.bincount(v % q, minlength=q) for v in (i, j, i + j, i - j)]


_DATA = {}


def data():
    if not _DATA:
        _DATA["X"] = np.loadtxt(f"{OUT}/data/grid_X.csv", delimiter=",")[:, :120]   # row levels 1-60, column levels 1-60
        _DATA["y"] = np.loadtxt(f"{OUT}/data/grid_y.csv")
    return _DATA["X"], _DATA["y"]


def check_graph_null_vector(args):
    """Run one flight with RREF null vectors and compare graph_null_vector with them at every step;
    then rerun the flight with graph_null_vector from the same seed and compare the end points."""
    k, seed = args
    X, y = data()
    prob = np.full(len(y), k / len(y))
    c = Cube(prob, np.c_[prob, X], np.random.default_rng(seed))
    m, worst, steps = c.p + 1, 0.0, 0
    while len(c.lst) >= m:
        cand = c.lst[:m]
        B = c.A[cand].T.copy()
        rref(B)
        u = null_vector(B)
        worst = max(worst, float(np.abs(u - graph_null_vector(cand)).max() / np.abs(u).max()))
        c.update(cand, u)
        steps += 1
    g = Cube(prob, np.c_[prob, X], np.random.default_rng(seed))
    g.flight(graph=True)
    same = sorted(c.sample) == sorted(g.sample) and sorted(c.lst) == sorted(g.lst)
    return {"k": k, "seed": seed, "steps": steps, "max_rel_diff": worst, "same_end_point": same,
            "max_abs_pi_diff": float(np.abs(c.pi - g.pi).max())}


def one_draw(args):
    k, seed = args
    X, y = data()
    N = len(y)
    prob = np.full(N, k / N)
    rng = np.random.default_rng(seed)
    c = Cube(prob, np.c_[prob, X], rng)
    t0 = time.time()
    c.flight()
    J, nul = len(c.fractional()), c.nullity()
    out = []
    for name in ("bs", "continued"):
        d = c.copy()
        steps = d.continue_flight() if name == "continued" else 0
        J_after = len(d.fractional())
        s = d.landing()
        rc, cc, ac, dc = counts(s)
        out.append({"k": k, "seed": seed, "landing": name, "flight_end_fractional": J,
                    "flight_end_nullity": nul, "continuation_steps": steps,
                    "fractional_before_landing": J_after, "size": len(s),
                    "maxdev_rows": float(np.abs(rc - k / q).max()), "maxdev_cols": float(np.abs(cc - k / q).max()),
                    "maxdev_S4": float(max(np.abs(c_ - k / q).max() for c_ in (rc, cc, ac, dc))),
                    "estimate": float(y[s].mean()), "seconds": time.time() - t0,
                    "sample": " ".join(map(str, s))})
    return out


def decompose():
    """MSE = E(mbar^2) + E(rbar^2) + 2E(mbar rbar), with y - Ybar = m + rho, m the row-plus-column main effects."""
    Y = np.loadtxt(f"{OUT}/data/volcano61.csv", delimiter=",")
    y, Ybar = Y.reshape(-1), Y.mean()
    m = ((Y.mean(axis=1, keepdims=True) - Ybar) + (Y.mean(axis=0, keepdims=True) - Ybar)).reshape(-1)
    rho = y - Ybar - m
    rows = list(csv.DictReader(open(f"{OUT}/results/cube_audit_S2_draws.csv")))
    out = []
    for k in (61, 183):
        srs = (1 - k / y.size) * y.var(ddof=1) / k
        for name in ("bs", "continued"):
            S = [np.array(r["sample"].split(), dtype=int) for r in rows if int(r["k"]) == k and r["landing"] == name]
            em = np.array([m[s_].mean() for s_ in S])
            er = np.array([rho[s_].mean() for s_ in S])
            est = np.array([float(r["estimate"]) for r in rows if int(r["k"]) == k and r["landing"] == name])
            # the parts add up to the squared error of every draw, and the total is the reported design effect
            assert np.abs((em + er) ** 2 - (est - Ybar) ** 2).max() < 1e-9
            row = {"k": k, "landing": name, "B": len(S),
                   "bias": float(est.mean() - Ybar), "bias_se": float(est.std(ddof=1) / np.sqrt(len(S))),
                   "max_drawwise_identity_error": float(np.abs(em ** 2 + er ** 2 + 2 * em * er - (est - Ybar) ** 2).max())}
            for part, v in (("main", em ** 2), ("residual", er ** 2), ("cross", 2 * em * er), ("total", (em + er) ** 2)):
                row[f"deff_{part}"] = v.mean() / srs
                row[f"deff_{part}_se"] = v.std(ddof=1) / np.sqrt(len(S)) / srs
            out.append(row)
    with open(f"{OUT}/results/cube_audit_S2_decomposition.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    for r in out:
        print({a: (round(b, 4) if isinstance(b, float) else b) for a, b in r.items()})


if __name__ == "__main__" and "--decompose" in sys.argv:
    decompose()
elif __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    y = np.loadtxt(f"{OUT}/data/grid_y.csv")
    N, Ybar = len(y), y.mean()
    jobs = [(k, int(s)) for k in (61, 183)
            for s in np.random.SeedSequence(2026 + k).generate_state(B)]
    draws, t0 = [], time.time()
    with Pool(procs) as pool:
        for n, res in enumerate(pool.imap(one_draw, jobs, chunksize=4), 1):
            draws += res
            if n % 200 == 0:
                print(f"{n}/{len(jobs)} draws, {time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    with open(f"{OUT}/results/cube_audit_S2_draws.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(draws[0]))
        w.writeheader(); w.writerows(draws)
    summary = []
    for k in (61, 183):
        srs = (1 - k / N) * y.var(ddof=1) / k
        for name in ("bs", "continued"):
            d = [r for r in draws if r["k"] == k and r["landing"] == name]
            se2 = np.array([(r["estimate"] - Ybar) ** 2 for r in d])
            exact_r = np.array([r["maxdev_rows"] < 1e-9 for r in d])
            exact_c = np.array([r["maxdev_cols"] < 1e-9 for r in d])
            summary.append({"k": k, "landing": name, "B": len(d), "deff": se2.mean() / srs,
                            "deff_se": se2.std(ddof=1) / np.sqrt(len(d)) / srs,
                            "share_exact_rows": exact_r.mean(), "share_exact_cols": exact_c.mean(),
                            "share_exact_S2": (exact_r & exact_c).mean(),
                            "share_exact_S4": np.mean([r["maxdev_S4"] < 1e-9 for r in d]),
                            "mean_maxdev_rows": np.mean([r["maxdev_rows"] for r in d]),
                            "mean_maxdev_cols": np.mean([r["maxdev_cols"] for r in d]),
                            "size_min": min(r["size"] for r in d), "size_max": max(r["size"] for r in d),
                            "flight_end_fractional_mean": np.mean([r["flight_end_fractional"] for r in d]),
                            "flight_end_fractional_min": min(r["flight_end_fractional"] for r in d),
                            "flight_end_fractional_max": max(r["flight_end_fractional"] for r in d),
                            "share_nullity_positive": np.mean([r["flight_end_nullity"] > 0 for r in d]),
                            "flight_end_nullity_mean": np.mean([r["flight_end_nullity"] for r in d]),
                            "fractional_before_landing_max": max(r["fractional_before_landing"] for r in d)})
    with open(f"{OUT}/results/cube_audit_S2.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader(); w.writerows(summary)
    for r in summary:
        print({a: (round(b, 4) if isinstance(b, float) else b) for a, b in r.items()})
    decompose()
