"""Is the bound (4) attained by some exactly balanced design where coset mixtures fail?  (Section 6, Appendix A)

For every configuration in which the linear program (8) over coset mixtures stays above the bound (twelve
configurations; results/minimax_lp.csv), all samples of size k balanced on S are enumerated, up to
translation (samples containing 0).  Averaging over translations, a minimax exactly balanced design may be
taken translation invariant; if every balanced sample is a flat, such a design is a coset design
C(h, nu, 1), so the dual certificates of certify_lp.py bound every exactly balanced design.

The only configuration with balanced samples that are not flats is 2^7 balanced on all one- and two-variable
margins with k = 16 (orthogonal arrays OA(16,7,2,2) with distinct runs).  For it, c/enum_oa16.c lists the
weight profiles S_j = sum_{wt(u)=j} J_u^2 of all such arrays (J_u the J-characteristics), and the program
min_nu max_j sum_c nu_c S_j(c) / (k^2 C(7,j)) over j = 3..7 is solved; it is the minimax problem over all
translation-invariant exactly balanced designs after averaging over coordinate permutations, which preserve
S and map each weight class of frequencies onto itself.  Its optimum is verified in rational arithmetic by a
primal mixture and a dual solution.  Writes results/attainment.json.
"""
import itertools, json, math, os, subprocess, tempfile
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

OUT = os.path.dirname(os.path.abspath(__file__))


def balanced_samples(q, m, k, strength):
    """All samples of size k of F_q^m containing 0 and balanced on every margin of at most `strength`
    variables (backtracking over points in increasing code order)."""
    N = q ** m
    dig = np.array(list(itertools.product(range(q), repeat=m)))          # code -> digits
    pairs = list(itertools.combinations(range(m), 2)) if strength >= 2 else []
    cap1, cap2 = k // q, k // q ** 2
    c1 = np.zeros((m, q), dtype=int)
    c2 = np.zeros((len(pairs), q, q), dtype=int)
    sols, cur = [], [0]

    def add(x, s):
        d = dig[x]
        c1[np.arange(m), d] += s
        ok = (c1[np.arange(m), d] <= cap1).all()
        for t, (i, j) in enumerate(pairs):
            c2[t, d[i], d[j]] += s
            ok = ok and c2[t, d[i], d[j]] <= cap2
        return ok

    def rec(last):
        if len(cur) == k:
            sols.append(list(cur))
            return
        for x in range(last + 1, N - (k - len(cur)) + 1):
            if add(x, 1):
                cur.append(x)
                rec(x)
                cur.pop()
            add(x, -1)
    add(0, 1)
    rec(0)
    return dig, sols


def is_flat(dig, q, s):
    """s contains 0: it is a flat iff it is a linear subspace (closed under addition and scaling)."""
    pts = {tuple(dig[x]) for x in s}
    for a, b in itertools.combinations(pts, 2):
        for c in range(1, q):
            if tuple((np.array(a) + c * np.array(b)) % q) not in pts:
                return False
    return True


def oa16_profiles():
    exe = os.path.join(tempfile.mkdtemp(), "enum_oa16")
    subprocess.run(["gcc", "-O2", "-o", exe, f"{OUT}/c/enum_oa16.c"], check=True)
    prof = {}
    # t = 0 (third column = product of the first two) is not needed: for any two columns at most one other
    # column is plus or minus their product (distinct columns are orthogonal), so every array has a triple
    # with t in {1, 2, 3}, and t = 3 becomes t = 1 when the third column changes sign
    for t in (1, 2):
        res = subprocess.run([exe, str(t)], check=True, capture_output=True, text=True).stdout.splitlines()
        for line in res[1:]:
            f = list(map(int, line.split()))
            key = (f[1],) + tuple(f[2:])
            prof[key] = prof.get(key, 0) + f[0]
    return prof


def exact_lp(A, bound):
    """min z s.t. A nu <= z, sum nu = 1 (A: classes x columns, Fractions), solved in floats and certified in
    rational arithmetic: a rounded primal mixture nu gives the upper bound max_j (A nu)_j, a rounded dual
    distribution beta on the classes gives the lower bound min_c sum_j beta_j A_jc; they must coincide."""
    Af = np.array([[float(x) for x in row] for row in A])
    nj, nc = Af.shape
    res = linprog(np.r_[np.zeros(nc), 1.0], A_ub=np.c_[Af, -np.ones(nj)], b_ub=np.zeros(nj),
                  A_eq=np.r_[np.ones(nc), 0.0].reshape(1, -1), b_eq=[1.0],
                  bounds=[(0, None)] * nc + [(None, None)], method="highs")
    for den in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 8, 10 ** 10):
        nu = [F(max(x, 0.0)).limit_denominator(den) for x in res.x[:nc]]
        nu = [x / sum(nu) for x in nu]
        upper = max(sum(A[j][c] * nu[c] for c in range(nc)) for j in range(nj))
        beta = [F(max(-x, 0.0)).limit_denominator(den) for x in res.ineqlin.marginals]
        beta = [x / sum(beta) for x in beta]
        lower = min(sum(beta[j] * A[j][c] for j in range(nj)) for c in range(nc))
        if upper == lower:
            break
    assert upper == lower, (float(upper), float(lower))
    return upper, {"z_float": float(res.x[-1]), "nu": {str(c): str(x) for c, x in enumerate(nu) if x},
                   "beta": [str(x) for x in beta], "ratio_to_bound": float(upper / bound)}


if __name__ == "__main__":
    configs = [(2, 5, 2, 1), (2, 5, 4, 1), (2, 6, 2, 1), (2, 6, 4, 1), (2, 7, 2, 1), (2, 7, 4, 1),
               (2, 6, 8, 2), (2, 7, 8, 2), (3, 4, 3, 1), (3, 5, 3, 1), (3, 4, 9, 2)]
    out = []
    for q, m, k, st in configs:
        dig, sols = balanced_samples(q, m, k, st)
        nflat = sum(is_flat(dig, q, s) for s in sols)
        out.append({"q": q, "m": m, "k": k, "margins": "one-way" if st == 1 else "one- and two-way",
                    "balanced_samples_containing_0": len(sols), "flats": nflat,
                    "balanced_samples": len(sols) * q ** m // k, "all_flats": nflat == len(sols)})
        print(out[-1], flush=True)
    # 2^7, one- and two-variable margins, k = 16
    prof = oa16_profiles()
    n = [math.comb(7, j) for j in range(8)]
    keys = sorted(prof)
    A = [[F(key[1 + j], 256 * n[j]) for key in keys] for j in range(3, 8)]
    bound = (F(128, 16) - 1) / 99                                         # (N/k - 1)/|R_S|, |R_S| = 99
    assert bound == F(7, 99)
    flat_cols = [i for i, key in enumerate(keys) if key[0] == 1]
    z_all, info_all = exact_lp(A, bound)
    z_flat, info_flat = exact_lp([[row[i] for i in flat_cols] for row in A], bound)
    out.append({"q": 2, "m": 7, "k": 16, "margins": "one- and two-way", "distinct_profiles": len(keys),
                "profiles_of_flats": len(flat_cols), "arrays_counted": sum(prof.values()),
                "bound": str(bound), "minimax_flats": str(z_flat), "minimax_flats_ratio": float(z_flat / bound),
                "minimax_all_designs": str(z_all), "minimax_all_ratio": float(z_all / bound),
                "lp_all": info_all, "lp_flats": info_flat,
                "profiles": [{"flat": key[0], "S": list(key[1:]), "count": prof[key]} for key in keys]})
    print({k_: v for k_, v in out[-1].items() if k_ not in ("profiles", "lp_all", "lp_flats")})
    json.dump(out, open(f"{OUT}/results/attainment.json", "w"), indent=1)
