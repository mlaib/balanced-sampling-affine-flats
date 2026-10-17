"""Standalone check of the exact certificates of Table 1, Section 6 and Appendix A.

It uses the Python standard library only, and none of the solver, the enumeration programs or the other
modules of the repository; it reads results/lp_certificates.json and results/attainment.json.

1. Program (8), for the 31 configurations of Table 1 and the archive.  The admissible directions are listed
   afresh as their annihilators D = V^perp: the subspaces of F_q^m of dimension m - h, by their reduced
   row-echelon bases, whose nonzero vectors all have Hamming weight above the order s of the margins
   (s = 1 for one-variable margins, s = 2 with two-variable margins).  For a mixture nu that is uniform within
   classes of directions with the same weight distribution, c_nu(u) at a frequency u of weight j is
   sum_classes P(class) A_j(D)/n_j, n_j = C(m,j)(q-1)^j.  Where (8) reaches the bound, this is checked to
   equal (Q-1)/|R_S| for every weight j > s.  Where it stays above, the stored primal mixture is checked to
   give max_j c_nu = zeta, and the stored dual distribution beta, constant on the weight classes of R_S, to
   give beta(D minus 0) >= zeta for every admissible D, one by one; so zeta is the minimax value over all
   coset designs C(h, nu, 1).  The weight distributions and class sizes of the stored mixtures must occur
   among the directions found.
2. Attainment over all designs.  For the eleven configurations in which every balanced sample containing 0
   is a flat, that is a linear subspace V, their number must equal the number of admissible directions.  For
   OA(16,7,2,2): every stored profile (S_0, ..., S_7) must satisfy S_0 = k^2, S_1 = S_2 = 0 and Parseval's
   identity sum_j S_j = N k, and its inverse MacWilliams transform A_i = 2^-m sum_j S_j K_i(j), with K_i the
   Krawtchouk polynomials, must count ordered pairs of units at distance i: k pairs at distance 0, even
   nonnegative integers otherwise.  The profiles of flats must be exactly the vectors k^2 A_j(D) of the
   admissible D.  The stored primal mixture and dual distribution on the weight classes 3 to 7 must give
   the same value over all profiles, and over the profiles of flats.
Exits with status 1 if any check fails.
"""
import itertools, json, math, os, sys
from fractions import Fraction as F

OUT = os.path.dirname(os.path.abspath(__file__))
failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)
        print("FAILED:", msg)


def subspaces(q, m, r):
    """All r-dimensional subspaces of F_q^m (q prime), each as the set of its vectors."""
    for piv in itertools.combinations(range(m), r):
        free = [(i, c) for i in range(r) for c in range(piv[i] + 1, m) if c not in piv]
        for vals in itertools.product(range(q), repeat=len(free)):
            B = [[0] * m for _ in range(r)]
            for i in range(r):
                B[i][piv[i]] = 1
            for (i, c), v in zip(free, vals):
                B[i][c] = v
            yield {tuple(sum(a * b[c] for a, b in zip(coef, B)) % q for c in range(m))
                   for coef in itertools.product(range(q), repeat=r)}


def weight_distribution(D, m):
    A = [0] * (m + 1)
    for u in D:
        A[sum(1 for x in u if x)] += 1
    return tuple(A)


def admissible(q, m, k, s):
    """Weight distributions of the annihilators of the admissible directions of dimension log_q k."""
    h = round(math.log(k, q))
    assert q ** h == k
    out = []
    for D in subspaces(q, m, m - h):
        A = weight_distribution(D, m)
        if all(A[j] == 0 for j in range(1, s + 1)):
            out.append(A)
    return out


def residual_weights(q, m, s):
    return [j for j in range(s + 1, m + 1)], [math.comb(m, j) * (q - 1) ** j for j in range(m + 1)]


def check_lp(cert, attained):
    q, m, k = cert["q"], cert["m"], cert["k"]
    s = 1 if cert["balanced_on"] == "one-way" else 2
    dirs = admissible(q, m, k, s)
    js, n = residual_weights(q, m, s)
    RS = sum(n[j] for j in js)
    Q = F(q ** m, k)
    bound = (Q - 1) / RS
    name = f"q={q} m={m} k={k} {cert['balanced_on']}"
    check(F(cert["bound"]) == bound, f"{name}: bound")
    mix = cert["mixture"] if attained else cert["primal_mixture"]
    P = [F(c["probability_of_class"]) for c in mix]
    check(all(p >= 0 for p in P) and sum(P) == 1, f"{name}: mixture is a probability distribution")
    for c in mix:
        W = tuple(c["weight_distribution"])
        check(dirs.count(W) == c["n_directions_in_class"], f"{name}: class {W} has {dirs.count(W)} directions")
    cnu = {j: sum(p * F(c["weight_distribution"][j], n[j]) for p, c in zip(P, mix)) for j in js}
    if attained:
        check(all(v == bound for v in cnu.values()), f"{name}: c_nu equals the bound on R_S")
        return name, len(dirs), bound, bound
    zeta = F(cert["minimax_value"])
    check(max(cnu.values()) == zeta, f"{name}: primal mixture gives max c_nu = zeta")
    mass = {int(j): F(v) for j, v in cert["dual_mass_by_hamming_weight"].items()}
    check(sorted(mass) == js and all(v >= 0 for v in mass.values()) and sum(mass.values()) == 1,
          f"{name}: dual is a probability distribution on R_S")
    worst = min(sum(mass[j] * F(A[j], n[j]) for j in js) for A in dirs)
    check(worst >= zeta, f"{name}: beta(D minus 0) >= zeta for every admissible D")
    check(cert["n_admissible"] == len(dirs), f"{name}: number of admissible directions")
    return name, len(dirs), bound, zeta


def krawtchouk(m, i, j):
    """K_i(j) for q = 2."""
    return sum((-1) ** l * math.comb(j, l) * math.comb(m - j, i - l) for l in range(i + 1))


def check_oa16(entry, cert):
    q, m, k = 2, 7, 16
    N = q ** m
    js, n = residual_weights(q, m, 2)
    flats = {tuple(k * k * a for a in A) for A in admissible(q, m, k, 2)}
    profiles = [(p["flat"], tuple(p["S"])) for p in entry["profiles"]]
    check(len(set(profiles)) == len(profiles) == entry["distinct_profiles"], "OA16: profiles are distinct")
    for flat, S in profiles:
        check(S[0] == k * k and S[1] == S[2] == 0 and sum(S) == N * k, f"OA16: identities for {S}")
        A = [F(sum(S[j] * krawtchouk(m, i, j) for j in range(m + 1)), N) for i in range(m + 1)]
        check(A[0] == k and sum(A) == k * k and all(a.denominator == 1 and a >= 0 and a % 2 == 0 for a in A[1:]),
              f"OA16: distance distribution of {S}")
        if flat:
            check(S in flats, f"OA16: flat profile {S} is k^2 A(D) for an admissible D")
    check({S for flat, S in profiles if flat} == flats, "OA16: every admissible D gives a flat profile")
    check(sum(p["count"] for p in entry["profiles"]) == entry["arrays_counted"], "OA16: counts")
    bound = F(N, k) - 1
    bound /= sum(n[j] for j in js)
    check(F(entry["bound"]) == bound, "OA16: bound")
    vecs = {"all": [S for _, S in profiles], "flats": [S for flat, S in profiles if flat]}
    keys = sorted(profiles)
    flat_keys = [key for key in keys if key[0] == 1]
    values = {}
    for which, lp, cols in (("all", entry["lp_all"], keys), ("flats", entry["lp_flats"], flat_keys)):
        a = lambda S, j: F(S[j], k * k * n[j])
        nu = {cols[int(i)][1]: F(p) for i, p in lp["nu"].items()}
        check(all(p >= 0 for p in nu.values()) and sum(nu.values()) == 1, f"OA16 {which}: primal is a distribution")
        upper = max(sum(p * a(S, j) for S, p in nu.items()) for j in js)
        beta = [F(b) for b in lp["beta"]]
        check(all(b >= 0 for b in beta) and sum(beta) == 1 and len(beta) == len(js), f"OA16 {which}: dual is a distribution")
        lower = min(sum(b * a(S, j) for b, j in zip(beta, js)) for S in vecs[which])
        check(upper == lower, f"OA16 {which}: primal value {upper} equals dual value {lower}")
        values[which] = upper
    check(values["all"] == F(entry["minimax_all_designs"]) and values["flats"] == F(entry["minimax_flats"]),
          "OA16: stored minimax values")
    check(values["flats"] == F(cert["minimax_value"]), "OA16: flats give the minimax value of program (8)")
    return len(set(S for _, S in profiles)), len(flats), bound, values


if __name__ == "__main__":
    lp = json.load(open(f"{OUT}/results/lp_certificates.json"))
    n_dirs = {}
    for cert in lp["certificates"]:
        name, nd, bound, value = check_lp(cert, True)
        n_dirs[name] = nd
        print(f"{name}: {nd} admissible directions, bound {bound} attained")
    for cert in lp["nonattainment_certificates"]:
        name, nd, bound, value = check_lp(cert, False)
        n_dirs[name] = nd
        print(f"{name}: {nd} admissible directions, minimax {value} = {value / bound} x bound {bound}")
    att = json.load(open(f"{OUT}/results/attainment.json"))
    for e in att:
        name = f"q={e['q']} m={e['m']} k={e['k']} {e['margins']}"
        if "profiles" in e:
            cert = next(c for c in lp["nonattainment_certificates"]
                        if (c["q"], c["m"], c["k"], c["balanced_on"]) == (2, 7, 16, "one- and two-way"))
            nvec, nflat, bound, values = check_oa16(e, cert)
            print(f"{name}: {len(e['profiles'])} (flat, profile) pairs, {nvec} distinct profiles, {nflat} from flats; "
                  f"minimax {values['all']} over all profiles and {values['flats']} over flats, bound {bound}")
        else:
            check(e["all_flats"] and e["flats"] == e["balanced_samples_containing_0"] == n_dirs[name],
                  f"{name}: balanced samples containing 0 are the {n_dirs[name]} admissible subspaces")
            print(f"{name}: {e['balanced_samples_containing_0']} balanced samples containing 0, all subspaces")
    print("all checks passed" if not failures else f"{len(failures)} checks failed")
    sys.exit(1 if failures else 0)
