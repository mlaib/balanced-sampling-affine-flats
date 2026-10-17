"""Independent check of the enumeration of OA(16,7,2,2) in settle_attainment.py (Appendix A).

c/enum_oa16.c lists the balanced samples of size 16 of F_2^m (one- and two-variable margins) in a canonical
form: columns 1 and 2 fixed, column 3 with t, 4-t, 4-t and t entries +1 in the four blocks defined by
columns 1 and 2, its +1 entries first in each block, and columns 4 to m an increasing tuple of pairwise
orthogonal candidates.  c/enum_direct.c lists every balanced sample containing 0 by backtracking over the
units, without any reduction by symmetries.  The two are linked exactly.  A labelled array whose first three
columns are in canonical form with parameter t is a balanced sample together with an ordering of its runs;
a sample whose third column has parameter t has (t!(4-t)!)^4 such orderings, and each array listed by
enum_oa16 stands for the (m-3)! orders of its columns 4 to m.  So the number of balanced samples with
profile P and parameter t is n_t(P) = (m-3)! C_t(P) / (t!(4-t)!)^4, where C_t(P) is the count of enum_oa16.
A translation leaves P and tau = min(t, 4-t) unchanged, and every class of samples closed under
translations has 16/2^m of its members containing 0, so the count D_tau(P) of enum_direct must satisfy
    2^m D_tau(P) / 16 = sum_{t: min(t, 4-t) = tau} n_t(P).
This is checked for every profile and tau, for m = 5 and 6 with t = 0, ..., 4 and for m = 7 with t = 0, ..., 3, where
t = 4 follows from t = 0 by a change of sign of column 3 (as the runs for m = 5 and 6 confirm); the set of profiles
found with t in {1, 2}, the only values used by settle_attainment.py, must be the set found by the direct
enumeration.  About half an hour; writes results/crosscheck_oa16.json.
"""
import json, os, subprocess, tempfile
from fractions import Fraction as F
from math import factorial

OUT = os.path.dirname(os.path.abspath(__file__))


def run(exe, *args):
    out = subprocess.run([exe, *map(str, args)], check=True, capture_output=True, text=True).stdout.splitlines()
    head = dict(x.split("=") for x in out[0][2:].split())
    rows = [list(map(int, line.split())) for line in out[1:]]
    return head, rows


if __name__ == "__main__":
    tmp = tempfile.mkdtemp()
    exe = {}
    for name in ("enum_oa16", "enum_direct"):
        exe[name] = os.path.join(tmp, name)
        subprocess.run(["gcc", "-O2", "-o", exe[name], f"{OUT}/c/{name}.c"], check=True)
    pi = {t: (factorial(t) * factorial(4 - t)) ** 4 for t in range(5)}
    report = []
    for m in (5, 6, 7):
        ts = (0, 1, 2, 3, 4) if m < 7 else (0, 1, 2, 3)
        C, canon_counts = {}, {}
        for t in ts:
            head, rows = run(exe["enum_oa16"], t, m)
            C[t] = {(r[1], tuple(r[2:])): r[0] for r in rows}
            canon_counts[t] = {"tuples": int(head["tuples"]), "with_distinct_runs": int(head["with_distinct_runs"]),
                               "profiles": int(head["distinct_profiles"])}
        if 4 not in C and 0 in C:
            C[4] = C[0]                                    # a sign change of column 3 maps t to 4 - t
        head, rows = run(exe["enum_direct"], m)
        D = {(r[1], r[2], tuple(r[3:])): r[0] for r in rows}
        assert sum(D.values()) == int(head["samples_containing_0"])
        profiles_direct = {(fl, S) for fl, tau, S in D}
        profiles_canon = set(C[1]) | set(C[2])
        ok_sets = profiles_direct == profiles_canon
        ok_sym = C[3] == C[1]
        checked, ok_counts = [], True
        for tau, group in ((0, (0, 4)), (1, (1, 3)), (2, (2,))):
            if not all(t in C for t in group):
                continue
            for P in profiles_direct | profiles_canon:
                lhs = F(2 ** m, 16) * D.get((P[0], tau, P[1]), 0)
                rhs = sum(F(C[t].get(P, 0) * factorial(m - 3), pi[t]) for t in group)
                ok_counts &= lhs == rhs
            checked.append(tau)
        nflat0 = sum(v for (fl, tau, S), v in D.items() if fl)
        rec = {"m": m, "t_enumerated": list(ts), "canonical": {str(t): canon_counts[t] for t in ts},
               "direct_samples_containing_0": sum(D.values()),
               "balanced_samples": sum(D.values()) * 2 ** m // 16,
               "flats_containing_0": nflat0,
               "direct_by_tau": {str(tau): sum(v for (fl, t, S), v in D.items() if t == tau) for tau in range(3)},
               "flat_profile_pairs": len(profiles_direct), "distinct_profiles": len({S for fl, S in profiles_direct}),
               "profiles_t12_equal_direct": ok_sets, "counts_t3_equal_t1": ok_sym,
               "count_identity_checked_for_tau": checked, "count_identity_holds": ok_counts}
        if m == 7:
            att = json.load(open(f"{OUT}/results/attainment.json"))[-1]
            stored = {(p["flat"], tuple(p["S"])): p["count"] for p in att["profiles"]}
            rec["canonical_t12_equal_stored"] = stored == {P: C[1].get(P, 0) + C[2].get(P, 0) for P in profiles_canon}
        report.append(rec)
        print(rec, flush=True)
        assert ok_sets and ok_sym and ok_counts and rec.get("canonical_t12_equal_stored", True)
    json.dump(report, open(f"{OUT}/results/crosscheck_oa16.json", "w"), indent=1)
    print("all cross-checks passed")
