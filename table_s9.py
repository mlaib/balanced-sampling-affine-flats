"""Table S9 of the supplementary material: verification map of Proposition 6 (Appendix A).

For the twelve configurations in which the linear program (8) stays above the bound (4): the number of
admissible directions, the bound, the exact minimax value over coset designs and its ratio to the bound (both
certified in rational arithmetic, results/lp_certificates.json), the number of balanced samples that contain
the origin and how many of them are flats (results/attainment.json), and how completeness is established.
Writes results/tableS9.md.
"""
import json, os

OUT = os.path.dirname(os.path.abspath(__file__))
lp = json.load(open(f"{OUT}/results/lp_certificates.json"))["nonattainment_certificates"]
att = {(a["q"], a["m"], a["k"], a["margins"]): a for a in json.load(open(f"{OUT}/results/attainment.json"))}
MARGINS = {"one-way": "one", "one- and two-way": "one and two"}


def num(n):                                   # thousands separated by a space from five digits on, as in the text
    return f"{n:,}".replace(",", " ") if n >= 10000 else str(n)


def frac(s):
    return f"${s}$" if "/" not in s else "${}/{}$".format(*s.split("/"))


rows = []
for c in lp:
    assert c["primal_verified"] and c["dual_verified"] and c["value_exact"]
    a = att[(c["q"], c["m"], c["k"], c["balanced_on"])]
    if "balanced_samples_containing_0" in a:
        n0, flats = a["balanced_samples_containing_0"], a["flats"]
        assert a["all_flats"] and flats == c["n_admissible"]
        how = "analytic, enumerated" if c["k"] <= 4 else "enumerated"
    else:                                      # OA(16,7,2,2): two enumerations agreeing pair by pair (crosscheck_oa16.py)
        x = next(x for x in json.load(open(f"{OUT}/results/crosscheck_oa16.json")) if x["m"] == c["m"])
        n0, flats, how = x["direct_samples_containing_0"], x["flats_containing_0"], "enumerated twice"
        assert flats == c["n_admissible"] and a["minimax_all_designs"] == c["minimax_value"]
    rows.append((c["balanced_on"] != "one-way", c["k"] == 16, c["q"], c["m"], c["k"],
                 f"| ${c['q']}^{c['m']}$ | {MARGINS[c['balanced_on']]} | {c['k']} | {c['n_admissible']} | "
                 f"{frac(c['bound'])} | {frac(c['minimax_value'])} | {c['ratio']:.3f} | {num(n0)} | {num(flats)} | {how} |"))
rows.sort()
lines = ["| Frame | Margins | $k$ | Directions | Bound | Minimax | Ratio | Samples with 0 | Flats | Completeness |", "|---|---|---|---|---|---|---|---|---|---|"] + [r[-1] for r in rows]
open(f"{OUT}/results/tableS9.md", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
