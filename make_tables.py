"""Tables 1-3 of the manuscript (Markdown) from results/, with Monte Carlo standard errors and exact/MC labels.

Writes results/table1.md to results/table4.md and results/grid_balance_shares.json (exact shares of
coset samples on the grid that are balanced on S2 and on S4).
"""
import csv, glob, json, math, os
import numpy as np
from affsamp import Space, admissible_directions, main_effect_set
OUT = os.path.dirname(os.path.abspath(__file__))

def rd(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))

def table1():
    keep = {("2", "7", "one-way"), ("2", "7", "one- and two-way"), ("3", "5", "one-way"), ("3", "4", "one- and two-way")}
    out = ["| Frame | Balanced on | $k$ | $\\lvert\\mathcal A_h(\\mathcal S)\\rvert$ | Bound $(Q-1)/\\lvert R_{\\mathcal S}\\rvert$ | Uniform $\\nu$ | Minimax $\\nu$ |",
           "|---|---|---|---|---|---|---|"]
    for r in rd(f"{OUT}/results/minimax_lp.csv"):
        if (r["q"], r["m"], r["balance"]) not in keep:
            continue
        frame = f"${r['q']}^{r['m']}$"
        bal = "one-variable margins" if r["balance"] == "one-way" else "one- and two-variable margins"
        if r["n_admissible"] == "0":
            out.append(f"| {frame} | {bal} | {r['k']} | 0 | | none | none |")
        else:
            out.append(f"| {frame} | {bal} | {r['k']} | {r['n_admissible']} | {float(r['flat_bound']):.4f} | "
                       f"{float(r['uniform_ratio']):.3f} | {float(r['minimax_ratio']):.3f} |")
    return "\n".join(out)

def table2():
    ex = rd(f"{OUT}/results/factorial_exact.csv")
    lat = {(r["p"], r["k"]): r for r in rd(f"{OUT}/results/factorial_lattice.csv")}
    cube = {}
    for p, m in [(2, 7), (3, 5)]:
        for r in rd(f"{OUT}/results/factorial_cube_p{p}_m{m}.csv"):
            cube[(str(p), r["k"], r["response"])] = r
    out = ["| Frame | $k$ | Response | $R^2_{\\mathrm{adj}}$ | Cube method, MC (s.e.) | Coset, uniform $\\nu$ | Coset, minimax $\\nu$ | Coset, $r=q$ | Lattice sampling (strength) | Cube: share exactly balanced, MC (s.e.) |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for r in ex:
        c = cube[(r["p"], r["k"], r["response"])]
        L = lat[(r["p"], r["k"])]
        vs = float(r["var_srs"]); B = int(c["B"]); ph = float(c["prop_exact"])
        rr = r.get("deff_coset-uniform r=2") or r["deff_coset-uniform r=3"]
        out.append(f"| ${r['p']}^{r['m']}$ | {r['k']} | {r['response']} | {float(r['adjR2']):.3f} | "
                   f"{float(c['mse'])/vs:.3f} ({float(c['mse_se'])/vs:.3f}) | {float(r['deff_coset-uniform']):.3f} | "
                   f"{float(r['deff_coset-minimax']):.3f} | {float(rr):.3f} | {float(L['deff_' + r['response']]):.3f} ({L['strength']}) | "
                   f"{ph:.3f} ({math.sqrt(ph*(1-ph)/B):.3f}) |")
    return "\n".join(out)

def grid_balance_shares():
    """Share of admissible lines (h=1) for each S_j whose flats are also balanced on S2 and on S4 (exact)."""
    q = 61; sp = Space(q, 2)
    ROWS, COLS, ANTI, DIAG = np.array([[1, 0]]), np.array([[0, 1]]), np.array([[1, 1]]), np.array([[1, q - 1]])
    S = {"S1 rows": [ROWS], "S2 rows+cols": [ROWS, COLS], "S3 +anti-diagonals": [ROWS, COLS, ANTI],
         "S4 +diagonals": [ROWS, COLS, ANTI, DIAG]}
    M2 = main_effect_set(sp, S["S2 rows+cols"]); M4 = main_effect_set(sp, S["S4 +diagonals"])
    shares = {}
    for name, fac in S.items():
        A = admissible_directions(sp, fac, 1)
        b2 = sum(not M2[sp.span(D)].any() for D in A) / len(A)
        b4 = sum(not M4[sp.span(D)].any() for D in A) / len(A)
        shares[name] = [b2, b4, len(A)]
    return shares

def table3(sh):
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    srs = {r["k"]: r["var"] for r in d["designs"] if r["design"] == "SRS"}
    rows = ["| Design | Balanced on | Evaluation | $k=61$ | $k=183$ | Residual worst-case factor | Share exactly balanced on $\\mathcal S_2$ / $\\mathcal S_4$ |",
            "|---|---|---|---|---|---|---|",
            "| Simple random sampling | none | exact | 1 | 1 | infinite | below $10^{-6}$ / below $10^{-6}$ |"]
    name_s = {"S1 rows": "$\\mathcal S_1$ (rows)", "S2 rows+cols": "$\\mathcal S_2$ (rows, columns)",
              "S3 +anti-diagonals": "$\\mathcal S_3$ ($\\mathcal S_2$ + anti-diagonals)", "S4 +diagonals": "$\\mathcal S_4$ ($\\mathcal S_3$ + diagonals)"}
    bal_s = {k: f"{v[0]:.3f} / {v[1]:.3f}" for k, v in sh.items()}
    ex = {}
    for r in d["designs"]:
        if r["design"] != "SRS":
            ex.setdefault((r["design"], r["strat"]), {})[r["k"]] = r
    for (des, strat), byk in ex.items():
        if des == "coset flat":
            name, wc = "Equal-weight coset design", "1"
        else:
            L = int(des.split("best-")[1].split()[0])
            name, wc = f"Coset design, best {L} lattices", f"{58 / L:.2f}"
        rows.append(f"| {name} | {name_s[strat]} | exact | {byk[61]['deff']:.4f} | {byk[183]['deff']:.4f} | {wc} | {bal_s[strat]} |")
    mc = {}
    for path in glob.glob(f"{OUT}/results/mc_grid/*.csv"):
        for r in rd(path):
            mc.setdefault(r["design"], {})[int(r["k"])] = r
    label = {"cubeS2": ("Cube method", name_s["S2 rows+cols"]), "cubeS4": ("Cube method", name_s["S4 +diagonals"]),
             "lpm2": ("Local pivotal method", "none (spatial spread)"), "lcubeS4": ("Local cube method", name_s["S4 +diagonals"] + " and spread")}
    for key in ["cubeS2", "cubeS4", "lpm2", "lcubeS4"]:
        cells, bal = [], []
        for k in (61, 183):
            r = mc.get(key, {}).get(k)
            if r is None:
                cells.append("[TO ADD]"); bal.append("[TO ADD]"); continue
            B = int(r["B"]); p2, p4 = float(r["prop_exact_S2"]), float(r["prop_exact_S4"])
            cells.append(f"{float(r['mse'])/srs[k]:.4f} ({float(r['mse_se'])/srs[k]:.4f})")
            share = lambda p: "none" if p == 0 else f"{p:.3f} ({math.sqrt(p*(1-p)/B):.3f})"
            bal.append(f"{share(p2)} / {share(p4)} at $k={k}$")
        B_ = int(mc.get(key, {}).get(61, {}).get("B", 0))
        B_ = f"{B_:,}".replace(",", " ") if B_ >= 10000 else str(B_)
        rows.append(f"| {label[key][0]} | {label[key][1]} | Monte Carlo, {B_} draws (s.e.) | {cells[0]} | {cells[1]} | infinite | {'; '.join(bal)} |")
    audit = f"{OUT}/results/cube_audit_S2.csv"
    if os.path.exists(audit):   # cube method with the flight phase continued to a vertex (audit_cube.py)
        au = {int(r["k"]): r for r in rd(audit) if r["landing"] == "continued"}
        share = lambda p, B: "none" if p == 0 else ("all" if p == 1 else f"{p:.3f} ({math.sqrt(p*(1-p)/B):.3f})")
        cells = [f"{float(au[k]['deff']):.4f} ({float(au[k]['deff_se']):.4f})" for k in (61, 183)]
        bal = [f"{share(float(au[k]['share_exact_S2']), int(au[k]['B']))} / {share(float(au[k]['share_exact_S4']), int(au[k]['B']))} at $k={k}$"
               for k in (61, 183)]
        rows.append(f"| Cube method, flight continued to a vertex | {name_s['S2 rows+cols']} | Monte Carlo, {int(au[61]['B'])} draws (s.e.) | "
                    f"{cells[0]} | {cells[1]} | finite, not computed | {'; '.join(bal)} |")
    return "\n".join(rows)


def table4():
    """Exact evaluation of the variance estimator of Proposition 4(b) for the r = 3 grid designs."""
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    name_s = {"S1 rows": "$\\mathcal S_1$", "S2 rows+cols": "$\\mathcal S_2$",
              "S3 +anti-diagonals": "$\\mathcal S_3$", "S4 +diagonals": "$\\mathcal S_4$"}
    rows = ["| Design ($k=183$, $r=3$) | Balanced on | Coefficient of variation of $v$ | Coverage of the nominal 95% interval | Share of samples with $v=0$ |",
            "|---|---|---|---|---|"]
    for r in d["variance_estimation"]:
        name = "Equal-weight coset design" if r["design"] == "flat" else f"Coset design, best {r['design'][4:]} lattices"
        rows.append(f"| {name} | {name_s[r['strat']]} | {r['cv_v']:.2f} | {100 * r['coverage_t95']:.1f}% | {100 * r['share_v_zero']:.3f}% |")
    return "\n".join(rows)

if __name__ == "__main__":
    sh = grid_balance_shares()
    json.dump(sh, open(f"{OUT}/results/grid_balance_shares.json", "w"), indent=1)
    open(f"{OUT}/results/table1.md", "w").write(table1() + "\n")
    open(f"{OUT}/results/table2.md", "w").write(table2() + "\n")
    open(f"{OUT}/results/table3.md", "w").write(table3(sh) + "\n")
    open(f"{OUT}/results/table4.md", "w").write(table4() + "\n")
    print(table1()); print(); print(table3(sh)); print(); print(table4())
