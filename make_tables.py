"""Tables 1-5 and A1 of the manuscript (Markdown) from results/, with Monte Carlo standard errors and exact/MC labels.

Writes results/table1.md to results/table5.md, results/tableA1.md and results/grid_balance_shares.json (exact
shares of grid coset samples that are balanced on S2 and on S4).
"""
import csv, glob, json, math, os
import numpy as np
from affsamp import Space, CosetDesign, admissible_directions, main_effect_set
OUT = os.path.dirname(os.path.abspath(__file__))
RESPONSES = ["additive", "mixed", "interaction", "higher-order", "noise"]
q61 = 61


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


def factorial_inputs():
    ex = rd(f"{OUT}/results/factorial_exact.csv")
    lat = {(r["p"], r["k"]): r for r in rd(f"{OUT}/results/factorial_lattice.csv")}
    cube = {}
    for p, m in [(2, 7), (3, 5)]:
        for r in rd(f"{OUT}/results/factorial_cube_p{p}_m{m}.csv"):
            cube[(str(p), r["k"], r["response"].replace(".", "-"))] = r     # R writes higher-order as higher.order
    cont = {(r["q"], r["k"], r["response"]): r for r in rd(f"{OUT}/results/cube_factorial.csv") if r["landing"] == "continued"}
    return ex, lat, cube, cont


def table2():
    ex, lat, cube, cont = factorial_inputs()
    out = ["| Frame | $k$ | Response | $R^2_{\\mathrm{adj}}$ | Cube method, MC (s.e.) | Cube, flight continued, MC (s.e.) | Coset, uniform $\\nu$ | Coset, minimax $\\nu$ | Coset, $r=q$ | Lattice sampling (strength) |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    order = {n: i for i, n in enumerate(RESPONSES)}
    for r in sorted(ex, key=lambda r: (r["p"], int(r["k"]), order[r["response"]])):
        c = cube[(r["p"], r["k"], r["response"])]
        cc = cont[(r["p"], r["k"], r["response"])]
        L = lat[(r["p"], r["k"])]
        vs = float(r["var_srs"])
        rr = r.get("deff_coset-uniform r=2") or r["deff_coset-uniform r=3"]
        out.append(f"| ${r['p']}^{r['m']}$ | {r['k']} | {r['response']} | {float(r['adjR2']):.3f} | "
                   f"{float(c['mse'])/vs:.3f} ({float(c['mse_se'])/vs:.3f}) | {float(cc['deff']):.3f} ({float(cc['deff_se']):.3f}) | "
                   f"{float(r['deff_coset-uniform']):.3f} | {float(r['deff_coset-minimax']):.3f} | {float(rr):.3f} | "
                   f"{float(L['deff_' + r['response']]):.3f} ({L['strength']}) |")
    return "\n".join(out)


def table3():
    """Properties of the factorial designs that do not depend on the response."""
    _, lat, cube, cont = factorial_inputs()
    info = json.load(open(f"{OUT}/results/factorial_designs.json"))
    wc = {(str(i["p"]), str(i["k"]), i["design"]): i["max_w_over_flat"] for i in info}
    out = ["| Frame | $k$ | Cube: exactly balanced, package / flight continued, MC (s.e.) | Cube: main-effect weight share $\\tau$, package / flight continued, MC | Worst-case factor: coset, uniform $\\nu$ | Coset, $r=q$ | Lattice sampling |",
           "|---|---|---|---|---|---|---|"]
    for p, m, ks in [("2", "7", ("8", "16", "32")), ("3", "5", ("9", "27", "81"))]:
        for k in ks:
            c, cc = cube[(p, k, "additive")], cont[(p, k, "additive")]
            B = int(c["B"]); ph = float(c["prop_exact"])
            rr = wc.get((p, k, "coset-uniform r=2"), wc.get((p, k, "coset-uniform r=3")))
            out.append(f"| ${p}^{m}$ | {k} | {ph:.3f} ({math.sqrt(ph * (1 - ph) / B):.3f}) / "
                       f"{float(cc['share_exact']):.3f} ({float(cc['share_exact_se']):.3f}) | "
                       f"{float(c['tau']):.4f} / {float(cc['tau']):.4f} | {wc[(p, k, 'coset-uniform')]:.3f} | {rr:.3f} | "
                       f"{float(lat[(p, k)]['worst_case_over_flat']):.2f} |")
    return "\n".join(out)


def grid_balance_shares():
    """Share of admissible lines (h=1) for each S_j whose flats are also balanced on S2 and on S4 (exact)."""
    sp = Space(q61, 2)
    S = grid_strata()
    M2 = main_effect_set(sp, S["S2 rows+cols"]); M4 = main_effect_set(sp, S["S4 +diagonals"])
    shares = {}
    for name, fac in S.items():
        A = admissible_directions(sp, fac, 1)
        b2 = sum(not M2[sp.span(D)].any() for D in A) / len(A)
        b4 = sum(not M4[sp.span(D)].any() for D in A) / len(A)
        shares[name] = [b2, b4, len(A)]
    return shares


def grid_strata():
    ROWS, COLS, ANTI, DIAG = np.array([[1, 0]]), np.array([[0, 1]]), np.array([[1, 1]]), np.array([[1, q61 - 1]])
    return {"S1 rows": [ROWS], "S2 rows+cols": [ROWS, COLS], "S3 +anti-diagonals": [ROWS, COLS, ANTI],
            "S4 +diagonals": [ROWS, COLS, ANTI, DIAG]}


def grid_tau_exact():
    """Main-effect weight share tau on S4 of the grid coset designs (independent of r): sum_{M_S4} c_nu / (Q - 1)."""
    sp = Space(q61, 2)
    S = grid_strata()
    M4 = main_effect_set(sp, S["S4 +diagonals"])
    tau = {}
    for name, fac in S.items():
        d = CosetDesign(sp, admissible_directions(sp, fac, 1), r=1)
        tau[("coset flat", name)] = float(d.coverage()[M4].sum() / (d.Q - 1))
    return tau


def tau_from_samples(samples, k):
    """Mean main-effect weight share on S4 of samples of the grid (unit code 61 i + j)."""
    T = q61 * q61 / k - 1
    vals = []
    for s in samples:
        i, j = np.divmod(s, q61)
        tm = sum(((np.bincount(v % q61, minlength=q61) - k / q61) ** 2).sum() for v in (i, j, i + j, i - j)) * q61 / k ** 2
        vals.append(tm / T)
    return float(np.mean(vals)), float(np.std(vals, ddof=1) / math.sqrt(len(vals)))


def table4(sh):
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    srs = {r["k"]: r["var"] for r in d["designs"] if r["design"] == "SRS"}
    tau_ex = grid_tau_exact()
    tau_srs = 240 / (q61 * q61 - 1)
    rows = ["| Design | Balanced on | Evaluation | $k=61$ | $k=183$ | Residual worst-case factor | Main-effect weight share $\\tau$ on $\\mathcal S_4$, $k=61$ / $k=183$ | Exactly balanced on $\\mathcal S_2$ / $\\mathcal S_4$, $k=61$ | Exactly balanced on $\\mathcal S_2$ / $\\mathcal S_4$, $k=183$ |",
            "|---|---|---|---|---|---|---|---|---|",
            f"| Simple random sampling | none | exact | 1 | 1 | infinite | {tau_srs:.4f} / {tau_srs:.4f} | below $10^{{-6}}$ / below $10^{{-6}}$ | below $10^{{-6}}$ / below $10^{{-6}}$ |"]
    name_s = {"S1 rows": "$\\mathcal S_1$ (rows)", "S2 rows+cols": "$\\mathcal S_2$ (rows, columns)",
              "S3 +anti-diagonals": "$\\mathcal S_3$ ($\\mathcal S_2$ + anti-diagonals)", "S4 +diagonals": "$\\mathcal S_4$ ($\\mathcal S_3$ + diagonals)"}
    bal_s = {k: f"{v[0]:.3f} / {v[1]:.3f}" for k, v in sh.items()}
    ex = {}
    for r in d["designs"]:
        if r["design"] != "SRS":
            ex.setdefault((r["design"], r["strat"]), {})[r["k"]] = r
    for (des, strat), byk in ex.items():
        if des == "coset flat":
            name, wc, t = "Equal-weight coset design", "1", tau_ex[("coset flat", strat)]
        else:
            L = int(des.split("best-")[1].split()[0])
            name, wc, t = f"Coset design, best {L} lattices", f"{58 / L:.2f}", 0.0
        ts = "0" if t == 0 else f"{t:.4f}"
        rows.append(f"| {name} | {name_s[strat]} | exact | {byk[61]['deff']:.4f} | {byk[183]['deff']:.4f} | {wc} | {ts} / {ts} | {bal_s[strat]} | {bal_s[strat]} |")
    mc = {}
    for path in glob.glob(f"{OUT}/results/mc_grid/*.csv") + glob.glob(f"{OUT}/results/mc_grid_v2/lpm2torus_*.csv"):
        for r in rd(path):
            mc.setdefault(r["design"], {})[int(r["k"])] = r
    v2 = {}
    for path in glob.glob(f"{OUT}/results/mc_grid_v2/*.csv"):
        for r in rd(path):
            v2.setdefault(r["design"], {})[int(r["k"])] = r
    label = {"cubeS2": ("Cube method", name_s["S2 rows+cols"]), "cubeS4": ("Cube method", name_s["S4 +diagonals"]),
             "lpm2": ("Local pivotal method", "none (spatial spread, planar distances)"),
             "lpm2torus": ("Local pivotal method", "none (spatial spread, wrapped distances)"),
             "lcubeS4": ("Local cube method", name_s["S4 +diagonals"] + " and spread")}
    share = lambda p, B: "none" if p == 0 else ("all" if p == 1 else f"{p:.3f} ({math.sqrt(p*(1-p)/B):.3f})")
    for key in ["cubeS2", "cubeS4", "lpm2", "lpm2torus", "lcubeS4"]:
        cells, bal, taus = [], [], []
        for k in (61, 183):
            r = mc.get(key, {}).get(k)
            if r is None:
                cells.append("[TO ADD]"); bal.append("[TO ADD]"); taus.append("[TO ADD]"); continue
            B = int(r["B"]); p2, p4 = float(r["prop_exact_S2"]), float(r["prop_exact_S4"])
            cells.append(f"{float(r['mse'])/srs[k]:.4f} ({float(r['mse_se'])/srs[k]:.4f})")
            bal.append(f"{share(p2, B)} / {share(p4, B)}")
            t = v2.get(key, {}).get(k)
            taus.append(f"{float(t['tau_S4']):.4f}" if t is not None else "[TO ADD]")
        B_ = int(mc.get(key, {}).get(61, {}).get("B", 0))
        B_ = f"{B_:,}".replace(",", " ") if B_ >= 10000 else str(B_)
        rows.append(f"| {label[key][0]} | {label[key][1]} | Monte Carlo, {B_} draws (s.e.) | {cells[0]} | {cells[1]} | infinite | {taus[0]} / {taus[1]} | {bal[0]} | {bal[1]} |")
    audit = f"{OUT}/results/cube_audit_S2.csv"
    au = {int(r["k"]): r for r in rd(audit) if r["landing"] == "continued"}
    draws = [r for r in rd(f"{OUT}/results/cube_audit_S2_draws.csv") if r["landing"] == "continued"]
    tc = {k: tau_from_samples([np.array(r["sample"].split(), dtype=int) for r in draws if int(r["k"]) == k], k)[0] for k in (61, 183)}
    cells = [f"{float(au[k]['deff']):.4f} ({float(au[k]['deff_se']):.4f})" for k in (61, 183)]
    bal = [f"{share(float(au[k]['share_exact_S2']), int(au[k]['B']))} / {share(float(au[k]['share_exact_S4']), int(au[k]['B']))}" for k in (61, 183)]
    rows.append(f"| Cube method, flight continued to a vertex | {name_s['S2 rows+cols']} | Monte Carlo, {int(au[61]['B'])} draws (s.e.) | "
                f"{cells[0]} | {cells[1]} | finite, not computed | {tc[61]:.4f} / {tc[183]:.4f} | {bal[0]} | {bal[1]} |")
    return "\n".join(rows)


def table5():
    """Exact evaluation of the variance estimators (varest_exact.py and the grid results of study_grid.py)."""
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    name_s = {"S1 rows": "$\\mathcal S_1$", "S2 rows+cols": "$\\mathcal S_2$", "S3 +anti-diagonals": "$\\mathcal S_3$", "S4 +diagonals": "$\\mathcal S_4$"}
    rows = ["| Population and design | $k$ | Estimator | Degrees of freedom | Coefficient of variation | Coverage of the nominal 95% interval | Share of samples with estimate 0 |",
            "|---|---|---|---|---|---|---|"]
    for r in d["variance_estimation"]:
        name = "equal-weight coset design" if r["design"] == "flat" else f"best {r['design'][4:]} lattices"
        rows.append(f"| Grid, {name}, {name_s[r['strat']]} | 183 | random groups | 2 | {r['cv_v']:.2f} | {100 * r['coverage_t95']:.1f}% | {100 * r['share_v_zero']:.3f}% |")
    ve = rd(f"{OUT}/results/varest_exact.csv")
    g = [r for r in ve if r["population"] == "61x61 grid" and r["estimator"] == "regression"][0]
    rows.append(f"| Grid, equal-weight coset design, $\\mathcal S_1$ | 183 | regression | {g['df_min']} | {float(g['cv']):.2f} | "
                f"{100 * float(g['coverage_t95']):.1f}% | 0 |")
    for pop, frame in (("2^7", "$2^7$"), ("3^5", "$3^5$")):
        for k in sorted({int(r["k"]) for r in ve if r["population"] == pop}):
            for est in ("random groups", "regression"):
                sub = [r for r in ve if r["population"] == pop and int(r["k"]) == k and r["estimator"] == est]
                if not sub:
                    continue
                cv = [float(r["cv"]) for r in sub]; cov = [100 * float(r["coverage_t95"]) for r in sub]
                df = (sub[0]["df_min"] if sub[0]["df_min"] == sub[0]["df_max"] else f"{sub[0]['df_min']} to {sub[0]['df_max']}")
                z = max(float(r["share_zero"]) for r in sub)
                rows.append(f"| {frame}, $C(h-1,\\nu,q)$, five responses | {k} | {est} | {df} | {min(cv):.2f} to {max(cv):.2f} | "
                            f"{min(cov):.1f}% to {max(cov):.1f}% | {'0' if z == 0 else f'{100 * z:.3f}%'} |")
    return "\n".join(rows)


def tableA1():
    """Replicated populations: median and 5%-95% range of ratios of design effects to the minimax coset design."""
    rows = rd(f"{OUT}/results/replicated_populations.csv")
    get = {(r["q"], r["k"], r["response"], r["design"], r["statistic"]): r for r in rows}
    P = rows[0]["P"]
    designs = [("coset-uniform", "Coset, uniform $\\nu$"), ("coset-replicated", "Coset, $r=q$"), ("lattice", "Lattice sampling"),
               ("cube", "Cube method"), ("cube-continued", "Cube, flight continued")]
    out = [f"| Frame | $k$ | Response | " + " | ".join(n for _, n in designs) + " |", "|---|---|---|" + "---|" * len(designs)]
    for q, m, ks in [("2", "7", ("8", "16", "32")), ("3", "5", ("9", "27", "81"))]:
        for k in ks:
            for resp in RESPONSES:
                cells = []
                for key, _ in designs:
                    r = get.get((q, k, resp, key, "ratio_to_minimax"))
                    cells.append("" if r is None else f"{float(r['q50']):.3f} ({float(r['q05']):.3f}, {float(r['q95']):.3f})")
                out.append(f"| ${q}^{m}$ | {k} | {resp} | " + " | ".join(cells) + " |")
    return "\n".join(out), P


if __name__ == "__main__":
    sh = grid_balance_shares()
    json.dump(sh, open(f"{OUT}/results/grid_balance_shares.json", "w"), indent=1)
    t = {"table1": table1(), "table2": table2(), "table3": table3(), "table4": table4(sh), "table5": table5()}
    if os.path.exists(f"{OUT}/results/replicated_populations.csv"):
        t["tableA1"] = tableA1()[0]
    for name, body in t.items():
        open(f"{OUT}/results/{name}.md", "w").write(body + "\n")
        print(body, "\n")
