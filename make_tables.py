"""Tables 2-5, A.1 and S1-S7 of the manuscript (Markdown) from results/, with Monte Carlo standard errors and
exact/MC labels, in the numbering of the shortened version of the manuscript (release v1.9).

Writes results/table2.md to results/table5.md, results/tableA.1.md, results/tableS1.md to results/tableS7.md and
results/grid_balance_shares.json (exact shares of grid coset samples that are balanced on S2 and on S4).  Table 1
of the manuscript (prior work) is not computed.  The functions keep the names of an earlier numbering:
  table1 -> Table 2 (minimax mixtures);
  table3_main -> Table 3 (factorial design effects: the cube method of both packages, the uniform and minimax coset
      designs and lattice sampling), table2 -> Table S3 (the same with the continued flight, both landings of
      samplecube and the replicated coset design), table3 -> Table S4 (properties that do not depend on the response);
  table4 -> the full grid table, which grid_tables splits into Table 4 (design effects), Table 5 (worst case over
      Y_1(sigma^2) for S4 and main-effect weight shares) and Table S5 (all worst cases and shares of exactly balanced
      samples); tableS2 -> Table S2 (implementation diagnostics);
  table_order -> Table S6 (order of the balancing variables; written only once the runs in reverse order exist),
  table5 -> Table S7 (variance estimators; Table S6 when there is no Table S6 for the order),
  tableA1 -> Table S1 (replicated populations), table_certificate -> Table A.1.
Release v1.8 wrote the earlier numbering (Tables 2-7, A.1, S1 and S2).
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
    ht = {(r["q"], r["k"], r["source"], r["response"]): r for r in rd(f"{OUT}/results/factorial_ht.csv")}
    return ex, lat, cube, cont, ht


def factorial_table(full):
    """Design effects of the Horvitz-Thompson mean; the cube variants from their stored samples (ht_factorial.py).
    full=False: Table 3, with the cube method of BalancedSampling and samplecube with the landing by linear
    programming, the uniform and minimax coset designs and lattice sampling.  full=True: Table S3, which adds the
    continued flight, samplecube with the landing by suppression of variables and the replicated design."""
    ex, lat, _, _, ht = factorial_inputs()
    has = lambda src: any(key[2] == src for key in ht)
    mc = [("Cube", "package")]
    if full:
        mc.append(("Cube, continued", "transcription, flight continued"))
    if has("samplecube, LP landing"):
        mc.append(("Cube, sampling" if not full else "Cube, sampling, LP landing", "samplecube, LP landing"))
    if full and has("samplecube, suppression"):
        mc.append(("Cube, sampling, suppression", "samplecube, suppression"))
    ex_cols = ["Uniform", "Minimax"] + (["$r=q$"] if full else []) + ["Lattice"]
    head = ["Frame", "$k$", "Response", "$R^2_{\\mathrm{adj}}$"] + [n for n, _ in mc] + ex_cols
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    order = {n: i for i, n in enumerate(RESPONSES)}
    for r in sorted(ex, key=lambda r: (r["p"], int(r["k"]), order[r["response"]])):
        cells = [f"${r['p']}^{r['m']}$", r["k"], r["response"], f"{float(r['adjR2']):.3f}"]
        for _, src in mc:
            c = ht[(r["p"], r["k"], src, r["response"])]
            cells.append(f"{float(c['deff_ht']):.3f} ({float(c['deff_ht_se']):.3f})")
        cells += [f"{float(r['deff_coset-uniform']):.3f}", f"{float(r['deff_coset-minimax']):.3f}"]
        if full:
            cells.append(f"{float(r.get('deff_coset-uniform r=2') or r['deff_coset-uniform r=3']):.3f}")
        L = lat[(r["p"], r["k"])]
        cells.append(f"{float(L['deff_' + r['response']]):.3f} ({L['strength']})")
        out.append("| " + " | ".join(str(c) for c in cells) + " |")
    return "\n".join(out)


def table3_main():
    return factorial_table(False)


def table2():
    return factorial_table(True)


def table_order():
    """Table S7: design effects for the additive response, and shares of exactly balanced samples, with the
    balancing variables in the given and in the reverse order (the inclusion probabilities first in both), for the
    cube method of BalancedSampling and samplecube with both landings (ht_factorial.py).  None if the runs in
    reverse order are missing."""
    ht = factorial_inputs()[4]
    impl = [("BalancedSampling", "package", "package, reversed order"),
            ("samplecube, LP landing", "samplecube, LP landing", "samplecube, LP landing, reversed order"),
            ("samplecube, suppression", "samplecube, suppression", "samplecube, suppression, reversed order")]
    impl = [x for x in impl if any(key[2] == x[2] for key in ht)]
    if not impl:
        return None
    head = ["Frame", "$k$"] + [f"{n}, design effect" for n, _, _ in impl] + [f"{n}, exactly balanced" for n, _, _ in impl]
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for p, m, ks in [("2", "7", ("8", "16", "32")), ("3", "5", ("9", "27", "81"))]:
        for k in ks:
            cells = [f"${p}^{m}$", k]
            for _, a, b in impl:
                ca, cb = ht[(p, k, a, "additive")], ht[(p, k, b, "additive")]
                cells.append(f"{float(ca['deff_ht']):.3f} ({float(ca['deff_ht_se']):.3f}) / {float(cb['deff_ht']):.3f} ({float(cb['deff_ht_se']):.3f})")
            for _, a, b in impl:
                ca, cb = ht[(p, k, a, "additive")], ht[(p, k, b, "additive")]
                cells.append(f"{float(ca['share_exact']):.3f} / {float(cb['share_exact']):.3f}")
            out.append("| " + " | ".join(str(c) for c in cells) + " |")
    return "\n".join(out)


def table3():
    """Properties of the factorial designs that do not depend on the response."""
    _, lat, cube, cont, ht = factorial_inputs()
    info = json.load(open(f"{OUT}/results/factorial_designs.json"))
    wc = {(str(i["p"]), str(i["k"]), i["design"]): i["max_w_over_flat"] for i in info}
    out = ["| Frame | $k$ | Exactly balanced | Size $k\\pm1$ | $\\tau$ | $W$, uniform | $W$, $r=q$ | $W$, lattice |",
           "|---|---|---|---|---|---|---|---|"]
    for p, m, ks in [("2", "7", ("8", "16", "32")), ("3", "5", ("9", "27", "81"))]:
        for k in ks:
            src = [ht[(p, k, s_, "additive")] for s_ in ("package", "transcription, flight continued", "samplecube, LP landing",
                                                         "samplecube, suppression")
                   if (p, k, s_, "additive") in ht]
            rr = wc.get((p, k, "coset-uniform r=2"), wc.get((p, k, "coset-uniform r=3")))
            out.append(f"| ${p}^{m}$ | {k} | " + " / ".join(f"{float(c['share_exact']):.3f} ({float(c['share_exact_se']):.3f})" for c in src) + " | "
                       + " / ".join(str(c["size_exceptions"]) for c in src) + " | "
                       + " / ".join(f"{float(c['tau']):.4f}" for c in src) + f" | {wc[(p, k, 'coset-uniform')]:.3f} | {rr:.3f} | "
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
    """Mean main-effect weight share on S4 of samples of the grid (unit code 61 i + j), with its standard error:
    t_M = sum over the four line stratifications of (q sum_b n_b^2 - |s|^2)/k^2, that is sum_{u in M_S4} |eta_s(u)|^2
    with eta_s(u) = k^-1 sum_{x in s} chi_u(x)."""
    T = q61 * q61 / k - 1
    vals = []
    for s in samples:
        i, j = np.divmod(s, q61)
        tm = sum(q61 * (np.bincount(v % q61, minlength=q61) ** 2).sum() - len(s) ** 2 for v in (i, j, i + j, i - j)) / k ** 2
        vals.append(tm / T)
    return float(np.mean(vals)), float(np.std(vals, ddof=1) / math.sqrt(len(vals)))


def grid_prop3_exact(C=1.0):
    """Worst case over Y_C(sigma^2) for S4 of the equal-weight coset designs of the grid, relative to the
    equal-weight design for S4 (Proposition 3(b)): (C max_{M_S4} w + max_{R_S4} w) / ((Q-1)/|R_S4|); it does not
    depend on r.  For flats w = c_nu."""
    sp = Space(q61, 2)
    S = grid_strata()
    M4 = main_effect_set(sp, S["S4 +diagonals"])
    R4 = ~M4; R4[0] = False
    out = {}
    for name, fac in S.items():
        d = CosetDesign(sp, admissible_directions(sp, fac, 1), r=1)
        c = d.coverage()
        out[name] = float((C * c[M4].max() + c[R4].max()) / ((d.Q - 1) / R4.sum()))
    return out, int(R4.sum()), int(M4[1:].sum())


def table4(sh):
    """Table 4.  Proposition 2 column: worst case over Y(sigma^2) for the design's own balancing set, relative to the
    equal-weight design for that set.  Proposition 3 column: worst case over Y_1(sigma^2) for S4 relative to the
    equal-weight design for S4, exact for translation-invariant designs and the plug-in lower bound
    1 + tau (|R|/d - 1) of Proposition 3(a), with its standard error, for the Monte Carlo designs."""
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    srs = {r["k"]: r["var"] for r in d["designs"] if r["design"] == "SRS"}
    tau_ex = grid_tau_exact()
    f3, nR, dS = grid_prop3_exact()
    tau_srs = dS / (q61 * q61 - 1)
    f3_srs = 2 * nR / (q61 * q61 - 1)                    # (1 + C)|R|/(N - 1) with C = 1
    lb = lambda t: 1 + t * (nR / dS - 1)
    rows = ["| Design | Balanced on | Evaluation | $k=61$ | $k=183$ | Worst case over $\\mathcal Y(\\sigma^2)$ for its own $\\mathcal S$ (Proposition 2) | Worst case over $\\mathcal Y_1(\\sigma^2)$ for $\\mathcal S_4$ (Proposition 3), $k=61$ / $k=183$ | Main-effect weight share $\\tau$ on $\\mathcal S_4$, $k=61$ / $k=183$ | Exactly balanced on $\\mathcal S_2$ / $\\mathcal S_4$, $k=61$ | Exactly balanced on $\\mathcal S_2$ / $\\mathcal S_4$, $k=183$ |",
            "|---|---|---|---|---|---|---|---|---|---|",
            f"| Simple random sampling | none | exact | 1 | 1 | 1 | {f3_srs:.3f} / {f3_srs:.3f} | {tau_srs:.4f} / {tau_srs:.4f} | below $10^{{-6}}$ / below $10^{{-6}}$ | below $10^{{-6}}$ / below $10^{{-6}}$ |"]
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
            w3 = "1" if abs(f3[strat] - 1) < 1e-12 else f"{f3[strat]:.3f}"
        else:
            L = int(des.split("best-")[1].split()[0])
            name, wc, t, w3 = f"Coset design, best {L} lattices", f"{58 / L:.2f}", 0.0, f"{58 / L:.2f}"
        ts = "0" if t == 0 else f"{t:.4f}"
        rows.append(f"| {name} | {name_s[strat]} | exact | {byk[61]['deff']:.4f} | {byk[183]['deff']:.4f} | {wc} | {w3} / {w3} | {ts} / {ts} | {bal_s[strat]} | {bal_s[strat]} |")
    mc = {}
    for path in glob.glob(f"{OUT}/results/mc_grid_v3/*.csv"):
        for r in rd(path):
            mc.setdefault(r["design"], {})[int(r["k"])] = r
    label = {"cubeS2": ("Cube method", name_s["S2 rows+cols"], "infinite"),
             "cubeS4": ("Cube method", name_s["S4 +diagonals"], "infinite"),
             "lpm2": ("Local pivotal method", "none (spatial spread, planar distances)", "not computed"),
             "lpm2torus": ("Local pivotal method", "none (spatial spread, wrapped distances)", "not computed"),
             "lcubeS4": ("Local cube method", name_s["S4 +diagonals"] + " and spread", "infinite")}
    share = lambda p, B: "none" if p == 0 else ("all" if p == 1 else f"{p:.3f} ({math.sqrt(p*(1-p)/B):.3f})")
    for key in ["cubeS2", "cubeS4", "lpm2", "lpm2torus", "lcubeS4"]:
        cells, bal, taus, w3 = [], [], [], []
        for k in (61, 183):
            r = mc.get(key, {}).get(k)
            if r is None:
                cells.append("[TO ADD]"); bal.append("[TO ADD]"); taus.append("[TO ADD]"); w3.append("[TO ADD]"); continue
            B = int(r["B"]); p2, p4 = float(r["prop_exact_S2"]), float(r["prop_exact_S4"])
            cells.append(f"{float(r['mse_ht'])/srs[k]:.4f} ({float(r['mse_ht_se'])/srs[k]:.4f})")
            bal.append(f"{share(p2, B)} / {share(p4, B)}")
            t, tse = float(r["tau_S4"]), float(r["tau_S4_se"])
            taus.append(f"{t:.4f} ({tse:.4f})" if tse >= 5e-5 else f"{t:.4f} ({tse:.5f})")
            w3.append(f"$\\ge{lb(t):.3f}$ ({tse * (nR / dS - 1):.4f})")
        B_ = int(mc.get(key, {}).get(61, {}).get("B", 0))
        B_ = f"{B_:,}".replace(",", " ") if B_ >= 10000 else str(B_)
        ev = f"Monte Carlo, {B_} draws (s.e.)"
        if key.startswith(("cube", "lcube")):           # the cube and local cube methods can miss the size by one unit
            exc = [mc.get(key, {}).get(k, {}).get("size_exceptions", "[TO ADD]") for k in (61, 183)]
            ev += f"; size other than $k$ in {exc[0]} / {exc[1]}"
        rows.append(f"| {label[key][0]} | {label[key][1]} | {ev} | {cells[0]} | {cells[1]} | {label[key][2]} | {w3[0]} / {w3[1]} | {taus[0]} / {taus[1]} | {bal[0]} | {bal[1]} |")
        if key == "cubeS4":                              # the same algorithm with QR directions (cube_grid_fixed.py)
            fx = {int(r["k"]): r for r in rd(f"{OUT}/results/cube_grid_fixed.csv")} if os.path.exists(f"{OUT}/results/cube_grid_fixed.csv") else {}
            if fx:
                c_ = [f"{float(fx[k]['deff_ht']):.4f} ({float(fx[k]['deff_ht_se']):.4f})" for k in (61, 183)]
                b_ = [f"{share(float(fx[k]['share_exact_S2']), int(fx[k]['B']))} / {share(float(fx[k]['share_exact_S4']), int(fx[k]['B']))}" for k in (61, 183)]
                tse_ = [float(fx[k]["tau_S4_se"]) for k in (61, 183)]
                t_ = [f"{float(fx[k]['tau_S4']):.4f} ({e:.4f})" if e >= 5e-5 else f"{float(fx[k]['tau_S4']):.4f} ({e:.5f})" for k, e in zip((61, 183), tse_)]
                w_ = [f"$\\ge{lb(float(fx[k]['tau_S4'])):.3f}$ ({e * (nR / dS - 1):.4f})" for k, e in zip((61, 183), tse_)]
                rows.append(f"| Cube method, directions by QR factorisation | {name_s['S4 +diagonals']} | Monte Carlo, {int(fx[61]['B'])} draws (s.e.); "
                            f"size other than $k$ in {fx[61]['size_exceptions']} / {fx[183]['size_exceptions']} | {c_[0]} | {c_[1]} | infinite | "
                            f"{w_[0]} / {w_[1]} | {t_[0]} / {t_[1]} | {b_[0]} | {b_[1]} |")
            else:
                rows.append(f"| Cube method, directions by QR factorisation | {name_s['S4 +diagonals']} | [TO ADD] | [TO ADD] | [TO ADD] | infinite | [TO ADD] | [TO ADD] | [TO ADD] | [TO ADD] |")
    audit = f"{OUT}/results/cube_audit_S2.csv"
    au = {int(r["k"]): r for r in rd(audit) if r["landing"] == "continued"}
    draws = [r for r in rd(f"{OUT}/results/cube_audit_S2_draws.csv") if r["landing"] == "continued"]
    tc = {k: tau_from_samples([np.array(r["sample"].split(), dtype=int) for r in draws if int(r["k"]) == k], k) for k in (61, 183)}
    cells = [f"{float(au[k]['deff']):.4f} ({float(au[k]['deff_se']):.4f})" for k in (61, 183)]
    bal = [f"{share(float(au[k]['share_exact_S2']), int(au[k]['B']))} / {share(float(au[k]['share_exact_S4']), int(au[k]['B']))}" for k in (61, 183)]
    taus = [f"{tc[k][0]:.4f} ({tc[k][1]:.4f})" if tc[k][1] >= 5e-5 else f"{tc[k][0]:.4f} ({tc[k][1]:.5f})" for k in (61, 183)]
    w3 = [f"$\\ge{lb(tc[k][0]):.3f}$ ({tc[k][1] * (nR / dS - 1):.4f})" for k in (61, 183)]
    rows.append(f"| Cube method, flight continued to a vertex | {name_s['S2 rows+cols']} | Monte Carlo, {int(au[61]['B'])} draws (s.e.) | "
                f"{cells[0]} | {cells[1]} | finite, not computed | {w3[0]} / {w3[1]} | {taus[0]} / {taus[1]} | {bal[0]} | {bal[1]} |")
    return "\n".join(rows)


def grid_tables(sh):
    """Tables 5 and 6: the rows of the full grid table (table4) split into the design effects (Table 5) and the worst
    cases, main-effect weight shares and shares of exactly balanced samples (Table 6), with short row labels.  The
    cube and local cube methods of the package balanced on S4, whose row reduction loses the sample size, are not
    in these tables; tableS2 reports them as implementation diagnostics."""
    lines = table4(sh).split("\n")[2:]
    short = {"$\\mathcal S_1$ (rows)": "$\\mathcal S_1$", "$\\mathcal S_2$ (rows, columns)": "$\\mathcal S_2$",
             "$\\mathcal S_3$ ($\\mathcal S_2$ + anti-diagonals)": "$\\mathcal S_3$",
             "$\\mathcal S_4$ ($\\mathcal S_3$ + diagonals)": "$\\mathcal S_4$",
             "$\\mathcal S_4$ ($\\mathcal S_3$ + diagonals) and spread": "$\\mathcal S_4$"}
    t5 = ["| Design | Evaluation | $k=61$ | $k=183$ |", "|---|---|---|---|"]
    t6 = ["| Design | $W_2$ | $W_3$, $k=61$ / $k=183$ | $\\tau$, $k=61$ / $k=183$ | Exactly balanced on $\\mathcal S_2$ / $\\mathcal S_4$, $k=61$ | Same, $k=183$ |",
          "|---|---|---|---|---|---|"]
    rows5, rows6 = {}, {}
    for line in lines:
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        name, bal_on, ev, d61, d183, w2, w3, tau, b61, b183 = c
        if name == "Simple random sampling":
            label = name
        elif name.startswith("Local pivotal method"):
            label = name + (", planar distances" if "planar" in bal_on else ", wrapped distances")
        else:
            name = {"Cube method, directions by QR factorisation": "Cube method with QR directions",
                    "Cube method, flight continued to a vertex": "Cube method, flight continued"}.get(name, name)
            label = f"{name}, {short[bal_on]}"
        if label in ("Cube method, $\\mathcal S_4$", "Local cube method, $\\mathcal S_4$"):
            continue
        ev = ev.replace("Monte Carlo, ", "MC, ").replace(" draws (s.e.)", " draws")
        assert "size other than $k$ in" not in ev or ev.endswith("in 0 / 0"), (label, ev)
        ev = ev.split("; size other than")[0]
        rows5[label] = f"| {label} | {ev} | {d61} | {d183} |"
        rows6[label] = f"| {label} | {w2} | {w3} | {tau} | {b61} | {b183} |"
    sg = f"{OUT}/results/samplecube_grid.csv"
    if os.path.exists(sg):                                  # samplecube of the package sampling (samplecube_grid.py)
        scg = {(r["design"], int(r["k"])): r for r in rd(sg)}
        bc = {(r["design"], int(r["k"])): r for r in rd(f"{OUT}/results/bound_centred.csv") if r["frame"] == "61x61"}
        share = lambda p_, se_: "none" if float(p_) == 0 else ("all" if float(p_) == 1 else f"{float(p_):.3f} ({float(se_):.3f})")
        for des, bal in (("scubeS2", "$\\mathcal S_2$"), ("scubeS4", "$\\mathcal S_4$")):
            if (des, 61) not in scg or (des, 183) not in scg:
                continue
            label = f"Cube method (sampling), {bal}"
            r61, r183 = scg[(des, 61)], scg[(des, 183)]
            assert int(r61["size_exceptions"]) == 0 and int(r183["size_exceptions"]) == 0
            rows5[label] = (f"| {label} | MC, {r61['B']} draws | {float(r61['deff_ht']):.4f} ({float(r61['deff_ht_se']):.4f}) | "
                            f"{float(r183['deff_ht']):.4f} ({float(r183['deff_ht_se']):.4f}) |")
            own = "share_exact_S2" if des == "scubeS2" else "share_exact_S4"
            w2 = "finite, not computed" if float(r61[own]) == 1 and float(r183[own]) == 1 else "infinite"
            name = "samplecube S2" if des == "scubeS2" else "samplecube S4"
            b61, b183 = bc[(name, 61)], bc[(name, 183)]
            w3 = f"$\\ge{float(b61['bound']):.3f}$ ({float(b61['bound_se']):.4f}) / $\\ge{float(b183['bound']):.3f}$ ({float(b183['bound_se']):.4f})"
            fmt_t = lambda b_: f"{float(b_['tau']):.4f} ({float(b_['tau_se']):.4f})" if float(b_["tau_se"]) >= 5e-5 else f"{float(b_['tau']):.4f} ({float(b_['tau_se']):.5f})"
            rows6[label] = (f"| {label} | {w2} | {w3} | {fmt_t(b61)} / {fmt_t(b183)} | "
                            f"{share(r61['share_exact_S2'], r61['share_exact_S2_se'])} / {share(r61['share_exact_S4'], r61['share_exact_S4_se'])} | "
                            f"{share(r183['share_exact_S2'], r183['share_exact_S2_se'])} / {share(r183['share_exact_S4'], r183['share_exact_S4_se'])} |")
    order = ["Simple random sampling"] + [f"Equal-weight coset design, $\\mathcal S_{j}$" for j in (1, 2, 3, 4)] + \
            [f"Coset design, best {L} lattices, $\\mathcal S_4$" for L in (4, 8, 16)] + \
            ["Cube method, $\\mathcal S_2$", "Cube method, flight continued, $\\mathcal S_2$", "Cube method (sampling), $\\mathcal S_2$",
             "Cube method with QR directions, $\\mathcal S_4$", "Cube method (sampling), $\\mathcal S_4$",
             "Local pivotal method, planar distances", "Local pivotal method, wrapped distances"]
    order = [x for x in order if x in rows5]
    assert sorted(order) == sorted(rows5) == sorted(rows6), set(order) ^ set(rows5)
    t5 += [rows5[x] for x in order]
    t6 += [rows6[x] for x in order]
    return "\n".join(t5), "\n".join(t6)


def grid_tables_v2(sh):
    """Table 4 (design effects), Table 5 (worst case over Y_1(sigma^2) for S4 and the main-effect weight share tau)
    and Table S5 (the full table: both worst cases, tau and the shares of exactly balanced samples)."""
    t4, tS5 = grid_tables(sh)
    t5 = ["| Design | $W$, $k=61$ / $k=183$ | $\\tau$, $k=61$ / $k=183$ |", "|---|---|---|"]
    for line in tS5.split("\n")[2:]:
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        t5.append(f"| {c[0]} | {c[2]} | {c[3]} |")
    return t4, "\n".join(t5), tS5


def tableS2():
    """Table S2: the cube and local cube methods of BalancedSampling 2.1.1 balanced on S4 (implementation
    diagnostics), with the cube method with QR directions for comparison: sizes, design effects of the
    Horvitz-Thompson mean and of the realised sample mean, tau on S4 and the lower bound (S1) of the supplement
    (Proposition 3(a) for designs of varying size, responses with mean zero) with C = 1."""
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    srs = {r["k"]: r["var"] for r in d["designs"] if r["design"] == "SRS"}
    mc = {}
    for path in glob.glob(f"{OUT}/results/mc_grid_v3/*.csv"):
        for r in rd(path):
            mc.setdefault(r["design"], {})[int(r["k"])] = r
    fx = {int(r["k"]): r for r in rd(f"{OUT}/results/cube_grid_fixed.csv")}
    bc = {(r["design"], int(r["k"])): r for r in rd(f"{OUT}/results/bound_centred.csv") if r["frame"] == "61x61"}
    f4 = lambda x, se: f"{x:.4f} ({se:.4f})"
    rows = ["| Method | $k$ | Draws | Size $\\neq k$ | Sizes | Horvitz–Thompson | Sample mean | $\\tau$ | Bound (S1) |",
            "|---|---|---|---|---|---|---|---|---|"]
    for key, label, bkey in (("cubeS4", "Cube method, package", "cube S4 (package)"),
                             ("lcubeS4", "Local cube method, package", "local cube S4 (package)"),
                             ("fixed", "Cube method with QR directions", "cube S4, QR directions")):
        for k in (61, 183):
            b = bc[(bkey, k)]
            if key == "fixed":
                r = fx[k]
                ht, mean = f4(float(r["deff_ht"]), float(r["deff_ht_se"])), f4(float(r["deff_mean"]), float(r["deff_mean_se"]))
                B, exc = int(r["B"]), int(r["size_exceptions"])
            else:
                r = mc[key][k]
                ht = f4(float(r["mse_ht"]) / srs[k], float(r["mse_ht_se"]) / srs[k])
                mean = f4(float(r["mse"]) / srs[k], float(r["mse_se"]) / srs[k])
                B, exc = int(r["B"]), int(r["size_exceptions"])
            assert int(b["B"]) == B and int(b["size_exceptions"]) == exc
            sizes = f"{b['size_min']} to {b['size_max']}" if b["size_min"] != b["size_max"] else b["size_min"]
            tau, tse = float(b["tau"]), float(b["tau_se"])
            rows.append(f"| {label} | {k} | {B} | {exc} | {sizes} | {ht} | {mean} | {tau:.4f} ({tse:.4f}) | "
                        f"$\\ge{float(b['bound']):.3f}$ ({float(b['bound_se']):.4f}) |")
    return "\n".join(rows)


def table_certificate():
    """Table A.1: the 19 vectors (S_3, ..., S_7) of the balanced samples of 2^7 (one- and two-variable margins,
    k = 16) with the value 3 S_3 + 2 S_4 + S_5 of the dual certificate, which is at least 14 k^2 = 3584."""
    a = [x for x in json.load(open(f"{OUT}/results/attainment.json")) if x.get("k") == 16][0]
    vec = {}
    for p in a["profiles"]:
        vec.setdefault(tuple(p["S"][3:]), set()).add(p["flat"])
    rows = ["| $S_3$ | $S_4$ | $S_5$ | $S_6$ | $S_7$ | $3S_3+2S_4+S_5$ | Samples |", "|---|---|---|---|---|---|---|"]
    for v in sorted(vec):
        kind = {frozenset({1}): "flats", frozenset({0}): "not flats", frozenset({0, 1}): "both"}[frozenset(vec[v])]
        assert sum(v) == 7 * 256 and 3 * v[0] + 2 * v[1] + v[2] >= 14 * 256
        rows.append(f"| {v[0]} | {v[1]} | {v[2]} | {v[3]} | {v[4]} | {3 * v[0] + 2 * v[1] + v[2]} | {kind} |")
    return "\n".join(rows)


def table5():
    """Table 7: exhaustive evaluation of the variance estimators over the support (varest_exact.py and the grid
    results of study_grid.py), with relative bias, coefficient of variation and relative root mean squared error
    sqrt(var(v) + (E v - V)^2)/V of each estimator against the variance V of the design."""
    d = json.load(open(f"{OUT}/results/grid_exact.json"))
    name_s = {"S1 rows": "$\\mathcal S_1$", "S2 rows+cols": "$\\mathcal S_2$", "S3 +anti-diagonals": "$\\mathcal S_3$", "S4 +diagonals": "$\\mathcal S_4$"}
    rel_rmse = lambda E, cv, V: math.sqrt((cv * E) ** 2 + (E - V) ** 2) / V
    pct = lambda b: "0" if abs(b) < 1e-9 else f"${100 * b:.1f}\\%$"

    def rng(vals, fmt):
        a, b = fmt(min(vals)), fmt(max(vals))
        return a if a == b else f"{a} to {b}"
    rows = ["| Population and design | $k$ | Estimator | df | Bias | CV | RMSE | Coverage |",
            "|---|---|---|---|---|---|---|---|"]
    for r in d["variance_estimation"]:
        name = "equal-weight coset design" if r["design"] == "flat" else f"best {r['design'][4:]} lattices"
        rm = rel_rmse(r["E_v"], r["cv_v"], r["true_var"])
        rows.append(f"| Grid, {name}, {name_s[r['strat']]} | 183 | random groups | 2 | {pct(r['rel_bias'])} | {r['cv_v']:.2f} | "
                    f"{rm:.2f} | {100 * r['coverage_t95']:.1f}% |")
    ve = rd(f"{OUT}/results/varest_exact.csv")
    g = [r for r in ve if r["population"] == "61x61 grid" and r["estimator"] == "regression"][0]
    rm = rel_rmse(float(g["E_est"]), float(g["cv"]), float(g["true_var"]))
    rows.append(f"| Grid, equal-weight coset design, $\\mathcal S_1$ | 183 | regression | {g['df_min']} | {pct(float(g['rel_bias']))} | "
                f"{float(g['cv']):.2f} | {rm:.2f} | {100 * float(g['coverage_t95']):.1f}% |")
    for pop, frame in (("2^7", "$2^7$"), ("3^5", "$3^5$")):
        for k in sorted({int(r["k"]) for r in ve if r["population"] == pop}):
            for est in ("random groups", "regression"):
                sub = [r for r in ve if r["population"] == pop and int(r["k"]) == k and r["estimator"] == est]
                if not sub:
                    continue
                cv = [float(r["cv"]) for r in sub]; cov = [100 * float(r["coverage_t95"]) for r in sub]
                bias = [float(r["rel_bias"]) for r in sub]
                rm = [rel_rmse(float(r["E_est"]), float(r["cv"]), float(r["true_var"])) for r in sub]
                df = (sub[0]["df_min"] if sub[0]["df_min"] == sub[0]["df_max"] else f"{sub[0]['df_min']} to {sub[0]['df_max']}")
                bias_s = "0" if max(abs(b) for b in bias) < 1e-9 else f"{pct(min(bias))} to {pct(max(bias))}"
                rows.append(f"| {frame}, $C(h-1,\\nu,q)$, five responses | {k} | {est} | {df} | {bias_s} | "
                            f"{rng(cv, lambda x: f'{x:.2f}')} | {rng(rm, lambda x: f'{x:.2f}')} | "
                            f"{min(cov):.1f}% to {max(cov):.1f}% |")
    return "\n".join(rows)


def tableA1():
    """Replicated populations: median and 5%-95% range of ratios of design effects to the minimax coset design."""
    rows = rd(f"{OUT}/results/replicated_populations.csv")
    get = {(r["q"], r["k"], r["response"], r["design"], r["statistic"]): r for r in rows}
    P = rows[0]["P"]
    designs = [("coset-uniform", "Uniform"), ("coset-replicated", "$r=q$"), ("lattice", "Lattice"),
               ("cube", "Cube"), ("cube-continued", "Cube, continued")]
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
    t4, t5, tS5 = grid_tables_v2(sh)
    # file names follow the numbering of the manuscript, whose Table 1 (prior work) is not computed
    t = {"table2": table1(), "table3": table3_main(), "table4": t4, "table5": t5, "tableA.1": table_certificate(),
         "tableS2": tableS2(), "tableS3": table2(), "tableS4": table3(), "tableS5": tS5}
    if os.path.exists(f"{OUT}/results/replicated_populations.csv"):
        t["tableS1"] = tableA1()[0]
    to = table_order()                                     # supplementary tables are numbered in order of appearance
    if to is not None:
        t["tableS6"], t["tableS7"] = to, table5()
    else:
        t["tableS6"] = table5()
    for stale in ("table6", "table7") + (("tableS7",) if to is None else ()):   # numbering of release v1.8
        if os.path.exists(f"{OUT}/results/{stale}.md"):
            os.remove(f"{OUT}/results/{stale}.md")
    for name, body in t.items():
        open(f"{OUT}/results/{name}.md", "w").write(body + "\n")
        print(f"[{name}]\n{body}\n")
