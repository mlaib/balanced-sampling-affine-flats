"""Monte Carlo uncertainty of the cube quantiles of Table S1 (Supplementary material S3).

replicate_populations.py evaluates the cube method of BalancedSampling and the transcription with the flight
continued on one bank of 20 000 stored draws per configuration, the same bank for every population, so their
design effects in the 500 populations are conditional on that bank.  Here the bank is resampled: each bootstrap
replicate draws 20 000 whole draws with replacement from the bank (multinomial counts) and uses the same
resampled bank for every population.  For each replicate we recompute, for every configuration, response and
cube variant, the ratio of the design effect to that of the minimax coset design (exact) in each population,
and its 5%, 50% and 95% quantiles and share below 1 over the populations.  The standard deviation of these
statistics over the replicates is the Monte Carlo standard error that the shared bank adds to the quantiles of
Table S1; it is reported separately from the between-population spread.  The script also records the Monte Carlo
standard error of the ratio within each population (its standard deviation over the replicates), summarised by
its median over the populations, and the standard deviation of the ratio across the populations, so that the
two can be compared.
Writes results/bootstrap_bank.csv.  Usage: python3 bootstrap_bank.py [replicates] [P]
"""
import csv, os, sys
import numpy as np
from affsamp import Space
from study_factorial import responses
from replicate_populations import FRAMES, read_samples

OUT = os.path.dirname(os.path.abspath(__file__))
STATS = ("q05", "q50", "q95", "share_below_1")


def statistics(ratio):
    """ratio: replicates x populations -> dict of arrays (replicates)."""
    return {"q05": np.quantile(ratio, 0.05, axis=1), "q50": np.quantile(ratio, 0.5, axis=1),
            "q95": np.quantile(ratio, 0.95, axis=1), "share_below_1": (ratio < 1).mean(axis=1)}


if __name__ == "__main__":
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    P = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    store = np.load(f"{OUT}/results/replicated_populations.npz")
    rng = np.random.default_rng(20261017)
    rows = []
    for q, m, ks in FRAMES:
        sp = Space(q, m)
        N = sp.N
        pops = [responses(sp, 1000 + i) for i in range(P)]
        names = list(pops[0])
        Yall = np.array([[pop[n] for n in names] for pop in pops])            # P x responses x N
        ybar = Yall.mean(axis=2)
        cube = read_samples(f"{OUT}/results/cube_samples/factorial_p{q}_m{m}.txt.gz", N)
        cont = read_samples(f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz", N, select="continued")
        for k in ks:
            vsrs = (1 / k - 1 / N) * Yall.var(axis=2, ddof=1)                  # P x responses
            minimax = store[f"p{q}_k{k}_coset-minimax"]                        # P x responses (exact)
            for name, mats in (("cube", cube), ("cube-continued", cont)):
                S = mats[k]                                                    # draws x N, entries 1/k
                B = len(S)
                counts = rng.multinomial(B, np.full(B, 1.0 / B), size=R).astype(float)   # R x B, shared by all populations
                for j, yname in enumerate(names):
                    E2 = (S @ Yall[:, j, :].T - ybar[:, j][None, :]) ** 2      # draws x P
                    deff0 = E2.mean(axis=0) / vsrs[:, j]
                    assert np.allclose(deff0, store[f"p{q}_k{k}_{name}"][:, j], rtol=1e-10, atol=0)
                    deffb = (counts @ E2) / B / vsrs[:, j][None, :]           # R x P
                    st0 = statistics((deff0 / minimax[:, j])[None, :])
                    stb = statistics(deffb / minimax[:, j][None, :])
                    row = {"q": q, "m": m, "k": k, "response": yname, "design": name, "P": P, "B": B, "replicates": R}
                    for s in STATS:
                        row[s] = float(st0[s][0])
                        row[s + "_se"] = float(stb[s].std(ddof=1))
                    ratiob = deffb / minimax[:, j][None, :]
                    row["ratio_se_median"] = float(np.median(ratiob.std(axis=0, ddof=1)))
                    row["ratio_sd_between"] = float((deff0 / minimax[:, j]).std(ddof=1))
                    rows.append(row)
        print(f"{q}^{m} done", file=sys.stderr, flush=True)
    keys = list(rows[0])
    with open(f"{OUT}/results/bootstrap_bank.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({a: (f"{b:.6g}" if isinstance(b, float) else b) for a, b in r.items()})
    for r in rows:
        print(r["q"], r["k"], f'{r["response"]:12s} {r["design"]:15s}',
              " ".join(f'{s} {r[s]:.3f} ({r[s + "_se"]:.3f})' for s in STATS),
              "within-pop se %.4f between-pop sd %.4f" % (r["ratio_se_median"], r["ratio_sd_between"]))
