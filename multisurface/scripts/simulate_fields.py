"""Simulated Gaussian fields on the non-periodic 61 x 61 grid (PLAN.md, surfaces 2).

Matern correlation with smoothness nu in {0.5, 1.5, 2.5}, scaled so that the correlation is 0.05 at the practical
range rho in {5, 15, 45} cells; nugget share in {0, 0.25}; 20 fields per setting, exact simulation through the
eigendecomposition of the covariance matrix. Writes data/simulated_fields.npz.
"""
import itertools
import os

import numpy as np
from scipy.optimize import brentq
from scipy.special import gamma, kv

HERE = os.path.dirname(os.path.abspath(__file__))
q, REPS, SEED = 61, 20, 20261019
ii, jj = np.divmod(np.arange(q * q), q)
D = np.hypot(ii[:, None] - ii[None, :], jj[:, None] - jj[None, :])


def matern(d, nu, phi):
    t = np.where(d > 0, d / phi, 1.0)
    c = 2 ** (1 - nu) / gamma(nu) * t ** nu * kv(nu, t)
    return np.where(d > 0, c, 1.0)


rng = np.random.default_rng(SEED)
fields, settings = [], []
for nu, rho in itertools.product((0.5, 1.5, 2.5), (5, 15, 45)):
    phi = brentq(lambda p: matern(np.array(rho, float), nu, p) - 0.05, 1e-3, 1e3)
    lam, V = np.linalg.eigh(matern(D, nu, phi))
    L = V * np.sqrt(np.clip(lam, 0, None))
    for nug in (0.0, 0.25):
        for _ in range(REPS):
            z = L @ rng.standard_normal(q * q)
            e = rng.standard_normal(q * q)
            fields.append((np.sqrt(1 - nug) * z + np.sqrt(nug) * e).reshape(q, q))
            settings.append((nu, rho, nug))
    print(f"nu={nu} range={rho}: phi={phi:.3f}, smallest eigenvalue {lam.min():.2e}")
np.savez_compressed(os.path.join(HERE, "..", "data", "simulated_fields.npz"), fields=np.array(fields),
                    settings=np.array(settings))
print(len(fields), "fields")
