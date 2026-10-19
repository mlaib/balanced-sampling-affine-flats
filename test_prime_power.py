"""Supplement S1: a prime power reduces to its prime field (Section 4).

GF(4)^2 is the additive group GF(2)^4, with a + b w (w^2 = w + 1) coded as the bits (a, b). Two or three of the
five GF(4)-lines through 0 define the stratifications, as for rows and columns, and one diagonal, of a grid of side 4;
each is a GF(2)-subspace of dimension 2 with the same cosets. The test lists the GF(2)-subspaces of dimension 2 that
are transversal to all of them, checks that every coset of each is balanced, counts how many are GF(4)-lines, and
checks the variance of the coset designs (r = 1, 2) against enumeration of their support.
"""
import itertools
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import Space, CosetDesign, admissible_directions, nullspace_mod, brute_force_moments, is_balanced

F4 = [(0, 0), (1, 0), (0, 1), (1, 1)]                       # 0, 1, w, w + 1 = w^2


def mul(x, y):                                              # (a + b w)(c + d w), w^2 = w + 1
    a, b = x
    c, d = y
    return ((a * c + b * d) % 2, (a * d + b * c + b * d) % 2)


def gf4_line(slope):                                        # {(t, slope t)}, or the second axis for slope None
    pts = [(t, mul(slope, t)) if slope is not None else ((0, 0), t) for t in F4]
    return np.array([[*u, *v] for u, v in pts])


sp = Space(2, 4)
lines = {"first axis": gf4_line((0, 0)), "second axis": gf4_line(None), "diagonal": gf4_line((1, 0)),
         "w-line": gf4_line((0, 1)), "w2-line": gf4_line((1, 1))}


def is_gf4_linear(D):
    V = [v for v in itertools.product(range(2), repeat=4) if not (np.array(D) @ v % 2).any()]
    Vs = {tuple(v) for v in V}
    return all((*mul((v[0], v[1]), (0, 1)), *mul((v[2], v[3]), (0, 1))) in Vs for v in V)


rng = np.random.default_rng(4)
for names in (("first axis", "second axis"), ("first axis", "second axis", "diagonal")):
    S = [nullspace_mod(lines[n], 2) for n in names]         # bases of W^perp over GF(2)
    A = admissible_directions(sp, S, 2)                     # dual bases D, flats V = D^perp of size 4
    n_lin = sum(is_gf4_linear(D) for D in A)
    for D in A:
        for pr, s in CosetDesign(sp, [D], r=1).samples():
            assert is_balanced(sp, s, S)
    for r in (1, 2):
        d = CosetDesign(sp, A, r=r)
        for _ in range(20):
            y = rng.standard_normal(sp.N)
            assert abs(brute_force_moments(d, y)[1] - d.variance(y)) < 1e-12
    print(f"GF(4)^2 as GF(2)^4, {len(names)} GF(4)-line stratifications: {len(A)} admissible GF(2)-directions of "
          f"dimension 2, {n_lin} of them GF(4)-lines; every coset balanced; variances of C(2, uniform, r), r = 1, 2, "
          f"match enumeration for 20 responses each")
