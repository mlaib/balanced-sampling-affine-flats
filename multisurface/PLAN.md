# Multi-surface test of the best-lattice coset designs (v3, Section 5.2)

Written 2026-10-19, before any result was computed (then in a folder outside this repository; moved here unchanged, with the outcome appended, for version 4 of the manuscript).

## Question

On the 61x61 volcano block, the coset designs mixing the best L in {4, 8, 16} lattice directions admissible for
S4 (rows, columns, both wrapped diagonals) had design effects 0.010-0.012, against 0.041-0.098 for the local
pivotal method and 0.089-0.31 for the cube methods (Table 3). The paper calls the gain specific to that population.
Does it hold on other surfaces?

## Designs (all fixed before running; k = 61 and k = 183, equal inclusion probabilities)

Exact (spectral formula, `affsamp.py`):
- equal-weight coset designs for S2 and S4;
- best-L coset designs for S4, L = 4, 8, 16, with the ranking of `study_grid.py` (shortest wrapped dual vector,
  ties by increasing slope). No other selection rule is tried.

Monte Carlo (BalancedSampling 2.1.1; the same stored samples are evaluated on every surface):
- local pivotal method `lpm2`, planar coordinates;
- local cube method `lcube`, spread on planar coordinates and balanced on rows and columns (S2), the strongest
  spread-plus-balance competitor (S2 is totally unimodular, so the package's row-reduction defect does not occur;
  this is checked by recording sizes and stratum counts);
- cube method `cube` balanced on S2.

## Surfaces

1. Real elevation: public Terrarium terrain tiles (AWS `elevation-tiles-prod`), zoom 10 (pixel about
   153 m x cos(latitude), coarser than the 90 m SRTM source, so no upsampling). For each of 16 locations fixed
   below, the tile containing the point, cut into its 16 non-overlapping 61x61 blocks (offsets 6 + 61 i,
   i = 0..3). A block is excluded only if some pixel is below -50 m (sea bathymetry).
   Locations: Luxembourg Ardennes (50.05, 6.03); Swiss Alps (46.02, 7.75); Flevoland (52.55, 5.60);
   Iowa farmland (42.00, -93.50); Grand Canyon (36.10, -112.10); Everest (27.99, 86.93);
   Erg Chebbi dunes (31.15, -3.97); Pennsylvania ridge-and-valley (40.60, -77.60); Cairngorms (57.08, -3.67);
   Ethiopian highlands (9.00, 38.70); Amazon lowland (-3.50, -62.00); Tibetan plateau (33.00, 88.00);
   Black Forest (48.00, 8.00); Mount Fuji (35.36, 138.73); Loess Plateau (36.50, 109.50);
   Central Australia (-25.00, 131.00).
2. Simulated Gaussian fields on the non-periodic 61x61 grid (exact Cholesky): Matern smoothness
   nu in {0.5, 1.5, 2.5}, practical range in {5, 15, 45} cells, nugget share in {0, 0.25};
   20 fields per setting (360 fields), seed fixed.
3. Sanity check: the volcano block of the paper must reproduce Table 3.

## Decision rule (fixed in advance)

Primary design: best-16 (the most conservative lattice design, worst case 3.62 x equal weight).
Competitor per surface and k: the better of lpm2 and lcube-S2.
On the real elevation blocks, for each k:
- **holds**: best-16 beats the competitor on at least 80% of blocks AND the median ratio
  deff(best-16) / deff(competitor) is at most 0.5;
- **partial**: best-16 beats the competitor on more than 50% of blocks, but the holds condition fails;
- **fails**: otherwise.
The simulated fields show where (smoothness, range, nugget) any gain comes from; they do not change the verdict.

## Outcome (2026-10-19, after the runs)

See results/summary.md. Pre-registered verdict: k = 61 **holds** (best-16 better than the better of lpm2 and
lcube-S2 on 98.4% of the 256 real blocks, median ratio 0.468); k = 183 **partial** (60.9%, median 0.915).
Exploratory, post hoc, same surfaces (scripts/explore_systematic_cosets.py): choosing the 3 cosets systematically
(labels c, c+20, c+40) gives 86.3% wins and median 0.620 at k = 183, with worst case 11.2 x equal weight.
