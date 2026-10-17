/* All samples of size 16 of F_2^m (m = 7 in the paper; 4 <= m <= 7) balanced on every one- and two-variable
 * margin (strength-2 orthogonal arrays OA(16,m,2,2) with distinct runs), up to permutations of the runs and
 * of columns 4 to m, for settle_attainment.py.
 * Columns are +-1 vectors of length 16.  By a permutation of the runs, column 1 is (+^8,-^8) and column 2
 * is (+^4,-^4,+^4,-^4).  A third column orthogonal to both has t, 4-t, 4-t and t entries +1 in the four
 * blocks of four runs defined by columns 1 and 2, and permuting the runs within blocks puts its +1 entries
 * first in each block.  settle_attainment.py uses t = 1 and t = 2 (see there); t = 0, 3 and 4 serve the
 * cross-checks of crosscheck_oa16.py.  Columns 4 to m are then m-3 pairwise orthogonal balanced vectors
 * orthogonal to columns 1-3, taken as an increasing (m-3)-tuple (the profile below does not depend on the
 * order of the columns); no other pruning is done.  Usage: enum_oa16 t [m]  For every such array with
 * distinct runs, J_U = sum over runs of the product of the columns in U, and S_j = sum_{|U|=j} J_U^2; the
 * run set is a flat iff it is closed under addition of the m-bit run codes.  Output: the number of
 * candidate columns (stderr), then one line per distinct (flat, S_0..S_m) with its count. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define R 16
#define MAXC 2000
static uint16_t cand[MAXC]; static int nc = 0;
static uint64_t adj[MAXC][(MAXC + 63) / 64];
static int W, M = 7;
typedef struct { long long S[8]; int flat; long long count; } Rec;
static Rec recs[100000]; static int nrec = 0;
static long long ntuples = 0, ndistinct_runs = 0;
static uint16_t col[7];

static int dot(uint16_t a, uint16_t b) { return 16 - 2 * __builtin_popcount((a ^ b) & 0xFFFF); }

static void record(void) {
    int code[R];
    for (int r = 0; r < R; r++) {
        int c = 0;
        for (int i = 0; i < M; i++) if (col[i] >> r & 1) c |= 1 << i;   /* bit set = level 1 (value -1) */
        code[r] = c;
    }
    char seen[128]; memset(seen, 0, sizeof seen);
    for (int r = 0; r < R; r++) { if (seen[code[r]]) return; seen[code[r]] = 1; }
    ndistinct_runs++;
    long long S[8] = {0};
    for (int u = 0; u < (1 << M); u++) {
        int J = 0;
        for (int r = 0; r < R; r++) J += (__builtin_popcount(u & code[r]) & 1) ? -1 : 1;
        S[__builtin_popcount(u)] += (long long)J * J;
    }
    int flat = 1;
    for (int a = 0; a < R && flat; a++)
        for (int b = a + 1; b < R; b++) {
            int c = code[0] ^ code[a] ^ code[b];          /* translate so that run 0 is the origin */
            if (!seen[c]) { flat = 0; break; }
        }
    for (int t = 0; t < nrec; t++)
        if (recs[t].flat == flat && memcmp(recs[t].S, S, sizeof S) == 0) { recs[t].count++; return; }
    if (nrec >= 100000) { fprintf(stderr, "too many profiles\n"); exit(1); }
    memcpy(recs[nrec].S, S, sizeof S); recs[nrec].flat = flat; recs[nrec].count = 1; nrec++;
}

/* columns depth..M-1 from the candidates after index last that are orthogonal to all chosen columns */
static void extend(int depth, int last, const uint64_t *allowed) {
    if (depth == M) { ntuples++; record(); return; }
    uint64_t next[(MAXC + 63) / 64];
    for (int a = last + 1; a < nc; a++) {
        if (!(allowed[a / 64] >> (a % 64) & 1)) continue;
        col[depth] = cand[a];
        for (int w = 0; w < W; w++) next[w] = allowed[w] & adj[a][w];
        extend(depth + 1, a, next);
        if (depth == 3 && a % 100 == 0) fprintf(stderr, "  a=%d/%d tuples=%lld profiles=%d\n", a, nc, ntuples, nrec);
    }
}

int main(int argc, char **argv) {
    int tt = argc > 1 ? atoi(argv[1]) : 0;
    if (argc > 2) M = atoi(argv[2]);
    if (tt < 0 || tt > 4 || M < 4 || M > 7) { fprintf(stderr, "usage: enum_oa16 t [m], 0 <= t <= 4, 4 <= m <= 7\n"); return 1; }
    uint16_t c1 = 0xFF00, c2 = 0xF0F0;        /* bit r set = value -1 in run r */
    int plus[4] = {tt, 4 - tt, 4 - tt, tt};   /* +1 entries of column 3 in blocks rows 0-3, 4-7, 8-11, 12-15 */
    uint16_t c3 = 0;
    for (int bl = 0; bl < 4; bl++)
        for (int i = plus[bl]; i < 4; i++) c3 |= (uint16_t)(1u << (4 * bl + i));
    if (__builtin_popcount(c3) != 8 || dot(c3, c1) != 0 || dot(c3, c2) != 0) { fprintf(stderr, "bad c3\n"); return 1; }
    for (int v = 0; v < 65536; v++) {
        if (__builtin_popcount(v) != 8) continue;
        if (dot(v, c1) != 0 || dot(v, c2) != 0 || dot(v, c3) != 0) continue;
        cand[nc++] = (uint16_t)v;
    }
    W = (nc + 63) / 64;
    for (int a = 0; a < nc; a++)
        for (int b = 0; b < nc; b++)
            if (a != b && dot(cand[a], cand[b]) == 0) adj[a][b / 64] |= 1ULL << (b % 64);
    fprintf(stderr, "t=%d candidates: %d\n", tt, nc);
    col[0] = c1; col[1] = c2; col[2] = c3;
    uint64_t all[(MAXC + 63) / 64];
    memset(all, 0, sizeof all);
    for (int a = 0; a < nc; a++) all[a / 64] |= 1ULL << (a % 64);
    extend(3, -1, all);
    if (M == 7) printf("# t=%d tuples=%lld with_distinct_runs=%lld distinct_profiles=%d\n", tt, ntuples, ndistinct_runs, nrec);
    else printf("# t=%d m=%d tuples=%lld with_distinct_runs=%lld distinct_profiles=%d\n", tt, M, ntuples, ndistinct_runs, nrec);
    for (int r = 0; r < nrec; r++) {
        printf("%lld %d", recs[r].count, recs[r].flat);
        for (int j = 0; j <= M; j++) printf(" %lld", recs[r].S[j]);
        printf("\n");
    }
    return 0;
}
