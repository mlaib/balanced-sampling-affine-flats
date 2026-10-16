/* All samples of size 16 of F_2^7 balanced on every one- and two-variable margin (strength-2 orthogonal
 * arrays OA(16,7,2,2) with distinct runs), up to permutations of runs, for settle_attainment.py.
 * Columns are +-1 vectors of length 16.  By a permutation of the runs, column 1 is (+^8,-^8) and column 2
 * is (+^4,-^4,+^4,-^4).  A third column orthogonal to both has t, 4-t, 4-t and t entries +1 in the four
 * blocks of four runs defined by columns 1 and 2; changing its sign (a translation, which leaves the
 * profile unchanged) gives t <= 2, and permuting the runs within blocks puts its +1 entries first in each
 * block.  Columns 4-7 are then 4 pairwise orthogonal balanced vectors orthogonal to columns 1-3, taken as
 * an increasing 4-tuple (the profile below does not depend on the order of the columns).  Usage: enum_oa16 t  For every such array with distinct runs, J_U = sum over runs of the product
 * of the columns in U, and S_j = sum_{|U|=j} J_U^2; the run set is a flat iff it is closed under
 * addition of the 7-bit run codes.  Output: one line per distinct (flat, S_0..S_7) with its count. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define R 16
#define MAXC 2000
static uint16_t cand[MAXC]; static int nc = 0;
static uint64_t adj[MAXC][(MAXC + 63) / 64];
static int W;
typedef struct { long long S[8]; int flat; long long count; } Rec;
static Rec recs[100000]; static int nrec = 0;
static long long ntuples = 0, ndistinct_runs = 0;

static int dot(uint16_t a, uint16_t b) { return 16 - 2 * __builtin_popcount((a ^ b) & 0xFFFF); }

static void record(uint16_t col[7]) {
    int code[R];
    for (int r = 0; r < R; r++) {
        int c = 0;
        for (int i = 0; i < 7; i++) if (col[i] >> r & 1) c |= 1 << i;   /* bit set = level 1 (value -1) */
        code[r] = c;
    }
    char seen[128]; memset(seen, 0, sizeof seen);
    for (int r = 0; r < R; r++) { if (seen[code[r]]) return; seen[code[r]] = 1; }
    ndistinct_runs++;
    long long S[8] = {0};
    for (int u = 0; u < 128; u++) {
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

int main(int argc, char **argv) {
    int tt = argc > 1 ? atoi(argv[1]) : 0;
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
    uint16_t col[7]; col[0] = c1; col[1] = c2; col[2] = c3;
    uint64_t s2[(MAXC + 63) / 64], s3[(MAXC + 63) / 64];
    for (int a = 0; a < nc; a++) {
        col[3] = cand[a];
        for (int b = a + 1; b < nc; b++) {
            if (!(adj[a][b / 64] >> (b % 64) & 1)) continue;
            col[4] = cand[b];
            for (int w = 0; w < W; w++) s2[w] = adj[a][w] & adj[b][w];
            for (int c = b + 1; c < nc; c++) {
                if (!(s2[c / 64] >> (c % 64) & 1)) continue;
                col[5] = cand[c];
                for (int w = 0; w < W; w++) s3[w] = s2[w] & adj[c][w];
                for (int d = c + 1; d < nc; d++) {
                    if (!(s3[d / 64] >> (d % 64) & 1)) continue;
                    col[6] = cand[d];
                    ntuples++;
                    record(col);
                }
            }
        }
        if (a % 100 == 0) fprintf(stderr, "  a=%d/%d tuples=%lld profiles=%d\n", a, nc, ntuples, nrec);
    }
    printf("# t=%d tuples=%lld with_distinct_runs=%lld distinct_profiles=%d\n", tt, ntuples, ndistinct_runs, nrec);
    for (int r = 0; r < nrec; r++) {
        printf("%lld %d", recs[r].count, recs[r].flat);
        for (int j = 0; j < 8; j++) printf(" %lld", recs[r].S[j]);
        printf("\n");
    }
    return 0;
}
