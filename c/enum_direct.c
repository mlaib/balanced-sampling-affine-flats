/* Independent check of enum_oa16.c (crosscheck_oa16.py): all samples of size 16 of F_2^m (4 <= m <= 7)
 * that contain the origin and are balanced on every one- and two-variable margin, by backtracking over the
 * units in increasing order of their m-bit codes, with the count of every margin cell as a constraint and
 * without any reduction by symmetries.  A unit is added only if no one-variable cell then holds more than
 * 8 units and no two-variable cell more than 4; a branch is abandoned when some cell can no longer be
 * filled by the units that remain.  Each sample is also classified by tau = min(a, 4 - a), where a is the
 * number of its units with x_1 = x_2 = x_3 = 0: the third column of its array has a, 4-a, 4-a and a entries
 * +1 in the four blocks defined by the first two, as the parameter t of enum_oa16, and a translation
 * changes a into a or 4 - a.  Usage: enum_direct m.  Output: the number of samples, then one line per
 * distinct (flat, tau, S_0..S_m) with its count. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define K 16
static int M, N, NCELL;
static int cells[128][28], ncells_of;       /* cells of unit x: m one-variable and C(m,2) two-variable */
static int cap[98], cnt[98], avail[98][129]; /* avail[c][x]: units y >= x in cell c */
static int cur[K], len = 0;
typedef struct { long long S[8]; int flat, tau; long long count; } Rec;
static Rec recs[100000]; static int nrec = 0;
static long long nsol = 0;

static void record(void) {
    char seen[128]; memset(seen, 0, sizeof seen);
    for (int r = 0; r < K; r++) seen[cur[r]] = 1;
    long long S[8] = {0};
    for (int u = 0; u < N; u++) {
        int J = 0;
        for (int r = 0; r < K; r++) J += (__builtin_popcount(u & cur[r]) & 1) ? -1 : 1;
        S[__builtin_popcount(u)] += (long long)J * J;
    }
    int flat = 1;                                  /* cur[0] = 0: a flat iff closed under addition */
    for (int a = 0; a < K && flat; a++)
        for (int b = a + 1; b < K; b++)
            if (!seen[cur[a] ^ cur[b]]) { flat = 0; break; }
    int a = 0;
    for (int r = 0; r < K; r++) a += (cur[r] & 7) == 0;
    int tau = a < 4 - a ? a : 4 - a;
    nsol++;
    for (int t = 0; t < nrec; t++)
        if (recs[t].flat == flat && recs[t].tau == tau && memcmp(recs[t].S, S, sizeof S) == 0) { recs[t].count++; return; }
    if (nrec >= 100000) { fprintf(stderr, "too many profiles\n"); exit(1); }
    memcpy(recs[nrec].S, S, sizeof S); recs[nrec].flat = flat; recs[nrec].tau = tau; recs[nrec].count = 1; nrec++;
}

static int feasible(int next) {                    /* can every cell still be filled from units >= next? */
    for (int c = 0; c < NCELL; c++)
        if (cap[c] - cnt[c] > avail[c][next]) return 0;
    return 1;
}

static void rec(int last) {
    if (len == K) { record(); return; }
    for (int x = last + 1; x <= N - (K - len); x++) {
        int ok = 1;
        for (int i = 0; i < ncells_of; i++) if (cnt[cells[x][i]] >= cap[cells[x][i]]) { ok = 0; break; }
        if (!ok) continue;
        for (int i = 0; i < ncells_of; i++) cnt[cells[x][i]]++;
        cur[len++] = x;
        if (feasible(x + 1)) rec(x);
        len--;
        for (int i = 0; i < ncells_of; i++) cnt[cells[x][i]]--;
        if (len == 1) fprintf(stderr, "  second unit %d done, samples %lld\n", x, nsol);
    }
}

int main(int argc, char **argv) {
    M = argc > 1 ? atoi(argv[1]) : 7;
    if (M < 4 || M > 7) { fprintf(stderr, "usage: enum_direct m, 4 <= m <= 7\n"); return 1; }
    N = 1 << M;
    NCELL = 0;
    int id1[7][2], id2[7][7][4];
    for (int i = 0; i < M; i++) for (int b = 0; b < 2; b++) { id1[i][b] = NCELL; cap[NCELL++] = K / 2; }
    for (int i = 0; i < M; i++) for (int j = i + 1; j < M; j++)
        for (int c = 0; c < 4; c++) { id2[i][j][c] = NCELL; cap[NCELL++] = K / 4; }
    for (int x = 0; x < N; x++) {
        int n = 0;
        for (int i = 0; i < M; i++) cells[x][n++] = id1[i][x >> i & 1];
        for (int i = 0; i < M; i++) for (int j = i + 1; j < M; j++)
            cells[x][n++] = id2[i][j][2 * (x >> i & 1) + (x >> j & 1)];
        ncells_of = n;
    }
    for (int c = 0; c < NCELL; c++) avail[c][N] = 0;
    for (int x = N - 1; x >= 0; x--)
        for (int c = 0; c < NCELL; c++) {
            int in = 0;
            for (int i = 0; i < ncells_of; i++) if (cells[x][i] == c) in = 1;
            avail[c][x] = avail[c][x + 1] + in;
        }
    for (int i = 0; i < ncells_of; i++) cnt[cells[0][i]]++;   /* the origin */
    cur[len++] = 0;
    rec(0);
    printf("# m=%d samples_containing_0=%lld distinct_classes=%d\n", M, nsol, nrec);
    for (int r = 0; r < nrec; r++) {
        printf("%lld %d %d", recs[r].count, recs[r].flat, recs[r].tau);
        for (int j = 0; j <= M; j++) printf(" %lld", recs[r].S[j]);
        printf("\n");
    }
    return 0;
}
