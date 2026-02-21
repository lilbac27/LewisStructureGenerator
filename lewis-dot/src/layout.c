#include "layout.h"

#include <stdlib.h>
#include <string.h>

/* sin/cos lookup table (fixed-point, 8-bit fraction) for 12 radial positions. */
static const int16_t cos_tbl[12] = {
     256,  222,  128,    0, -128, -222,
    -256, -222, -128,    0,  128,  222
};

static const int16_t sin_tbl[12] = {
       0,  128,  222,  256,  222,  128,
       0, -128, -222, -256, -222, -128
};

/* Layout helper: render path-like molecules in a straight horizontal line. */
bool layout_linear_chain(const Molecule *mol, const LewisStructure *ls, int ax[MAX_ATOMS], int ay[MAX_ATOMS])
{
    if (mol->num_atoms < 3) return false;
    if (ls->num_bonds != mol->num_atoms - 1) return false;

    uint8_t deg[MAX_ATOMS];
    int8_t neigh[MAX_ATOMS][2];
    memset(deg, 0, sizeof(deg));
    memset(neigh, -1, sizeof(neigh));

    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        uint8_t a = ls->bonds[b].a;
        uint8_t c = ls->bonds[b].b;
        if (a >= mol->num_atoms || c >= mol->num_atoms) return false;
        if (deg[a] >= 2 || deg[c] >= 2) return false;
        neigh[a][deg[a]++] = (int8_t)c;
        neigh[c][deg[c]++] = (int8_t)a;
    }

    int endpoints = 0;
    int start = -1;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        if (deg[i] == 1) {
            endpoints++;
            if (start < 0) start = i;
        } else if (deg[i] != 2) {
            return false;
        }
    }
    if (endpoints != 2 || start < 0) return false;

    uint8_t order[MAX_ATOMS];
    int prev = -1;
    int cur = start;
    for (uint8_t k = 0; k < mol->num_atoms; k++) {
        order[k] = (uint8_t)cur;
        int next = -1;
        for (uint8_t ni = 0; ni < deg[cur]; ni++) {
            int cand = neigh[cur][ni];
            if (cand != prev) {
                next = cand;
                break;
            }
        }
        prev = cur;
        cur = next;
        if (k + 1 < mol->num_atoms && cur < 0) return false;
    }

    int step = (SCR_W - 80) / (mol->num_atoms - 1);
    if (step > BOND_LEN) step = BOND_LEN;
    if (step < 22) step = 22;
    int total_w = step * (mol->num_atoms - 1);
    int x0 = LEWIS_CENTER_X - total_w / 2;

    for (uint8_t k = 0; k < mol->num_atoms; k++) {
        uint8_t idx = order[k];
        ax[idx] = x0 + k * step;
        ay[idx] = LEWIS_CENTER_Y;
    }
    return true;
}

/* Layout helper: place atoms by graph distance from central atom. */
bool layout_tree_from_central(const Molecule *mol, const LewisStructure *ls, int ax[MAX_ATOMS], int ay[MAX_ATOMS])
{
    if (mol->num_atoms == 0) return false;

    int8_t dist[MAX_ATOMS];
    int8_t parent[MAX_ATOMS];
    uint8_t q[MAX_ATOMS];
    uint8_t qh = 0;
    uint8_t qt = 0;

    for (uint8_t i = 0; i < MAX_ATOMS; i++) {
        dist[i] = -1;
        parent[i] = -1;
    }

    dist[mol->central] = 0;
    q[qt++] = mol->central;

    while (qh < qt) {
        uint8_t u = q[qh++];
        for (uint8_t b = 0; b < ls->num_bonds; b++) {
            int v = -1;
            if (ls->bonds[b].a == u) v = ls->bonds[b].b;
            else if (ls->bonds[b].b == u) v = ls->bonds[b].a;
            if (v < 0 || v >= mol->num_atoms) continue;
            if (dist[(uint8_t)v] != -1) continue;
            dist[(uint8_t)v] = dist[u] + 1;
            parent[(uint8_t)v] = (int8_t)u;
            q[qt++] = (uint8_t)v;
        }
    }

    int max_dist = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        if (dist[i] < 0) return false;
        if (dist[i] > max_dist) max_dist = dist[i];
    }

    ax[mol->central] = LEWIS_CENTER_X;
    ay[mol->central] = LEWIS_CENTER_Y;

    /* First shell around central */
    uint8_t first[MAX_ATOMS];
    uint8_t n_first = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        if (dist[i] == 1) first[n_first++] = i;
    }
    for (uint8_t k = 0; k < n_first; k++) {
        int angle_idx = (k * 12) / (n_first ? n_first : 1);
        uint8_t node = first[k];
        ax[node] = LEWIS_CENTER_X + (int)(cos_tbl[angle_idx] * BOND_LEN / 256);
        ay[node] = LEWIS_CENTER_Y + (int)(sin_tbl[angle_idx] * BOND_LEN / 256);
    }

    /* Outer shells extend away from central, with slight sibling spreading. */
    for (int d = 2; d <= max_dist; d++) {
        for (uint8_t i = 0; i < mol->num_atoms; i++) {
            if (dist[i] != d) continue;
            int p = parent[i];
            if (p < 0) return false;

            int dx = ax[p] - LEWIS_CENTER_X;
            int dy = ay[p] - LEWIS_CENTER_Y;
            if (dx == 0 && dy == 0) dx = 1;

            int len = abs(dx) > abs(dy) ? abs(dx) : abs(dy);
            if (len == 0) len = 1;

            int bx = ax[p] + dx * BOND_LEN / len;
            int by = ay[p] + dy * BOND_LEN / len;

            int sib_count = 0;
            int sib_idx = 0;
            for (uint8_t j = 0; j < mol->num_atoms; j++) {
                if (dist[j] == d && parent[j] == p) {
                    if (j == i) sib_idx = sib_count;
                    sib_count++;
                }
            }

            if (sib_count > 1) {
                int pdx = -dy;
                int pdy = dx;
                int plen = abs(pdx) > abs(pdy) ? abs(pdx) : abs(pdy);
                if (plen == 0) plen = 1;
                int spread = (sib_idx * 2 - (sib_count - 1)) * 8;
                bx += pdx * spread / plen;
                by += pdy * spread / plen;
            }

            ax[i] = bx;
            ay[i] = by;
        }
    }

    return true;
}
