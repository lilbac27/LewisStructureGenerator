#include "lewis_engine.h"

#include <string.h>

typedef struct {
    uint8_t valence_pairs;
    uint8_t bond_pairs;
    uint8_t lone_pairs;
    const char *ep_geometry;
    const char *shape;
    const char *hybridization;
    const char *bond_angle;
} VseprRow;

static const VseprRow vsepr_rows[] = {
    {1, 1, 0, "Linear", "Linear", "s", "180"},
    {2, 2, 0, "Linear", "Linear", "sp", "180"},
    {2, 1, 1, "Linear", "Linear", "sp", "180"},
    {3, 3, 0, "Trigonal Planar", "Trigonal Planar", "sp2", "120"},
    {3, 2, 1, "Trigonal Planar", "Bent", "sp2", "<120"},
    {3, 1, 2, "Trigonal Planar", "Linear", "sp2", "180"},
    {4, 4, 0, "Tetrahedral", "Tetrahedral", "sp3", "109.5"},
    {4, 3, 1, "Tetrahedral", "Trigonal Pyramidal", "sp3", "<109.5"},
    {4, 2, 2, "Tetrahedral", "Bent", "sp3", "<109.5"},
    {4, 1, 3, "Tetrahedral", "Linear", "sp3", "180"},
    {5, 5, 0, "Trigonal Bipyramidal", "Trigonal Bipyramidal", "sp3d", "90, 120"},
    {5, 4, 1, "Trigonal Bipyramidal", "Seesaw", "sp3d", "<90, <120"},
    {5, 3, 2, "Trigonal Bipyramidal", "T-shaped", "sp3d", "<90"},
    {5, 2, 3, "Trigonal Bipyramidal", "Linear", "sp3d", "180"},
    {5, 1, 4, "Trigonal Bipyramidal", "Linear", "sp3d", "180"},
    {6, 6, 0, "Octahedral", "Octahedral", "sp3d2", "90"},
    {6, 5, 1, "Octahedral", "Square Pyramidal", "sp3d2", "<90"},
    {6, 4, 2, "Octahedral", "Square Planar", "sp3d2", "90"},
    {6, 3, 3, "Octahedral", "T-shaped", "sp3d2", "<90"},
    {6, 2, 4, "Octahedral", "Linear", "sp3d2", "180"},
    {6, 1, 5, "Octahedral", "Linear", "sp3d2", "180"},
    {7, 7, 0, "Pentagonal Bipyramidal", "Pentagonal Bipyramidal", "sp3d3", "72, 90"},
    {7, 6, 1, "Pentagonal Bipyramidal", "Pentagonal Pyramidal", "sp3d3", "<90"},
    {7, 5, 2, "Pentagonal Bipyramidal", "Pentagonal Planar", "sp3d3", "72"},
    {7, 4, 3, "Pentagonal Bipyramidal", "Seesaw", "sp3d3", "<90, <72"},
    {7, 3, 4, "Pentagonal Bipyramidal", "T-shaped", "sp3d3", "<90"},
    {7, 2, 5, "Pentagonal Bipyramidal", "Linear", "sp3d3", "180"},
    {7, 1, 6, "Pentagonal Bipyramidal", "Linear", "sp3d3", "180"},
};

static int bond_order_sum(const LewisStructure *ls, uint8_t atom_idx)
{
    int sum = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if (ls->bonds[b].a == atom_idx || ls->bonds[b].b == atom_idx) {
            sum += ls->bonds[b].order;
        }
    }
    return sum;
}

static int electrons_on_atom(const LewisStructure *ls, uint8_t atom_idx)
{
    return (ls->lone_pairs[atom_idx] * 2) + (bond_order_sum(ls, atom_idx) * 2);
}

static void recompute_formal_charges(const Molecule *mol, LewisStructure *ls)
{
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        int val = elements[mol->atoms[i].elem].valence;
        int lpe = ls->lone_pairs[i] * 2;
        int bnd_e = bond_order_sum(ls, i);
        ls->formal_charge[i] = (int8_t)(val - lpe - bnd_e);
    }
}

static int formal_charge_sum(const Molecule *mol, const LewisStructure *ls)
{
    int sum = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        sum += ls->formal_charge[i];
    }
    return sum;
}

static bool structures_equal(const Molecule *mol, const LewisStructure *a, const LewisStructure *b)
{
    if (a->num_bonds != b->num_bonds) {
        return false;
    }
    for (uint8_t i = 0; i < a->num_bonds; i++) {
        if (a->bonds[i].a != b->bonds[i].a) return false;
        if (a->bonds[i].b != b->bonds[i].b) return false;
        if (a->bonds[i].order != b->bonds[i].order) return false;
    }
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        if (a->lone_pairs[i] != b->lone_pairs[i]) return false;
    }
    return true;
}

static bool resonance_exists(const Molecule *mol, const LewisStructure *candidate)
{
    for (uint8_t i = 0; i < mol->num_res; i++) {
        if (structures_equal(mol, candidate, &mol->res[i])) {
            return true;
        }
    }
    return false;
}

static bool atom_has_h_neighbor(const Molecule *mol, const LewisStructure *ls, uint8_t atom_idx)
{
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        int other = -1;
        if (ls->bonds[b].a == atom_idx) other = ls->bonds[b].b;
        else if (ls->bonds[b].b == atom_idx) other = ls->bonds[b].a;
        if (other < 0) continue;
        if (mol->atoms[(uint8_t)other].elem == ELEM_H) {
            return true;
        }
    }
    return false;
}

/* Keep protonated terminal oxygens (X-O-H) out of the resonance swap set. */
static bool is_protonated_terminal_oxygen(const Molecule *mol, const LewisStructure *ls, uint8_t term_idx)
{
    if (mol->atoms[term_idx].elem != ELEM_O) return false;
    if (!atom_has_h_neighbor(mol, ls, term_idx)) return false;

    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if ((ls->bonds[b].a == mol->central && ls->bonds[b].b == term_idx) ||
            (ls->bonds[b].b == mol->central && ls->bonds[b].a == term_idx)) {
            return true;
        }
    }
    return false;
}

static int abs_int(int x)
{
    return (x < 0) ? -x : x;
}

static bool center_is_terminal_elem(uint8_t elem_idx)
{
    if (elem_idx == ELEM_H) return true;
    const Element *e = &elements[elem_idx];
    return (e->group == 17 || e->bond_cap <= 1);
}

static void score_structure(const Molecule *mol, const LewisStructure *ls, int *sum_abs_fc, int *nonzero_fc, int *abs_central_fc)
{
    int sum_abs = 0;
    int nonzero = 0;

    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        int fc = ls->formal_charge[i];
        if (fc != 0) {
            nonzero++;
        }
        sum_abs += abs_int(fc);
    }

    *sum_abs_fc = sum_abs;
    *nonzero_fc = nonzero;
    *abs_central_fc = abs_int(ls->formal_charge[mol->central]);
}

static bool candidate_is_better(int cand_sum_abs_fc,
                                int cand_nonzero_fc,
                                int cand_abs_central_fc,
                                uint8_t cand_count,
                                bool cand_terminal,
                                uint8_t cand_eneg,
                                uint8_t cand_period,
                                uint8_t cand_atomic_num,
                                int best_sum_abs_fc,
                                int best_nonzero_fc,
                                int best_abs_central_fc,
                                uint8_t best_count,
                                bool best_terminal,
                                uint8_t best_eneg,
                                uint8_t best_period,
                                uint8_t best_atomic_num)
{
    if (cand_sum_abs_fc != best_sum_abs_fc) return cand_sum_abs_fc < best_sum_abs_fc;
    if (cand_nonzero_fc != best_nonzero_fc) return cand_nonzero_fc < best_nonzero_fc;
    if (cand_abs_central_fc != best_abs_central_fc) return cand_abs_central_fc < best_abs_central_fc;
    if (cand_count != best_count) return cand_count < best_count;

    if (cand_terminal != best_terminal) return !cand_terminal;
    if (cand_eneg != best_eneg) return cand_eneg < best_eneg;
    if (cand_period != best_period) return cand_period < best_period;
    if (cand_atomic_num != best_atomic_num) return cand_atomic_num < best_atomic_num;
    return false;
}

static uint8_t gather_center_candidates(const Molecule *mol, uint8_t out_idx[MAX_ATOMS])
{
    uint8_t n = 0;

    /* Pass 1: non-terminal, non-hydrogen atoms. */
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        uint8_t elem_idx = mol->atoms[i].elem;
        if (elem_idx == ELEM_H) continue;
        if (center_is_terminal_elem(elem_idx)) continue;
        out_idx[n++] = i;
    }

    /* Pass 2: remaining non-hydrogen atoms (halogens/noble-like terminals). */
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        uint8_t elem_idx = mol->atoms[i].elem;
        if (elem_idx == ELEM_H) continue;
        if (!center_is_terminal_elem(elem_idx)) continue;
        out_idx[n++] = i;
    }

    /* All-hydrogen fallback (e.g., H2). */
    if (n == 0 && mol->num_atoms > 0) {
        out_idx[n++] = 0;
    }

    return n;
}

static int required_electrons(const Molecule *mol, uint8_t atom_idx, bool is_central)
{
    uint8_t elem_idx = mol->atoms[atom_idx].elem;
    const Element *e = &elements[elem_idx];

    if (elem_idx == ELEM_H || elem_idx == ELEM_HE) return 2;

    if (is_central && e->group == 2) return 4;   /* Be/Mg */
    if (is_central && e->group == 13) return 6;  /* B/Al */

    return 8;
}

static bool shell_satisfied(const Molecule *mol, uint8_t atom_idx, int electrons, bool is_central)
{
    uint8_t elem_idx = mol->atoms[atom_idx].elem;
    const Element *e = &elements[elem_idx];

    if (elem_idx == ELEM_H || elem_idx == ELEM_HE) return electrons == 2;

    if (electrons < required_electrons(mol, atom_idx, is_central)) {
        return false;
    }

    /* Period 2 atoms should not exceed octet */
    if (e->period <= 2 && electrons > 8) {
        return false;
    }

    return true;
}

static int bond_limit(const Molecule *mol, uint8_t atom_idx, bool is_central)
{
    uint8_t elem_idx = mol->atoms[atom_idx].elem;
    const Element *e = &elements[elem_idx];
    int limit = e->bond_cap;

    /* Allow ammonium-like cations for period-2 group-15 centers (e.g., NH4+). */
    if (is_central && e->period == 2 && e->group == 15 && mol->charge > 0 && limit < 4) {
        limit = 4;
    }

    /* Allow expanded-valence central atoms from period 3+ */
    if (is_central && e->period >= 3) {
        if (e->group == 15 && limit < 5) limit = 5;
        if (e->group == 16 && limit < 6) limit = 6;
        if (e->group == 17 && limit < 7) limit = 7;
    }

    return limit;
}

static bool add_single_bond(LewisStructure *ls, uint8_t a, uint8_t b, int *ve_pool, uint8_t remain[])
{
    if (*ve_pool < 2) return false;
    if (ls->num_bonds >= MAX_BONDS) return false;
    if (remain[a] == 0 || remain[b] == 0) return false;

    ls->bonds[ls->num_bonds].a = a;
    ls->bonds[ls->num_bonds].b = b;
    ls->bonds[ls->num_bonds].order = 1;
    ls->num_bonds++;

    remain[a]--;
    remain[b]--;
    *ve_pool -= 2;
    return true;
}

static bool build_skeleton(const Molecule *mol, LewisStructure *ls, int *ve_pool)
{
    uint8_t remain[MAX_ATOMS];
    bool connected[MAX_ATOMS];
    uint8_t backbone[MAX_ATOMS];
    uint8_t n_backbone = 0;
    uint8_t ordered[MAX_ATOMS];
    bool used[MAX_ATOMS];

    memset(connected, 0, sizeof(connected));
    memset(used, 0, sizeof(used));

    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        uint8_t elem_idx = mol->atoms[i].elem;
        remain[i] = (uint8_t)bond_limit(mol, i, i == mol->central);

        if (i == mol->central) {
            backbone[n_backbone++] = i;
            continue;
        }

        /* Keep highly terminal atoms off the backbone */
        if (elem_idx == ELEM_H) continue;
        if (elements[elem_idx].group == 17) continue;
        if (elements[elem_idx].bond_cap >= 3) {
            backbone[n_backbone++] = i;
        }
    }

    if (mol->num_atoms == 1) {
        return true;
    }

    /*
     * Triatomic special-case: the chosen central atom must connect to both
     * other atoms (A-X-A / A-X-B style), not a three-atom chain.
     */
    if (mol->num_atoms == 3) {
        uint8_t others[2];
        uint8_t n_others = 0;

        connected[mol->central] = true;
        for (uint8_t i = 0; i < mol->num_atoms; i++) {
            if (i == mol->central) continue;
            others[n_others++] = i;
        }

        if (n_others != 2) return false;
        if (!add_single_bond(ls, mol->central, others[0], ve_pool, remain)) return false;
        connected[others[0]] = true;
        if (!add_single_bond(ls, mol->central, others[1], ve_pool, remain)) return false;
        connected[others[1]] = true;
        return true;
    }

    if (n_backbone == 0) {
        return false;
    }

    /* Order backbone: central first, then higher capacity and lower EN */
    ordered[0] = mol->central;
    used[mol->central] = true;
    uint8_t n_ordered = 1;

    while (n_ordered < n_backbone) {
        int best = -1;
        int best_score = -32768;
        for (uint8_t bi = 0; bi < n_backbone; bi++) {
            uint8_t atom = backbone[bi];
            if (used[atom]) continue;
            const Element *e = &elements[mol->atoms[atom].elem];
            int score = (int)bond_limit(mol, atom, false) * 10 - (int)e->eneg;
            if (score > best_score) {
                best_score = score;
                best = atom;
            }
        }
        if (best < 0) return false;
        ordered[n_ordered++] = (uint8_t)best;
        used[(uint8_t)best] = true;
    }

    connected[ordered[0]] = true;
    for (uint8_t i = 1; i < n_ordered; i++) {
        uint8_t a = ordered[i - 1];
        uint8_t b = ordered[i];
        if (!add_single_bond(ls, a, b, ve_pool, remain)) {
            return false;
        }
        connected[a] = true;
        connected[b] = true;
    }

    /* Attach non-backbone atoms in two passes: heavy atoms first, then H. */
    for (uint8_t pass = 0; pass < 2; pass++) {
        bool target_h = (pass == 1);

        for (uint8_t i = 0; i < mol->num_atoms; i++) {
            if (connected[i]) continue;

            bool is_h = (mol->atoms[i].elem == ELEM_H);
            if (is_h != target_h) continue;
            if (remain[i] == 0) return false;

            int best_host = -1;
            int best_score = -32768;
            for (uint8_t j = 0; j < mol->num_atoms; j++) {
                if (!connected[j]) continue;
                if (i == j) continue;
                if (remain[j] == 0) continue;
                if (mol->atoms[j].elem == ELEM_H) continue;

                int score = (int)remain[j] * 10 - (int)elements[mol->atoms[j].elem].eneg;

                /* Count heavy-atom neighbors already attached to this host. */
                int heavy_neighbors = 0;
                for (uint8_t b = 0; b < ls->num_bonds; b++) {
                    int other = -1;
                    if (ls->bonds[b].a == j) other = ls->bonds[b].b;
                    else if (ls->bonds[b].b == j) other = ls->bonds[b].a;
                    if (other >= 0 && mol->atoms[(uint8_t)other].elem != ELEM_H) {
                        heavy_neighbors++;
                    }
                }

                if (!is_h) {
                    if (j == mol->central) score += 12;
                } else {
                    score -= heavy_neighbors * 8;
                    if (j == mol->central) score -= 4;
                }

                if (score > best_score) {
                    best_score = score;
                    best_host = j;
                }
            }

            if (best_host < 0) {
                for (uint8_t j = 0; j < mol->num_atoms; j++) {
                    if (!connected[j]) continue;
                    if (i == j) continue;
                    if (remain[j] == 0) continue;
                    best_host = j;
                    break;
                }
            }

            if (best_host < 0) return false;
            if (!add_single_bond(ls, (uint8_t)best_host, i, ve_pool, remain)) {
                return false;
            }
            connected[i] = true;
        }
    }

    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        if (!connected[i]) return false;
    }
    return true;
}

/*
 * Find central atom with tie-breakers:
 * - avoid H and strongly terminal atoms when possible
 * - then prefer lower electronegativity
 * - then prefer higher bond capacity / frequency
 */
static uint8_t find_central(const Molecule *mol)
{
    if (mol->num_atoms == 0) return 0;

    uint8_t counts[NUM_ELEMENTS];
    memset(counts, 0, sizeof(counts));
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        counts[mol->atoms[i].elem]++;
    }

    int best = -1;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        uint8_t elem_idx = mol->atoms[i].elem;
        if (elem_idx == ELEM_H) continue;

        const Element *cand = &elements[elem_idx];
        bool cand_terminal = (cand->group == 17 || cand->bond_cap <= 1);

        if (best < 0) {
            best = i;
            continue;
        }

        const Element *cur = &elements[mol->atoms[best].elem];
        bool cur_terminal = (cur->group == 17 || cur->bond_cap <= 1);

        if (cand_terminal != cur_terminal) {
            if (!cand_terminal) best = i;
            continue;
        }
        if (cand->eneg < cur->eneg) {
            best = i;
            continue;
        }
        if (cand->eneg == cur->eneg && cand->bond_cap > cur->bond_cap) {
            best = i;
            continue;
        }
        if (cand->eneg == cur->eneg &&
            cand->bond_cap == cur->bond_cap &&
            counts[elem_idx] > counts[mol->atoms[best].elem]) {
            best = i;
        }
    }

    if (best >= 0) return (uint8_t)best;
    return 0;
}

/*
 * Generate one Lewis structure.
 * Returns false when electron count/connectivity/octet constraints are invalid.
 */
static bool generate_structure(const Molecule *mol, LewisStructure *ls, InvalidReason *reason)
{
    memset(ls, 0, sizeof(*ls));
    *reason = INVALID_NONE;

    if (mol->num_atoms == 0) {
        *reason = INVALID_NO_ATOMS;
        return false;
    }
    if (mol->total_ve < 0) {
        *reason = INVALID_NEGATIVE_ELECTRONS;
        return false;
    }
    if ((mol->total_ve & 1) != 0) {
        *reason = INVALID_ODD_ELECTRONS;
        return false;
    }

    int ve_pool = mol->total_ve;

    if (!build_skeleton(mol, ls, &ve_pool)) {
        *reason = INVALID_SKELETON;
        return false;
    }

    /* Fill terminal atoms first */
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        if (i == mol->central) continue;

        int target = required_electrons(mol, i, false);
        int bonded_e = bond_order_sum(ls, i) * 2;
        int need = target - bonded_e;

        if (need > 0) {
            int pairs = need / 2;
            int max_pairs = ve_pool / 2;
            if (pairs > max_pairs) pairs = max_pairs;
            ls->lone_pairs[i] = (uint8_t)pairs;
            ve_pool -= pairs * 2;
        }
    }

    /* Remaining electrons go to central atom */
    if (ve_pool > 0) {
        ls->lone_pairs[mol->central] = (uint8_t)(ve_pool / 2);
        ve_pool -= ls->lone_pairs[mol->central] * 2;
    }

    if (ve_pool != 0) {
        *reason = INVALID_LEFTOVER_ELECTRONS;
        return false;
    }

    /* Promote central bonds to satisfy central shell */
    if (mol->num_atoms > 1) {
        int target_c = required_electrons(mol, mol->central, true);
        uint8_t next_bond = 0;

        for (int pass = 0; pass < MAX_BONDS * 3; pass++) {
            int central_e = electrons_on_atom(ls, mol->central);
            if (central_e >= target_c) break;

            bool promoted = false;
            for (uint8_t scan = 0; scan < ls->num_bonds; scan++) {
                uint8_t b = (next_bond + scan) % ls->num_bonds;
                if (!(ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central)) {
                    continue;
                }

                uint8_t term = (ls->bonds[b].a == mol->central) ? ls->bonds[b].b : ls->bonds[b].a;
                if (mol->atoms[term].elem == ELEM_H) continue;
                if (ls->bonds[b].order >= 3) continue;
                if (ls->lone_pairs[term] == 0) continue;

                ls->bonds[b].order++;
                ls->lone_pairs[term]--;
                next_bond = (b + 1) % ls->num_bonds;
                promoted = true;
                break;
            }
            if (!promoted) break;
        }
    }

    /* For period 3+ centers, use available lone pairs to reduce charge separation. */
    if (elements[mol->atoms[mol->central].elem].period >= 3) {
        recompute_formal_charges(mol, ls);
        for (int pass = 0; pass < MAX_BONDS; pass++) {
            if (ls->formal_charge[mol->central] <= 0) break;

            int best_bond = -1;
            int most_negative = 0;
            for (uint8_t b = 0; b < ls->num_bonds; b++) {
                if (!(ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central)) {
                    continue;
                }

                uint8_t term = (ls->bonds[b].a == mol->central) ? ls->bonds[b].b : ls->bonds[b].a;
                if (mol->atoms[term].elem == ELEM_H) continue;
                if (ls->bonds[b].order >= 3) continue;
                if (ls->lone_pairs[term] == 0) continue;
                if (ls->formal_charge[term] >= 0) continue;
                if (ls->formal_charge[term] < most_negative) {
                    most_negative = ls->formal_charge[term];
                    best_bond = b;
                }
            }

            if (best_bond < 0) break;
            uint8_t term = (ls->bonds[best_bond].a == mol->central) ? ls->bonds[best_bond].b : ls->bonds[best_bond].a;
            ls->bonds[best_bond].order++;
            ls->lone_pairs[term]--;
            recompute_formal_charges(mol, ls);
        }
    }

    recompute_formal_charges(mol, ls);

    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        int electrons = electrons_on_atom(ls, i);
        if (!shell_satisfied(mol, i, electrons, i == mol->central)) {
            *reason = INVALID_SHELL_RULE;
            return false;
        }
    }

    if (formal_charge_sum(mol, ls) != mol->charge) {
        *reason = INVALID_FORMAL_CHARGE_SUM;
        return false;
    }

    return true;
}

void generate_resonance(Molecule *mol)
{
    mol->num_res = 0;
    mol->cur_res = 0;
    mol->invalid_reason = INVALID_NONE;

    if (mol->num_atoms == 0) {
        mol->invalid_reason = INVALID_NO_ATOMS;
        return;
    }

    mol->total_ve = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        mol->total_ve += elements[mol->atoms[i].elem].valence;
    }
    mol->total_ve -= mol->charge;

    uint8_t candidates[MAX_ATOMS];
    uint8_t elem_counts[NUM_ELEMENTS];
    uint8_t n_candidates = gather_center_candidates(mol, candidates);
    memset(elem_counts, 0, sizeof(elem_counts));
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        elem_counts[mol->atoms[i].elem]++;
    }

    bool found_valid = false;
    uint8_t best_center = find_central(mol);
    LewisStructure best_ls;
    InvalidReason first_reason = INVALID_NONE;

    int best_sum_abs_fc = 0;
    int best_nonzero_fc = 0;
    int best_abs_central_fc = 0;
    uint8_t best_count = 0;
    bool best_terminal = false;
    uint8_t best_eneg = 0;
    uint8_t best_period = 0;
    uint8_t best_atomic_num = 0;

    for (uint8_t ci = 0; ci < n_candidates; ci++) {
        mol->central = candidates[ci];

        LewisStructure cand_ls;
        InvalidReason reason = INVALID_NONE;
        if (!generate_structure(mol, &cand_ls, &reason)) {
            if (first_reason == INVALID_NONE) {
                first_reason = reason;
            }
            continue;
        }

        int cand_sum_abs_fc = 0;
        int cand_nonzero_fc = 0;
        int cand_abs_central_fc = 0;
        score_structure(mol, &cand_ls, &cand_sum_abs_fc, &cand_nonzero_fc, &cand_abs_central_fc);

        const Element *cand_elem = &elements[mol->atoms[mol->central].elem];
        uint8_t cand_count = elem_counts[mol->atoms[mol->central].elem];
        bool cand_terminal = center_is_terminal_elem(mol->atoms[mol->central].elem);

        if (!found_valid ||
            candidate_is_better(cand_sum_abs_fc,
                                cand_nonzero_fc,
                                cand_abs_central_fc,
                                cand_count,
                                cand_terminal,
                                cand_elem->eneg,
                                cand_elem->period,
                                cand_elem->atomic_num,
                                best_sum_abs_fc,
                                best_nonzero_fc,
                                best_abs_central_fc,
                                best_count,
                                best_terminal,
                                best_eneg,
                                best_period,
                                best_atomic_num)) {
            found_valid = true;
            best_center = mol->central;
            memcpy(&best_ls, &cand_ls, sizeof(best_ls));

            best_sum_abs_fc = cand_sum_abs_fc;
            best_nonzero_fc = cand_nonzero_fc;
            best_abs_central_fc = cand_abs_central_fc;
            best_count = cand_count;
            best_terminal = cand_terminal;
            best_eneg = cand_elem->eneg;
            best_period = cand_elem->period;
            best_atomic_num = cand_elem->atomic_num;
        }
    }

    if (!found_valid) {
        mol->invalid_reason = (first_reason == INVALID_NONE) ? INVALID_SKELETON : first_reason;
        return;
    }

    mol->central = best_center;
    memcpy(&mol->res[0], &best_ls, sizeof(best_ls));
    mol->invalid_reason = INVALID_NONE;
    mol->num_res = 1;

    for (uint8_t seed_idx = 0; seed_idx < mol->num_res && mol->num_res < MAX_RESONANCE; seed_idx++) {
        LewisStructure seed;
        memcpy(&seed, &mol->res[seed_idx], sizeof(seed));

        for (uint8_t src = 0; src < seed.num_bonds && mol->num_res < MAX_RESONANCE; src++) {
            if (!(seed.bonds[src].a == mol->central || seed.bonds[src].b == mol->central)) {
                continue;
            }
            if (seed.bonds[src].order <= 1) {
                continue;
            }

            uint8_t src_term = (seed.bonds[src].a == mol->central) ? seed.bonds[src].b : seed.bonds[src].a;
            uint8_t src_elem = mol->atoms[src_term].elem;
            if (is_protonated_terminal_oxygen(mol, &seed, src_term)) continue;
            uint8_t shift = seed.bonds[src].order - 1;

            for (uint8_t dst = 0; dst < seed.num_bonds && mol->num_res < MAX_RESONANCE; dst++) {
                if (dst == src) continue;
                if (!(seed.bonds[dst].a == mol->central || seed.bonds[dst].b == mol->central)) {
                    continue;
                }

                uint8_t dst_term = (seed.bonds[dst].a == mol->central) ? seed.bonds[dst].b : seed.bonds[dst].a;
                if (mol->atoms[dst_term].elem != src_elem) continue;
                if (mol->atoms[dst_term].elem == ELEM_H) continue;
                if (is_protonated_terminal_oxygen(mol, &seed, dst_term)) continue;
                if (seed.bonds[dst].order >= seed.bonds[src].order) continue;
                if ((uint8_t)(seed.bonds[dst].order + shift) > 3) continue;
                if (seed.lone_pairs[dst_term] < shift) continue;

                LewisStructure cand;
                memcpy(&cand, &seed, sizeof(cand));

                cand.bonds[src].order = 1;
                cand.lone_pairs[src_term] += shift;

                cand.bonds[dst].order += shift;
                cand.lone_pairs[dst_term] -= shift;

                recompute_formal_charges(mol, &cand);
                if (formal_charge_sum(mol, &cand) != mol->charge) continue;

                bool valid = true;
                for (uint8_t i = 0; i < mol->num_atoms; i++) {
                    int electrons = electrons_on_atom(&cand, i);
                    if (!shell_satisfied(mol, i, electrons, i == mol->central)) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;
                if (resonance_exists(mol, &cand)) continue;

                memcpy(&mol->res[mol->num_res], &cand, sizeof(cand));
                mol->num_res++;
            }
        }
    }
}

static void fill_vsepr_fallback(VseprInfo *out)
{
    switch (out->valence_pairs) {
        case 1:
            out->ep_geometry = "Linear";
            out->hybridization = "s";
            out->bond_angle = "180";
            break;
        case 2:
            out->ep_geometry = "Linear";
            out->hybridization = "sp";
            out->bond_angle = "180";
            break;
        case 3:
            out->ep_geometry = "Trigonal Planar";
            out->hybridization = "sp2";
            out->bond_angle = "120";
            break;
        case 4:
            out->ep_geometry = "Tetrahedral";
            out->hybridization = "sp3";
            out->bond_angle = "109.5";
            break;
        case 5:
            out->ep_geometry = "Trigonal Bipyramidal";
            out->hybridization = "sp3d";
            out->bond_angle = "90, 120";
            break;
        case 6:
            out->ep_geometry = "Octahedral";
            out->hybridization = "sp3d2";
            out->bond_angle = "90";
            break;
        case 7:
            out->ep_geometry = "Pentagonal Bipyramidal";
            out->hybridization = "sp3d3";
            out->bond_angle = "72, 90";
            break;
        default:
            out->ep_geometry = "Unknown";
            out->hybridization = "Unknown";
            out->bond_angle = "N/A";
            break;
    }

    if (out->bond_pairs == 0) {
        out->shape = "No Bonded Atoms";
    } else if (out->bond_pairs <= 1) {
        out->shape = "Linear";
    } else {
        out->shape = out->ep_geometry;
    }
}

bool lewis_get_vsepr_info(const Molecule *mol, const LewisStructure *ls, VseprInfo *out)
{
    if (out == NULL) {
        return false;
    }

    memset(out, 0, sizeof(*out));

    if (mol == NULL || ls == NULL) {
        return false;
    }
    if (mol->num_atoms == 0 || mol->central >= mol->num_atoms) {
        return false;
    }

    uint8_t bond_pairs = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if (ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central) {
            bond_pairs++;
        }
    }

    out->bond_pairs = bond_pairs;
    out->lone_pairs = ls->lone_pairs[mol->central];
    out->valence_pairs = (uint8_t)(out->bond_pairs + out->lone_pairs);

    for (uint8_t i = 0; i < (uint8_t)(sizeof(vsepr_rows) / sizeof(vsepr_rows[0])); i++) {
        if (vsepr_rows[i].valence_pairs == out->valence_pairs &&
            vsepr_rows[i].bond_pairs == out->bond_pairs &&
            vsepr_rows[i].lone_pairs == out->lone_pairs) {
            out->ep_geometry = vsepr_rows[i].ep_geometry;
            out->shape = vsepr_rows[i].shape;
            out->hybridization = vsepr_rows[i].hybridization;
            out->bond_angle = vsepr_rows[i].bond_angle;
            return true;
        }
    }

    fill_vsepr_fallback(out);
    return true;
}

const char *invalid_reason_message(InvalidReason reason)
{
    switch (reason) {
        case INVALID_NONE:
            return "No error";
        case INVALID_NO_ATOMS:
            return "No atoms selected";
        case INVALID_NEGATIVE_ELECTRONS:
            return "Invalid charge for selected atoms";
        case INVALID_ODD_ELECTRONS:
            return "Odd electron count (radicals unsupported)";
        case INVALID_SKELETON:
            return "Cannot build a valid bond skeleton";
        case INVALID_LEFTOVER_ELECTRONS:
            return "Could not place all valence electrons";
        case INVALID_SHELL_RULE:
            return "Octet/duet shell constraints failed";
        case INVALID_FORMAL_CHARGE_SUM:
            return "Formal charge sum does not match ion charge";
        default:
            return "No valid structure for this composition";
    }
}
