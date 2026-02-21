#include "lewis_engine.h"

#include <string.h>

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
        if (e->group == 17 && limit < 4) limit = 4;
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

    mol->central = find_central(mol);
    mol->total_ve = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        mol->total_ve += elements[mol->atoms[i].elem].valence;
    }
    mol->total_ve -= mol->charge;

    if (!generate_structure(mol, &mol->res[0], &mol->invalid_reason)) {
        return;
    }

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
