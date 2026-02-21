#include <stdbool.h>
#include <stdio.h>

#include "../src/lewis_engine.h"
#include "../src/lewis_model.h"

#define ELEM_B_IDX   4
#define ELEM_P_IDX   14
#define ELEM_F_IDX   8
#define ELEM_CL_IDX  16

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

static int formal_charge_sum(const Molecule *mol, const LewisStructure *ls)
{
    int sum = 0;
    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        sum += ls->formal_charge[i];
    }
    return sum;
}

static int central_double_bond_count(const Molecule *mol, const LewisStructure *ls)
{
    int count = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        bool touches_central = (ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central);
        if (touches_central && ls->bonds[b].order == 2) {
            count++;
        }
    }
    return count;
}

static bool structures_equal(const Molecule *mol, const LewisStructure *a, const LewisStructure *b)
{
    if (a->num_bonds != b->num_bonds) return false;
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

static bool all_resonance_unique(const Molecule *mol)
{
    for (uint8_t i = 0; i < mol->num_res; i++) {
        for (uint8_t j = i + 1; j < mol->num_res; j++) {
            if (structures_equal(mol, &mol->res[i], &mol->res[j])) {
                return false;
            }
        }
    }
    return true;
}

static bool all_formal_charge_sums_match(const Molecule *mol)
{
    for (uint8_t i = 0; i < mol->num_res; i++) {
        if (formal_charge_sum(mol, &mol->res[i]) != mol->charge) {
            return false;
        }
    }
    return true;
}

static bool success_invariants(const Molecule *mol)
{
    if (mol->invalid_reason != INVALID_NONE) return false;
    if (mol->num_res == 0) return false;
    if (mol->central >= mol->num_atoms) return false;
    if (!all_formal_charge_sums_match(mol)) return false;
    if (!all_resonance_unique(mol)) return false;
    return true;
}

static void build_molecule(Molecule *mol, int8_t charge, const uint8_t atoms[], uint8_t count)
{
    molecule_reset(mol);
    mol->charge = charge;
    for (uint8_t i = 0; i < count; i++) {
        mol->atoms[i].elem = atoms[i];
    }
    mol->num_atoms = count;
}

static void build_and_generate(Molecule *mol, int8_t charge, const uint8_t atoms[], uint8_t count)
{
    build_molecule(mol, charge, atoms, count);
    generate_resonance(mol);
}

static bool test_co2(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_C, ELEM_O, ELEM_O };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_C) return false;

    LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 4) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;
    return true;
}

static bool test_nitrate(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_N, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 3) return false;
    if (mol.atoms[mol.central].elem != ELEM_N) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 4) return false;
        if (central_double_bond_count(&mol, ls) != 1) return false;
    }
    return true;
}

static bool test_sulfate(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_S, ELEM_O, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -2, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != MAX_RESONANCE) return false;
    if (mol.atoms[mol.central].elem != ELEM_S) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 6) return false;
        if (central_double_bond_count(&mol, ls) != 2) return false;
    }
    return true;
}

static bool test_ammonium(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_N, ELEM_H, ELEM_H, ELEM_H, ELEM_H };
    build_and_generate(&mol, 1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_N) return false;

    LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 4) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;
    if (ls->formal_charge[mol.central] != 1) return false;
    return true;
}

static bool test_carbonate(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_C, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -2, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 3) return false;
    if (mol.atoms[mol.central].elem != ELEM_C) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 4) return false;
        if (central_double_bond_count(&mol, ls) != 1) return false;
    }
    return true;
}

static bool test_phosphate(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_P_IDX, ELEM_O, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -3, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 4) return false;
    if (mol.atoms[mol.central].elem != ELEM_P_IDX) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 5) return false;
        if (central_double_bond_count(&mol, ls) != 1) return false;
    }
    return true;
}

static bool test_bf3_incomplete_octet(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_B_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_B_IDX) return false;

    LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 3) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;
    return true;
}

static bool test_sf6_expanded_valence(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_S, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_S) return false;
    if (bond_order_sum(&mol.res[0], mol.central) != 6) return false;
    return true;
}

static bool test_pcl5_expanded_valence(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_P_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_P_IDX) return false;
    if (bond_order_sum(&mol.res[0], mol.central) != 5) return false;
    return true;
}

static bool test_no_atoms_failure(void)
{
    Molecule mol;
    molecule_reset(&mol);
    generate_resonance(&mol);

    if (mol.num_res != 0) return false;
    if (mol.invalid_reason != INVALID_NO_ATOMS) return false;
    return true;
}

static bool test_negative_electrons_failure(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_H };
    build_and_generate(&mol, 2, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (mol.num_res != 0) return false;
    if (mol.total_ve >= 0) return false;
    if (mol.invalid_reason != INVALID_NEGATIVE_ELECTRONS) return false;
    return true;
}

static bool test_skeleton_failure(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_HE, ELEM_HE };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (mol.num_res != 0) return false;
    if (mol.invalid_reason != INVALID_SKELETON) return false;
    return true;
}

static bool test_odd_electron_failure(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_N, ELEM_O };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (mol.num_res != 0) return false;
    if (mol.invalid_reason != INVALID_ODD_ELECTRONS) return false;
    return true;
}

int main(void)
{
    struct {
        const char *name;
        bool (*fn)(void);
    } tests[] = {
        { "CO2", test_co2 },
        { "NO3-", test_nitrate },
        { "CO3^2-", test_carbonate },
        { "SO4^2-", test_sulfate },
        { "PO4^3-", test_phosphate },
        { "NH4+", test_ammonium },
        { "BF3 (incomplete octet)", test_bf3_incomplete_octet },
        { "SF6 (expanded valence)", test_sf6_expanded_valence },
        { "PCl5 (expanded valence)", test_pcl5_expanded_valence },
        { "No-atoms failure", test_no_atoms_failure },
        { "Negative-electron failure", test_negative_electrons_failure },
        { "Skeleton failure", test_skeleton_failure },
        { "Odd-electron rejection", test_odd_electron_failure },
    };

    int passed = 0;
    int total = (int)(sizeof(tests) / sizeof(tests[0]));

    for (int i = 0; i < total; i++) {
        bool ok = tests[i].fn();
        if (ok) {
            printf("[PASS] %s\n", tests[i].name);
            passed++;
        } else {
            printf("[FAIL] %s\n", tests[i].name);
        }
    }

    printf("\n%d/%d tests passed\n", passed, total);
    return (passed == total) ? 0 : 1;
}
