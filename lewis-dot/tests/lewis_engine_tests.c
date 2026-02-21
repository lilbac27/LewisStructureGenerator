#include <stdbool.h>
#include <stdio.h>

#include "../src/lewis_engine.h"
#include "../src/lewis_model.h"

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

static void build_molecule(Molecule *mol, int8_t charge, const uint8_t atoms[], uint8_t count)
{
    molecule_reset(mol);
    mol->charge = charge;
    for (uint8_t i = 0; i < count; i++) {
        mol->atoms[i].elem = atoms[i];
    }
    mol->num_atoms = count;
}

static bool test_co2(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_C, ELEM_O, ELEM_O };
    build_molecule(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));
    generate_resonance(&mol);

    if (mol.num_res < 1) return false;
    if (mol.invalid_reason != INVALID_NONE) return false;
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
    build_molecule(&mol, -1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));
    generate_resonance(&mol);

    if (mol.num_res < 3) return false;
    if (mol.invalid_reason != INVALID_NONE) return false;
    if (mol.atoms[mol.central].elem != ELEM_N) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (formal_charge_sum(&mol, ls) != -1) return false;
        if (bond_order_sum(ls, mol.central) != 4) return false;
        if (central_double_bond_count(&mol, ls) != 1) return false;
    }
    return true;
}

static bool test_sulfate(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_S, ELEM_O, ELEM_O, ELEM_O, ELEM_O };
    build_molecule(&mol, -2, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));
    generate_resonance(&mol);

    if (mol.num_res < 2) return false;
    if (mol.invalid_reason != INVALID_NONE) return false;
    if (mol.atoms[mol.central].elem != ELEM_S) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (formal_charge_sum(&mol, ls) != -2) return false;
    }
    return true;
}

static bool test_ammonium(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_N, ELEM_H, ELEM_H, ELEM_H, ELEM_H };
    build_molecule(&mol, 1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));
    generate_resonance(&mol);

    if (mol.num_res != 1) return false;
    if (mol.invalid_reason != INVALID_NONE) return false;
    if (mol.atoms[mol.central].elem != ELEM_N) return false;

    LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 4) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;
    if (formal_charge_sum(&mol, ls) != 1) return false;
    return true;
}

static bool test_odd_electron_failure(void)
{
    Molecule mol;
    const uint8_t atoms[] = { ELEM_N, ELEM_O };
    build_molecule(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));
    generate_resonance(&mol);

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
        { "SO4^2-", test_sulfate },
        { "NH4+", test_ammonium },
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
