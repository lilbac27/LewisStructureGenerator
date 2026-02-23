#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "../src/lewis_engine.h"
#include "../src/lewis_model.h"

#define ELEM_B_IDX   4
#define ELEM_P_IDX   14
#define ELEM_F_IDX   8
#define ELEM_CL_IDX  16
#define ELEM_SE_IDX  23
#define ELEM_I_IDX   32
#define ELEM_XE_IDX  33

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

static int central_bond_count_by_order(const Molecule *mol, const LewisStructure *ls, uint8_t order)
{
    int count = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        bool touches_central = (ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central);
        if (touches_central && ls->bonds[b].order == order) {
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

static bool vsepr_info_matches(const VseprInfo *info,
                               uint8_t valence_pairs,
                               uint8_t bond_pairs,
                               uint8_t lone_pairs,
                               const char *ep_geometry,
                               const char *shape,
                               const char *hybridization)
{
    if (info->valence_pairs != valence_pairs) return false;
    if (info->bond_pairs != bond_pairs) return false;
    if (info->lone_pairs != lone_pairs) return false;
    if (info->ep_geometry == NULL || strcmp(info->ep_geometry, ep_geometry) != 0) return false;
    if (info->shape == NULL || strcmp(info->shape, shape) != 0) return false;
    if (info->hybridization == NULL || strcmp(info->hybridization, hybridization) != 0) return false;
    return true;
}

static bool test_vsepr_co2(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_C, ELEM_O, ELEM_O };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    return vsepr_info_matches(&info, 2, 2, 0, "Linear", "Linear", "sp");
}

static bool test_vsepr_nitrate(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_N, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    return vsepr_info_matches(&info, 3, 3, 0, "Trigonal Planar", "Trigonal Planar", "sp2");
}

static bool test_vsepr_ammonium(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_N, ELEM_H, ELEM_H, ELEM_H, ELEM_H };
    build_and_generate(&mol, 1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    return vsepr_info_matches(&info, 4, 4, 0, "Tetrahedral", "Tetrahedral", "sp3");
}

static bool test_vsepr_water(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_O, ELEM_H, ELEM_H };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    return vsepr_info_matches(&info, 4, 2, 2, "Tetrahedral", "Bent", "sp3");
}

static bool test_vsepr_pcl5(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_P_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    return vsepr_info_matches(&info, 5, 5, 0, "Trigonal Bipyramidal", "Trigonal Bipyramidal", "sp3d");
}

static bool test_vsepr_sf6(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_S, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    return vsepr_info_matches(&info, 6, 6, 0, "Octahedral", "Octahedral", "sp3d2");
}

static bool test_vsepr_h2_no_null(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_H, ELEM_H };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    if (!vsepr_info_matches(&info, 1, 1, 0, "Linear", "Linear", "s")) return false;
    return true;
}

static bool test_vsepr_invalid_guard(void)
{
    Molecule mol;
    LewisStructure ls;
    VseprInfo info;

    memset(&ls, 0, sizeof(ls));
    if (lewis_get_vsepr_info(NULL, NULL, &info)) return false;

    molecule_reset(&mol);
    if (lewis_get_vsepr_info(&mol, &ls, &info)) return false;
    if (lewis_get_vsepr_info(&mol, &ls, NULL)) return false;

    mol.num_atoms = 1;
    mol.central = 2;
    mol.atoms[0].elem = ELEM_C;
    if (lewis_get_vsepr_info(&mol, &ls, &info)) return false;

    return true;
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

static bool test_h2_basic(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_H, ELEM_H };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;

    const LewisStructure *ls = &mol.res[0];
    if (ls->num_bonds != 1) return false;
    if (ls->bonds[0].order != 1) return false;
    if (ls->lone_pairs[0] != 0 || ls->lone_pairs[1] != 0) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 1, 1, 0, "Linear", "Linear", "s")) return false;
    return true;
}

static bool test_o2_basic(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_O, ELEM_O };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;

    const LewisStructure *ls = &mol.res[0];
    if (ls->num_bonds != 1) return false;
    if (ls->bonds[0].order != 2) return false;
    if (ls->lone_pairs[mol.central] != 2) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 3, 1, 2, "Trigonal Planar", "Linear", "sp2")) return false;
    return true;
}

static bool test_n2_basic(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_N, ELEM_N };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;

    const LewisStructure *ls = &mol.res[0];
    if (ls->num_bonds != 1) return false;
    if (ls->bonds[0].order != 3) return false;
    if (ls->lone_pairs[mol.central] != 1) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 2, 1, 1, "Linear", "Linear", "sp")) return false;
    return true;
}

static bool test_ch4_basic(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_C, ELEM_H, ELEM_H, ELEM_H, ELEM_H };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_C) return false;

    const LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 4) return false;
    if (central_bond_count_by_order(&mol, ls, 1) != 4) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 4, 4, 0, "Tetrahedral", "Tetrahedral", "sp3")) return false;
    return true;
}

static bool test_nh3_basic(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_N, ELEM_H, ELEM_H, ELEM_H };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_N) return false;

    const LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 3) return false;
    if (central_bond_count_by_order(&mol, ls, 1) != 3) return false;
    if (ls->lone_pairs[mol.central] != 1) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 4, 3, 1, "Tetrahedral", "Trigonal Pyramidal", "sp3")) return false;
    return true;
}

static bool test_no2_minus_resonance(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_N, ELEM_O, ELEM_O };
    build_and_generate(&mol, -1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 2) return false;
    if (mol.atoms[mol.central].elem != ELEM_N) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 3) return false;
        if (central_double_bond_count(&mol, ls) != 1) return false;
        if (central_bond_count_by_order(&mol, ls, 1) != 1) return false;
        if (ls->lone_pairs[mol.central] != 1) return false;
    }

    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    if (!vsepr_info_matches(&info, 3, 2, 1, "Trigonal Planar", "Bent", "sp2")) return false;
    return true;
}

static bool test_clo3_minus_oxyanion(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_CL_IDX, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 3) return false;
    if (mol.atoms[mol.central].elem != ELEM_CL_IDX) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 5) return false;
        if (central_double_bond_count(&mol, ls) != 2) return false;
        if (central_bond_count_by_order(&mol, ls, 1) != 1) return false;
        if (ls->lone_pairs[mol.central] != 1) return false;
        if (formal_charge_sum(&mol, ls) != -1) return false;
    }

    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    if (!vsepr_info_matches(&info, 4, 3, 1, "Tetrahedral", "Trigonal Pyramidal", "sp3")) return false;
    return true;
}

static bool test_clo4_minus_oxyanion(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_CL_IDX, ELEM_O, ELEM_O, ELEM_O, ELEM_O };
    build_and_generate(&mol, -1, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 4) return false;
    if (mol.atoms[mol.central].elem != ELEM_CL_IDX) return false;

    for (uint8_t i = 0; i < mol.num_res; i++) {
        const LewisStructure *ls = &mol.res[i];
        if (bond_order_sum(ls, mol.central) != 7) return false;
        if (central_double_bond_count(&mol, ls) != 3) return false;
        if (central_bond_count_by_order(&mol, ls, 1) != 1) return false;
        if (ls->lone_pairs[mol.central] != 0) return false;
        if (formal_charge_sum(&mol, ls) != -1) return false;
    }

    if (!lewis_get_vsepr_info(&mol, &mol.res[mol.cur_res], &info)) return false;
    if (!vsepr_info_matches(&info, 4, 4, 0, "Tetrahedral", "Tetrahedral", "sp3")) return false;
    return true;
}

static bool test_xef2_hypervalent(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_XE_IDX, ELEM_F_IDX, ELEM_F_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_XE_IDX) return false;

    const LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 2) return false;
    if (central_bond_count_by_order(&mol, ls, 1) != 2) return false;
    if (ls->lone_pairs[mol.central] != 3) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 5, 2, 3, "Trigonal Bipyramidal", "Linear", "sp3d")) return false;
    return true;
}

static bool test_xef4_hypervalent(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_XE_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_XE_IDX) return false;

    const LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 4) return false;
    if (central_bond_count_by_order(&mol, ls, 1) != 4) return false;
    if (ls->lone_pairs[mol.central] != 2) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 6, 4, 2, "Octahedral", "Square Planar", "sp3d2")) return false;
    return true;
}

static bool test_if7_vsepr_supported(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_I_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX, ELEM_F_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_I_IDX) return false;

    const LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 7) return false;
    if (central_bond_count_by_order(&mol, ls, 1) != 7) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 7, 7, 0, "Pentagonal Bipyramidal", "Pentagonal Bipyramidal", "sp3d3")) return false;
    return true;
}

static bool test_cse2_center_and_geometry(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_C, ELEM_SE_IDX, ELEM_SE_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_C) return false;

    LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 4) return false;
    if (central_double_bond_count(&mol, ls) != 2) return false;
    if (ls->lone_pairs[mol.central] != 0) return false;
    if (formal_charge_sum(&mol, ls) != 0) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 2, 2, 0, "Linear", "Linear", "sp")) return false;
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

static bool test_icl5_expanded_valence(void)
{
    Molecule mol;
    VseprInfo info;
    const uint8_t atoms[] = { ELEM_I_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX, ELEM_CL_IDX };
    build_and_generate(&mol, 0, atoms, (uint8_t)(sizeof(atoms) / sizeof(atoms[0])));

    if (!success_invariants(&mol)) return false;
    if (mol.num_res != 1) return false;
    if (mol.atoms[mol.central].elem != ELEM_I_IDX) return false;

    LewisStructure *ls = &mol.res[0];
    if (bond_order_sum(ls, mol.central) != 5) return false;
    if (ls->lone_pairs[mol.central] != 1) return false;

    int central_bonds = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        bool touches_central = (ls->bonds[b].a == mol.central || ls->bonds[b].b == mol.central);
        if (!touches_central) continue;
        central_bonds++;
        if (ls->bonds[b].order != 1) return false;
    }
    if (central_bonds != 5) return false;

    if (!lewis_get_vsepr_info(&mol, ls, &info)) return false;
    if (!vsepr_info_matches(&info, 6, 5, 1, "Octahedral", "Square Pyramidal", "sp3d2")) return false;
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
        { "H2 basic", test_h2_basic },
        { "O2 basic", test_o2_basic },
        { "N2 basic", test_n2_basic },
        { "CH4 basic", test_ch4_basic },
        { "NH3 basic", test_nh3_basic },
        { "CO2", test_co2 },
        { "CSe2 center/geometry", test_cse2_center_and_geometry },
        { "NO2- resonance", test_no2_minus_resonance },
        { "ClO3- oxyanion", test_clo3_minus_oxyanion },
        { "ClO4- oxyanion", test_clo4_minus_oxyanion },
        { "NO3-", test_nitrate },
        { "CO3^2-", test_carbonate },
        { "SO4^2-", test_sulfate },
        { "PO4^3-", test_phosphate },
        { "NH4+", test_ammonium },
        { "BF3 (incomplete octet)", test_bf3_incomplete_octet },
        { "SF6 (expanded valence)", test_sf6_expanded_valence },
        { "PCl5 (expanded valence)", test_pcl5_expanded_valence },
        { "ICl5 (expanded valence)", test_icl5_expanded_valence },
        { "XeF2 hypervalent", test_xef2_hypervalent },
        { "XeF4 hypervalent", test_xef4_hypervalent },
        { "IF7 VSEPR supported", test_if7_vsepr_supported },
        { "VSEPR CO2", test_vsepr_co2 },
        { "VSEPR NO3-", test_vsepr_nitrate },
        { "VSEPR NH4+", test_vsepr_ammonium },
        { "VSEPR H2O", test_vsepr_water },
        { "VSEPR PCl5", test_vsepr_pcl5 },
        { "VSEPR SF6", test_vsepr_sf6 },
        { "VSEPR H2 no-null", test_vsepr_h2_no_null },
        { "VSEPR invalid-guard", test_vsepr_invalid_guard },
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
