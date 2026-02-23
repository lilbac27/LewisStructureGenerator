#include "lewis_model.h"

#include <string.h>

const Element elements[NUM_ELEMENTS] = {
    /* sym   name          Z   val bc  eneg per grp color */
    { "H",  "Hydrogen",     1,  1, 1,  22,  1,  1, COL_BLACK },    /* 0  */
    { "He", "Helium",       2,  2, 0,   0,  1, 18, COL_BLACK },    /* 1  */
    { "Li", "Lithium",      3,  1, 1,  10,  2,  1, COL_BLACK },    /* 2  */
    { "Be", "Beryllium",    4,  2, 2,  16,  2,  2, COL_BLACK },    /* 3  */
    { "B",  "Boron",        5,  3, 3,  20,  2, 13, COL_BLACK },    /* 4  */
    { "C",  "Carbon",       6,  4, 4,  26,  2, 14, COL_BLACK },    /* 5  */
    { "N",  "Nitrogen",     7,  5, 3,  30,  2, 15, COL_BLACK },    /* 6  */
    { "O",  "Oxygen",       8,  6, 2,  34,  2, 16, COL_BLACK },    /* 7  */
    { "F",  "Fluorine",     9,  7, 1,  40,  2, 17, COL_BLACK },    /* 8  */
    { "Ne", "Neon",        10,  8, 0,   0,  2, 18, COL_BLACK },    /* 9  */
    { "Na", "Sodium",      11,  1, 1,   9,  3,  1, COL_BLACK },    /* 10 */
    { "Mg", "Magnesium",   12,  2, 2,  13,  3,  2, COL_BLACK },    /* 11 */
    { "Al", "Aluminum",    13,  3, 3,  16,  3, 13, COL_BLACK },    /* 12 */
    { "Si", "Silicon",     14,  4, 4,  19,  3, 14, COL_BLACK },    /* 13 */
    { "P",  "Phosphorus",  15,  5, 5,  22,  3, 15, COL_BLACK },    /* 14 */
    { "S",  "Sulfur",      16,  6, 6,  26,  3, 16, COL_BLACK },    /* 15 */
    { "Cl", "Chlorine",    17,  7, 1,  32,  3, 17, COL_BLACK },    /* 16 */
    { "Ar", "Argon",       18,  8, 0,   0,  3, 18, COL_BLACK },    /* 17 */
    { "K",  "Potassium",   19,  1, 1,   8,  4,  1, COL_BLACK },    /* 18 */
    { "Ca", "Calcium",     20,  2, 2,  10,  4,  2, COL_BLACK },    /* 19 */
    { "Ga", "Gallium",     31,  3, 3,  18,  4, 13, COL_BLACK },    /* 20 */
    { "Ge", "Germanium",   32,  4, 4,  20,  4, 14, COL_BLACK },    /* 21 */
    { "As", "Arsenic",     33,  5, 5,  22,  4, 15, COL_BLACK },    /* 22 */
    { "Se", "Selenium",    34,  6, 6,  26,  4, 16, COL_BLACK },    /* 23 */
    { "Br", "Bromine",     35,  7, 1,  30,  4, 17, COL_BLACK },    /* 24 */
    { "Kr", "Krypton",     36,  8, 2,  30,  4, 18, COL_BLACK },    /* 25 */
    { "Rb", "Rubidium",    37,  1, 1,   8,  5,  1, COL_BLACK },    /* 26 */
    { "Sr", "Strontium",   38,  2, 2,  10,  5,  2, COL_BLACK },    /* 27 */
    { "In", "Indium",      49,  3, 3,  18,  5, 13, COL_BLACK },    /* 28 */
    { "Sn", "Tin",         50,  4, 4,  20,  5, 14, COL_BLACK },    /* 29 */
    { "Sb", "Antimony",    51,  5, 5,  21,  5, 15, COL_BLACK },    /* 30 */
    { "Te", "Tellurium",   52,  6, 6,  21,  5, 16, COL_BLACK },    /* 31 */
    { "I",  "Iodine",      53,  7, 1,  27,  5, 17, COL_BLACK },    /* 32 */
    { "Xe", "Xenon",       54,  8, 4,  26,  5, 18, COL_BLACK },    /* 33 */
};

/* Map (period, group) -> element index. ELEM_NONE = empty cell. */
uint8_t pt_grid[PT_ROWS][PT_COLS];

void init_pt_grid(void)
{
    memset(pt_grid, ELEM_NONE, sizeof(pt_grid));
    for (uint8_t i = 0; i < NUM_ELEMENTS; i++) {
        uint8_t r = elements[i].period - 1;
        uint8_t c = elements[i].group - 1;
        if (r < PT_ROWS && c < PT_COLS) {
            pt_grid[r][c] = i;
        }
    }
}

void molecule_reset(Molecule *mol)
{
    memset(mol, 0, sizeof(*mol));
    mol->invalid_reason = INVALID_NONE;
}
