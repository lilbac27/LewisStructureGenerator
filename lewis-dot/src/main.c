/*
 * Lewis Dot Structure Generator for TI-84 Plus CE
 * Generates and displays Lewis dot structures for small molecules.
 * Supports resonance structures with [left]/[right] navigation.
 *
 * Controls (Periodic Table screen):
 *   Arrow keys  — move cursor on periodic table
 *   [enter]     — add highlighted element to molecule
 *   [del]       — remove last added atom
 *   [alpha]     — cycle charge (0, +1, +2, -1, -2)
 *   [2nd]       — generate Lewis structure
 *   [mode]      — quit
 *
 * Controls (Lewis Structure screen):
 *   [left]/[right] — cycle resonance structures
 *   [clear]        — return to periodic table
 */

#include <graphx.h>
#include <keypadc.h>
#include <sys/timers.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

/* ── Screen layout constants ── */
#define SCR_W           320
#define SCR_H           240

/* Top info bar */
#define INFO_Y          0
#define INFO_H          30

/* Selected atoms bar */
#define SEL_Y           32
#define SEL_H           43

/* Periodic table area */
#define PT_Y            80
#define PT_CELL_W       17
#define PT_CELL_H       16
#define PT_COLS         18
#define PT_ROWS         5
#define PT_LEFT_PAD     2

/* Lewis structure screen */
#define LEWIS_CENTER_X  160
#define LEWIS_CENTER_Y  135
#define BOND_LEN        55
#define DOT_R           2
#define DOT_DIST        14

/* Molecule limits */
#define MAX_ATOMS       12
#define MAX_HEAVY       6
#define MAX_BONDS       12
#define MAX_RESONANCE   6

/* Frame rate target */
#define TARGET_FPS      20
#define FRAME_TICKS     (32768 / TARGET_FPS)

/* ── Color palette indices (default TI palette) ── */
#define COL_BLACK       0x00
#define COL_WHITE       0xFF
#define COL_RED         0xE0
#define COL_BLUE        0x10
#define COL_GREEN       0x04
#define COL_YELLOW      0xE7
#define COL_CYAN        0x1F
#define COL_MAGENTA     0xE3
#define COL_ORANGE      0xE4
#define COL_GRAY        0xB5
#define COL_DKGRAY      0x6B

/* ── Element info card constants ── */
#define CARD_X          56
#define CARD_Y          84
#define CARD_W          86
#define CARD_H          80

/* ── Element data ── */
typedef struct {
    char     symbol[3];
    char     name[14];
    uint8_t  atomic_num;   /* atomic number (Z) */
    uint8_t  valence;      /* valence electrons */
    uint8_t  bond_cap;     /* typical bonding capacity */
    uint8_t  eneg;         /* electronegativity * 10 (integer) */
    uint8_t  period;       /* 1-based */
    uint8_t  group;        /* 1-based, 1-18 */
    uint8_t  color;        /* display color */
} Element;

/* Index constants for quick lookup */
#define ELEM_NONE  0xFF
#define NUM_ELEMENTS 34

static const Element elements[NUM_ELEMENTS] = {
    /* sym   name          Z   val bc  eneg per grp color */
    { "H",  "Hydrogen",     1,  1, 1,  22,  1,  1, COL_CYAN   },    /* 0  */
    { "He", "Helium",       2,  2, 0,   0,  1, 18, COL_MAGENTA},    /* 1  */
    { "Li", "Lithium",      3,  1, 1,  10,  2,  1, COL_ORANGE },    /* 2  */
    { "Be", "Beryllium",    4,  2, 2,  16,  2,  2, COL_ORANGE },    /* 3  */
    { "B",  "Boron",        5,  3, 3,  20,  2, 13, COL_YELLOW },    /* 4  */
    { "C",  "Carbon",       6,  4, 4,  26,  2, 14, COL_DKGRAY },    /* 5  */
    { "N",  "Nitrogen",     7,  5, 3,  30,  2, 15, COL_BLUE   },    /* 6  */
    { "O",  "Oxygen",       8,  6, 2,  34,  2, 16, COL_RED    },    /* 7  */
    { "F",  "Fluorine",     9,  7, 1,  40,  2, 17, COL_GREEN  },    /* 8  */
    { "Ne", "Neon",        10,  8, 0,   0,  2, 18, COL_MAGENTA},    /* 9  */
    { "Na", "Sodium",      11,  1, 1,   9,  3,  1, COL_ORANGE },    /* 10 */
    { "Mg", "Magnesium",   12,  2, 2,  13,  3,  2, COL_ORANGE },    /* 11 */
    { "Al", "Aluminum",    13,  3, 3,  16,  3, 13, COL_YELLOW },    /* 12 */
    { "Si", "Silicon",     14,  4, 4,  19,  3, 14, COL_YELLOW },    /* 13 */
    { "P",  "Phosphorus",  15,  5, 5,  22,  3, 15, COL_BLUE   },    /* 14 */
    { "S",  "Sulfur",      16,  6, 6,  26,  3, 16, COL_YELLOW },    /* 15 */
    { "Cl", "Chlorine",    17,  7, 1,  32,  3, 17, COL_GREEN  },    /* 16 */
    { "Ar", "Argon",       18,  8, 0,   0,  3, 18, COL_MAGENTA},    /* 17 */
    { "K",  "Potassium",   19,  1, 1,   8,  4,  1, COL_ORANGE },    /* 18 */
    { "Ca", "Calcium",     20,  2, 2,  10,  4,  2, COL_ORANGE },    /* 19 */
    { "Ga", "Gallium",     31,  3, 3,  18,  4, 13, COL_YELLOW },    /* 20 */
    { "Ge", "Germanium",   32,  4, 4,  20,  4, 14, COL_YELLOW },    /* 21 */
    { "As", "Arsenic",     33,  5, 5,  22,  4, 15, COL_BLUE   },    /* 22 */
    { "Se", "Selenium",    34,  6, 6,  26,  4, 16, COL_YELLOW },    /* 23 */
    { "Br", "Bromine",     35,  7, 1,  30,  4, 17, COL_GREEN  },    /* 24 */
    { "Kr", "Krypton",     36,  8, 2,  30,  4, 18, COL_MAGENTA},    /* 25 */
    { "Rb", "Rubidium",    37,  1, 1,   8,  5,  1, COL_ORANGE },    /* 26 */
    { "Sr", "Strontium",   38,  2, 2,  10,  5,  2, COL_ORANGE },    /* 27 */
    { "In", "Indium",      49,  3, 3,  18,  5, 13, COL_YELLOW },    /* 28 */
    { "Sn", "Tin",         50,  4, 4,  20,  5, 14, COL_YELLOW },    /* 29 */
    { "Sb", "Antimony",    51,  5, 5,  21,  5, 15, COL_BLUE   },    /* 30 */
    { "Te", "Tellurium",   52,  6, 6,  21,  5, 16, COL_YELLOW },    /* 31 */
    { "I",  "Iodine",      53,  7, 1,  27,  5, 17, COL_GREEN  },    /* 32 */
    { "Xe", "Xenon",       54,  8, 4,  26,  5, 18, COL_MAGENTA},    /* 33 */
};

/* Map (period, group) -> element index. ELEM_NONE = empty cell. */
static uint8_t pt_grid[PT_ROWS][PT_COLS];

static void init_pt_grid(void)
{
    memset(pt_grid, ELEM_NONE, sizeof(pt_grid));
    for (uint8_t i = 0; i < NUM_ELEMENTS; i++) {
        uint8_t r = elements[i].period - 1;
        uint8_t c = elements[i].group  - 1;
        if (r < PT_ROWS && c < PT_COLS)
            pt_grid[r][c] = i;
    }
}

/* ── Molecule representation ── */
typedef struct {
    uint8_t  elem;         /* index into elements[] */
} Atom;

typedef struct {
    uint8_t  a, b;         /* atom indices */
    uint8_t  order;        /* 1=single, 2=double, 3=triple */
} Bond;

typedef struct {
    uint8_t  lone_pairs[MAX_ATOMS]; /* lone pairs per atom */
    Bond     bonds[MAX_BONDS];
    uint8_t  num_bonds;
    int8_t   formal_charge[MAX_ATOMS];
} LewisStructure;

typedef struct {
    Atom     atoms[MAX_ATOMS];
    uint8_t  num_atoms;
    int8_t   charge;       /* overall molecular charge */

    /* Generated structures */
    LewisStructure res[MAX_RESONANCE];
    uint8_t  num_res;
    uint8_t  cur_res;      /* currently displayed resonance form */

    uint8_t  central;      /* index of central atom */
    int      total_ve;     /* total valence electrons */
} Molecule;

static Molecule mol;

/* ── Periodic table cursor ── */
static uint8_t cur_row = 0, cur_col = 0;

/* ── Sine/cosine table (fixed-point, 8-bit fraction) ── */
/* sin(angle)*256 for angles 0..11 (30-degree steps mapped to atom positions) */
/* We precompute for up to 12 radial positions */
static const int16_t cos_tbl[12] = {
     256,  222,  128,    0, -128, -222,
    -256, -222, -128,    0,  128,  222
};
static const int16_t sin_tbl[12] = {
       0,  128,  222,  256,  222,  128,
       0, -128, -222, -256, -222, -128
};

/* ── Helper: safe print (clips to screen) ── */
static void safe_print(const char *s, int x, int y)
{
    if (y < 0 || y >= SCR_H - 8) return;
    if (x < 0 || x >= SCR_W) return;
    /* Truncate string if it would go off-screen */
    int max_chars = (SCR_W - x) / 8;
    if (max_chars <= 0) return;
    int len = (int)strlen(s);
    if (len <= max_chars) {
        gfx_PrintStringXY(s, x, y);
    } else {
        char buf[42];
        if (max_chars > 41) max_chars = 41;
        memcpy(buf, s, max_chars);
        buf[max_chars] = '\0';
        gfx_PrintStringXY(buf, x, y);
    }
}

/* ── Helper: integer to string ── */
static char *int_to_str(int val, char *buf)
{
    if (val == 0) { buf[0] = '0'; buf[1] = '\0'; return buf; }
    char tmp[12];
    int i = 0;
    bool neg = false;
    if (val < 0) { neg = true; val = -val; }
    while (val > 0) { tmp[i++] = '0' + (val % 10); val /= 10; }
    int j = 0;
    if (neg) buf[j++] = '-';
    while (i > 0) buf[j++] = tmp[--i];
    buf[j] = '\0';
    return buf;
}

static uint8_t text_color_for_bg(uint8_t bg)
{
    (void)bg;
    return COL_BLACK;
}

/* Safe append helpers for fixed-size UI buffers */
static void append_str(char *dst, size_t cap, const char *src)
{
    size_t len = strlen(dst);
    if (len + 1 >= cap) return;
    strncat(dst, src, cap - len - 1);
}

static void append_int(char *dst, size_t cap, int val)
{
    char nb[12];
    int_to_str(val, nb);
    append_str(dst, cap, nb);
}

/* ── Move cursor to nearest valid cell in a direction ── */
static void move_cursor(int dr, int dc)
{
    int nr = (int)cur_row + dr;
    int nc = (int)cur_col + dc;

    /* Clamp */
    if (nr < 0) nr = 0;
    if (nr >= PT_ROWS) nr = PT_ROWS - 1;
    if (nc < 0) nc = 0;
    if (nc >= PT_COLS) nc = PT_COLS - 1;

    /* If target cell is valid, go there */
    if (pt_grid[nr][nc] != ELEM_NONE) {
        cur_row = nr;
        cur_col = nc;
        return;
    }

    /* Otherwise search in the movement direction for nearest valid cell */
    for (int step = 1; step < PT_COLS; step++) {
        int tr = nr + dr * step;
        int tc = nc + dc * step;
        if (tr < 0 || tr >= PT_ROWS || tc < 0 || tc >= PT_COLS) break;
        if (pt_grid[tr][tc] != ELEM_NONE) {
            cur_row = tr;
            cur_col = tc;
            return;
        }
    }

    /* If direction search fails, scan the target row/col for any element */
    if (dr != 0) {
        /* Moved vertically — scan the target row for closest column */
        int best = -1, best_dist = 999;
        for (int c = 0; c < PT_COLS; c++) {
            if (pt_grid[nr][c] != ELEM_NONE) {
                int d = abs(c - (int)cur_col);
                if (d < best_dist) { best_dist = d; best = c; }
            }
        }
        if (best >= 0) { cur_row = nr; cur_col = best; }
    } else {
        /* Moved horizontally — scan target col for closest row */
        int best = -1, best_dist = 999;
        for (int r = 0; r < PT_ROWS; r++) {
            if (pt_grid[r][nc] != ELEM_NONE) {
                int d = abs(r - (int)cur_row);
                if (d < best_dist) { best_dist = d; best = r; }
            }
        }
        if (best >= 0) { cur_row = best; cur_col = nc; }
    }
}

/* ── Draw the periodic table selector screen ── */
static void draw_periodic_table(void)
{
    gfx_FillScreen(COL_WHITE);

    uint8_t sel_elem = pt_grid[cur_row][cur_col];

    /* ── Top info bar ── */
    gfx_SetColor(COL_BLUE);
    gfx_FillRectangle(0, INFO_Y, SCR_W, INFO_H);
    gfx_SetTextFGColor(COL_WHITE);
    gfx_SetTextBGColor(COL_BLUE);
    safe_print("Lewis Dot Structure Generator", 36, 4);
    if (sel_elem != ELEM_NONE) {
        const Element *e = &elements[sel_elem];
        char buf[32];
        strcpy(buf, "Val e-: ");
        char nb[4];
        int_to_str(e->valence, nb);
        strcat(buf, nb);
        strcat(buf, "  Bonds: ");
        int_to_str(e->bond_cap, nb);
        strcat(buf, nb);
        safe_print(buf, 60, 18);
    }

    /* ── Selected atoms bar ── */
    gfx_SetColor(COL_GRAY);
    gfx_FillRectangle(0, SEL_Y, SCR_W, SEL_H);
    gfx_SetTextFGColor(COL_BLACK);
    gfx_SetTextBGColor(COL_GRAY);

    if (mol.num_atoms > 0) {
        /* Draw formula summary: e.g. [C] [O] [O] */
        int x = 4;
        gfx_SetTextScale(1, 1);
        for (uint8_t i = 0; i < mol.num_atoms && x < SCR_W - 24; i++) {
            const Element *e = &elements[mol.atoms[i].elem];
            gfx_SetColor(e->color);
            gfx_FillRectangle(x, SEL_Y + 3, 18, 14);
            gfx_SetColor(COL_BLACK);
            gfx_Rectangle(x, SEL_Y + 3, 18, 14);
            gfx_SetTextFGColor(text_color_for_bg(e->color));
            gfx_SetTextBGColor(e->color);
            int tx = x + (18 - (int)strlen(e->symbol) * 8) / 2;
            int ty = SEL_Y + 6;
            safe_print(e->symbol, tx, ty);
            x += 20;
        }

        /* Total valence electrons */
        gfx_SetTextFGColor(COL_BLACK);
        gfx_SetTextBGColor(COL_GRAY);
        int total_ve = 0;
        for (uint8_t i = 0; i < mol.num_atoms; i++)
            total_ve += elements[mol.atoms[i].elem].valence;
        total_ve -= mol.charge;

        char buf[32];
        strcpy(buf, "VE: ");
        char nb[6];
        int_to_str(total_ve, nb);
        strcat(buf, nb);
        safe_print(buf, 4, SEL_Y + 22);

        /* Charge display */
        if (mol.charge != 0) {
            strcpy(buf, "Charge: ");
            if (mol.charge > 0) strcat(buf, "+");
            int_to_str(mol.charge, nb);
            strcat(buf, nb);
            safe_print(buf, 80, SEL_Y + 22);
        }

        /* Atom count */
        strcpy(buf, "Atoms: ");
        int_to_str(mol.num_atoms, nb);
        strcat(buf, nb);
        safe_print(buf, 180, SEL_Y + 22);
    } else {
        safe_print("Press [enter] to add atoms", 4, SEL_Y + 6);
        safe_print("[2nd] generate  [mode] quit", 4, SEL_Y + 22);
    }

    /* ── Charge indicator ── */
    {
        char cbuf[12] = "Chg:0";
        if (mol.charge > 0) {
            cbuf[4] = '+';
            int_to_str(mol.charge, cbuf + 5);
        } else if (mol.charge < 0) {
            int_to_str(mol.charge, cbuf + 4);
        }
        gfx_SetTextFGColor(COL_BLACK);
        gfx_SetTextBGColor(COL_GRAY);
        safe_print(cbuf, SCR_W - 56, SEL_Y + 22);
    }

    /* ── Periodic table grid ── */
    int pt_total_w = PT_COLS * PT_CELL_W;
    int pt_x0 = (SCR_W - pt_total_w) / 2; /* center the table */

    for (int r = 0; r < PT_ROWS; r++) {
        for (int c = 0; c < PT_COLS; c++) {
            uint8_t ei = pt_grid[r][c];
            if (ei == ELEM_NONE) continue;

            int cx = pt_x0 + c * PT_CELL_W;
            int cy = PT_Y + r * PT_CELL_H;

            /* Cell background */
            gfx_SetColor(elements[ei].color);
            gfx_FillRectangle(cx + 1, cy + 1, PT_CELL_W - 2, PT_CELL_H - 2);

            /* Highlight cursor */
            if (r == cur_row && c == cur_col) {
                gfx_SetColor(COL_WHITE);
                gfx_Rectangle(cx, cy, PT_CELL_W, PT_CELL_H);
                gfx_Rectangle(cx + 1, cy + 1, PT_CELL_W - 2, PT_CELL_H - 2);
            } else {
                gfx_SetColor(COL_BLACK);
                gfx_Rectangle(cx, cy, PT_CELL_W, PT_CELL_H);
            }

            /* Element symbol */
            gfx_SetTextFGColor(COL_BLACK);
            gfx_SetTextBGColor(elements[ei].color);
            int tx = cx + (PT_CELL_W - (int)strlen(elements[ei].symbol) * 8) / 2;
            int ty = cy + (PT_CELL_H - 8) / 2;
            if (tx >= 0 && tx < SCR_W && ty >= 0 && ty < SCR_H)
                safe_print(elements[ei].symbol, tx, ty);
        }
    }

    /* ── Element info card (in the empty center gap of the table) ── */
    if (sel_elem != ELEM_NONE) {
        const Element *e = &elements[sel_elem];

        /* Card border — colored outer frame, white interior for readability */
        gfx_SetColor(e->color);
        gfx_FillRectangle(CARD_X, CARD_Y, CARD_W, CARD_H);
        gfx_SetColor(COL_WHITE);
        gfx_FillRectangle(CARD_X + 2, CARD_Y + 2, CARD_W - 4, CARD_H - 4);
        gfx_SetColor(COL_BLACK);
        gfx_Rectangle(CARD_X, CARD_Y, CARD_W, CARD_H);

        /* Atomic number — top-left corner */
        gfx_SetTextFGColor(COL_DKGRAY);
        gfx_SetTextBGColor(COL_WHITE);
        {
            char abuf[4];
            int_to_str(e->atomic_num, abuf);
            safe_print(abuf, CARD_X + 5, CARD_Y + 5);
        }

        /* Big symbol — centered, 3x scale, in element color */
        gfx_SetTextScale(3, 3);
        {
            int sym_w = (int)strlen(e->symbol) * 24; /* 8px * 3 scale */
            int sx = CARD_X + (CARD_W - sym_w) / 2;
            int sy = CARD_Y + 16;
            gfx_SetTextFGColor(e->color);
            gfx_SetTextBGColor(COL_WHITE);
            safe_print(e->symbol, sx, sy);
        }
        gfx_SetTextScale(1, 1);

        /* Element name — centered below symbol */
        {
            int name_w = (int)strlen(e->name) * 8;
            int nx = CARD_X + (CARD_W - name_w) / 2;
            gfx_SetTextFGColor(COL_BLACK);
            gfx_SetTextBGColor(COL_WHITE);
            safe_print(e->name, nx, CARD_Y + CARD_H - 22);
        }

        /* Valence electrons — bottom line of card */
        {
            char vbuf[12] = "e-: ";
            char nb[4];
            int_to_str(e->valence, nb);
            strcat(vbuf, nb);
            int vw = (int)strlen(vbuf) * 8;
            gfx_SetTextFGColor(COL_DKGRAY);
            gfx_SetTextBGColor(COL_WHITE);
            safe_print(vbuf, CARD_X + (CARD_W - vw) / 2, CARD_Y + CARD_H - 11);
        }
    }

    /* ── Help text at bottom ── */
    gfx_SetTextFGColor(COL_DKGRAY);
    gfx_SetTextBGColor(COL_WHITE);
    safe_print("[enter]add [del]undo [alpha]chg [2nd]go", 4, SCR_H - 10);
}

/* ── Lewis structure algorithm ── */

/*
 * Helpers for Lewis generation.
 */
static int bond_order_sum(const LewisStructure *ls, uint8_t atom_idx)
{
    int sum = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if (ls->bonds[b].a == atom_idx || ls->bonds[b].b == atom_idx)
            sum += ls->bonds[b].order;
    }
    return sum;
}

static int electrons_on_atom(const LewisStructure *ls, uint8_t atom_idx)
{
    return (ls->lone_pairs[atom_idx] * 2) + (bond_order_sum(ls, atom_idx) * 2);
}

static void recompute_formal_charges(LewisStructure *ls)
{
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        int val = elements[mol.atoms[i].elem].valence;
        int lpe = ls->lone_pairs[i] * 2;
        int bnd_e = bond_order_sum(ls, i);
        ls->formal_charge[i] = (int8_t)(val - lpe - bnd_e);
    }
}

static int formal_charge_sum(const LewisStructure *ls)
{
    int sum = 0;
    for (uint8_t i = 0; i < mol.num_atoms; i++)
        sum += ls->formal_charge[i];
    return sum;
}

static bool structures_equal(const LewisStructure *a, const LewisStructure *b)
{
    if (a->num_bonds != b->num_bonds) return false;
    for (uint8_t i = 0; i < a->num_bonds; i++) {
        if (a->bonds[i].a != b->bonds[i].a) return false;
        if (a->bonds[i].b != b->bonds[i].b) return false;
        if (a->bonds[i].order != b->bonds[i].order) return false;
    }
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        if (a->lone_pairs[i] != b->lone_pairs[i]) return false;
    }
    return true;
}

static bool resonance_exists(const LewisStructure *candidate)
{
    for (uint8_t i = 0; i < mol.num_res; i++) {
        if (structures_equal(candidate, &mol.res[i]))
            return true;
    }
    return false;
}

static bool atom_has_h_neighbor(const LewisStructure *ls, uint8_t atom_idx)
{
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        int other = -1;
        if (ls->bonds[b].a == atom_idx) other = ls->bonds[b].b;
        else if (ls->bonds[b].b == atom_idx) other = ls->bonds[b].a;
        if (other < 0) continue;
        if (mol.atoms[(uint8_t)other].elem == 0)
            return true;
    }
    return false;
}

/* Keep protonated terminal oxygens (X-O-H) out of the resonance swap set. */
static bool is_protonated_terminal_oxygen(const LewisStructure *ls, uint8_t term_idx)
{
    if (mol.atoms[term_idx].elem != 7) return false; /* O terminal */
    if (!atom_has_h_neighbor(ls, term_idx)) return false;

    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if ((ls->bonds[b].a == mol.central && ls->bonds[b].b == term_idx) ||
            (ls->bonds[b].b == mol.central && ls->bonds[b].a == term_idx)) {
            return true;
        }
    }
    return false;
}

static int required_electrons(uint8_t atom_idx, bool is_central)
{
    uint8_t elem_idx = mol.atoms[atom_idx].elem;
    const Element *e = &elements[elem_idx];

    if (elem_idx == 0 || elem_idx == 1) return 2; /* H, He */

    /* Common electron-deficient centers */
    if (is_central && e->group == 2) return 4;   /* Be/Mg */
    if (is_central && e->group == 13) return 6;  /* B/Al */

    return 8;
}

static bool shell_satisfied(uint8_t atom_idx, int electrons, bool is_central)
{
    uint8_t elem_idx = mol.atoms[atom_idx].elem;
    const Element *e = &elements[elem_idx];

    if (elem_idx == 0 || elem_idx == 1) return electrons == 2;

    if (electrons < required_electrons(atom_idx, is_central))
        return false;

    /* Period 2 atoms should not exceed octet */
    if (e->period <= 2 && electrons > 8)
        return false;

    return true;
}

static int bond_limit(uint8_t atom_idx, bool is_central)
{
    uint8_t elem_idx = mol.atoms[atom_idx].elem;
    const Element *e = &elements[elem_idx];
    int limit = e->bond_cap;

    /* Allow ammonium-like cations for period-2 group-15 centers (e.g., NH4+). */
    if (is_central && e->period == 2 && e->group == 15 && mol.charge > 0 && limit < 4) {
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

static bool build_skeleton(LewisStructure *ls, int *ve_pool)
{
    uint8_t remain[MAX_ATOMS];
    bool connected[MAX_ATOMS];
    uint8_t backbone[MAX_ATOMS];
    uint8_t n_backbone = 0;
    uint8_t ordered[MAX_ATOMS];
    bool used[MAX_ATOMS];

    memset(connected, 0, sizeof(connected));
    memset(used, 0, sizeof(used));

    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        uint8_t elem_idx = mol.atoms[i].elem;
        remain[i] = (uint8_t)bond_limit(i, i == mol.central);

        if (i == mol.central) {
            backbone[n_backbone++] = i;
            continue;
        }

        /* Keep highly terminal atoms off the backbone */
        if (elem_idx == 0) continue; /* H */
        if (elements[elem_idx].group == 17) continue; /* halogens */
        if (elements[elem_idx].bond_cap >= 3)
            backbone[n_backbone++] = i;
    }

    if (mol.num_atoms == 1)
        return true;

    if (n_backbone == 0)
        return false;

    /* Order backbone: central first, then higher capacity and lower EN */
    ordered[0] = mol.central;
    used[mol.central] = true;
    uint8_t n_ordered = 1;

    while (n_ordered < n_backbone) {
        int best = -1;
        int best_score = -32768;
        for (uint8_t bi = 0; bi < n_backbone; bi++) {
            uint8_t atom = backbone[bi];
            if (used[atom]) continue;
            const Element *e = &elements[mol.atoms[atom].elem];
            int score = (int)bond_limit(atom, false) * 10 - (int)e->eneg;
            if (score > best_score) {
                best_score = score;
                best = atom;
            }
        }
        if (best < 0) return false;
        ordered[n_ordered++] = (uint8_t)best;
        used[best] = true;
    }

    connected[ordered[0]] = true;
    for (uint8_t i = 1; i < n_ordered; i++) {
        uint8_t a = ordered[i - 1];
        uint8_t b = ordered[i];
        if (!add_single_bond(ls, a, b, ve_pool, remain))
            return false;
        connected[a] = true;
        connected[b] = true;
    }

    /* Attach non-backbone atoms in two passes: heavy atoms first, then H. */
    for (uint8_t pass = 0; pass < 2; pass++) {
        bool target_h = (pass == 1);

        for (uint8_t i = 0; i < mol.num_atoms; i++) {
            if (connected[i]) continue;

            bool is_h = (mol.atoms[i].elem == 0);
            if (is_h != target_h) continue;
            if (remain[i] == 0) return false;

            int best_host = -1;
            int best_score = -32768;
            for (uint8_t j = 0; j < mol.num_atoms; j++) {
                if (!connected[j]) continue;
                if (i == j) continue;
                if (remain[j] == 0) continue;
                if (mol.atoms[j].elem == 0) continue; /* avoid H hosts */

                int score = (int)remain[j] * 10 - (int)elements[mol.atoms[j].elem].eneg;

                /* Count heavy-atom neighbors already attached to this host. */
                int heavy_neighbors = 0;
                for (uint8_t b = 0; b < ls->num_bonds; b++) {
                    int other = -1;
                    if (ls->bonds[b].a == j) other = ls->bonds[b].b;
                    else if (ls->bonds[b].b == j) other = ls->bonds[b].a;
                    if (other >= 0 && mol.atoms[(uint8_t)other].elem != 0)
                        heavy_neighbors++;
                }

                if (!is_h) {
                    /* Encourage heavy atoms to cluster on the likely functional center. */
                    if (j == mol.central) score += 12;
                } else {
                    /* Keep hydrogens off heavily substituted centers when alternatives exist. */
                    score -= heavy_neighbors * 8;
                    if (j == mol.central) score -= 4;
                }

                if (score > best_score) {
                    best_score = score;
                    best_host = j;
                }
            }

            if (best_host < 0) {
                for (uint8_t j = 0; j < mol.num_atoms; j++) {
                    if (!connected[j]) continue;
                    if (i == j) continue;
                    if (remain[j] == 0) continue;
                    best_host = j;
                    break;
                }
            }

            if (best_host < 0) return false;
            if (!add_single_bond(ls, (uint8_t)best_host, i, ve_pool, remain))
                return false;
            connected[i] = true;
        }
    }

    for (uint8_t i = 0; i < mol.num_atoms; i++) {
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
static uint8_t find_central(void)
{
    if (mol.num_atoms == 0) return 0;

    uint8_t counts[NUM_ELEMENTS];
    memset(counts, 0, sizeof(counts));
    for (uint8_t i = 0; i < mol.num_atoms; i++)
        counts[mol.atoms[i].elem]++;

    int best = -1;
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        uint8_t elem_idx = mol.atoms[i].elem;
        if (elem_idx == 0) continue; /* skip H */

        const Element *cand = &elements[elem_idx];
        bool cand_terminal = (cand->group == 17 || cand->bond_cap <= 1);

        if (best < 0) {
            best = i;
            continue;
        }

        const Element *cur = &elements[mol.atoms[best].elem];
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
            counts[elem_idx] > counts[mol.atoms[best].elem]) {
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
static bool generate_structure(LewisStructure *ls)
{
    memset(ls, 0, sizeof(*ls));

    if (mol.num_atoms == 0) return false;
    if (mol.total_ve < 0) return false;
    if ((mol.total_ve & 1) != 0) return false; /* radical species unsupported */

    int ve_pool = mol.total_ve;

    if (!build_skeleton(ls, &ve_pool))
        return false;

    /* Fill terminal atoms first */
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        if (i == mol.central) continue;

        int target = required_electrons(i, false);
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
        ls->lone_pairs[mol.central] = (uint8_t)(ve_pool / 2);
        ve_pool -= ls->lone_pairs[mol.central] * 2;
    }

    if (ve_pool != 0)
        return false;

    /* Promote central bonds to satisfy central shell */
    if (mol.num_atoms > 1) {
        int target_c = required_electrons(mol.central, true);
        uint8_t next_bond = 0;

        for (int pass = 0; pass < MAX_BONDS * 3; pass++) {
            int central_e = electrons_on_atom(ls, mol.central);
            if (central_e >= target_c) break;

            bool promoted = false;
            for (uint8_t scan = 0; scan < ls->num_bonds; scan++) {
                uint8_t b = (next_bond + scan) % ls->num_bonds;
                if (!(ls->bonds[b].a == mol.central || ls->bonds[b].b == mol.central))
                    continue;

                uint8_t term = (ls->bonds[b].a == mol.central) ? ls->bonds[b].b : ls->bonds[b].a;
                if (mol.atoms[term].elem == 0) continue;   /* H can't multiple bond */
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
    if (elements[mol.atoms[mol.central].elem].period >= 3) {
        recompute_formal_charges(ls);
        for (int pass = 0; pass < MAX_BONDS; pass++) {
            if (ls->formal_charge[mol.central] <= 0) break;

            int best_bond = -1;
            int most_negative = 0;
            for (uint8_t b = 0; b < ls->num_bonds; b++) {
                if (!(ls->bonds[b].a == mol.central || ls->bonds[b].b == mol.central))
                    continue;

                uint8_t term = (ls->bonds[b].a == mol.central) ? ls->bonds[b].b : ls->bonds[b].a;
                if (mol.atoms[term].elem == 0) continue;
                if (ls->bonds[b].order >= 3) continue;
                if (ls->lone_pairs[term] == 0) continue;
                if (ls->formal_charge[term] >= 0) continue;
                if (ls->formal_charge[term] < most_negative) {
                    most_negative = ls->formal_charge[term];
                    best_bond = b;
                }
            }

            if (best_bond < 0) break;
            uint8_t term = (ls->bonds[best_bond].a == mol.central) ? ls->bonds[best_bond].b : ls->bonds[best_bond].a;
            ls->bonds[best_bond].order++;
            ls->lone_pairs[term]--;
            recompute_formal_charges(ls);
        }
    }

    recompute_formal_charges(ls);

    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        int electrons = electrons_on_atom(ls, i);
        if (!shell_satisfied(i, electrons, i == mol.central))
            return false;
    }

    if (formal_charge_sum(ls) != mol.charge)
        return false;

    return true;
}

/*
 * Generate resonance structures by exploring multiple-bond shifts between
 * equivalent terminal atoms around the same central atom.
 */
static void generate_resonance(void)
{
    mol.num_res = 0;
    mol.cur_res = 0;

    if (mol.num_atoms == 0) return;

    mol.central = find_central();
    mol.total_ve = 0;
    for (uint8_t i = 0; i < mol.num_atoms; i++)
        mol.total_ve += elements[mol.atoms[i].elem].valence;
    mol.total_ve -= mol.charge;

    if (!generate_structure(&mol.res[0]))
        return;

    mol.num_res = 1;

    for (uint8_t seed_idx = 0; seed_idx < mol.num_res && mol.num_res < MAX_RESONANCE; seed_idx++) {
        LewisStructure seed;
        memcpy(&seed, &mol.res[seed_idx], sizeof(seed));

        for (uint8_t src = 0; src < seed.num_bonds && mol.num_res < MAX_RESONANCE; src++) {
            if (!(seed.bonds[src].a == mol.central || seed.bonds[src].b == mol.central))
                continue;
            if (seed.bonds[src].order <= 1)
                continue;

            uint8_t src_term = (seed.bonds[src].a == mol.central) ? seed.bonds[src].b : seed.bonds[src].a;
            uint8_t src_elem = mol.atoms[src_term].elem;
            if (is_protonated_terminal_oxygen(&seed, src_term)) continue; /* keep X-O-H fixed */
            uint8_t shift = seed.bonds[src].order - 1;

            for (uint8_t dst = 0; dst < seed.num_bonds && mol.num_res < MAX_RESONANCE; dst++) {
                if (dst == src) continue;
                if (!(seed.bonds[dst].a == mol.central || seed.bonds[dst].b == mol.central))
                    continue;

                uint8_t dst_term = (seed.bonds[dst].a == mol.central) ? seed.bonds[dst].b : seed.bonds[dst].a;
                if (mol.atoms[dst_term].elem != src_elem) continue;
                if (mol.atoms[dst_term].elem == 0) continue; /* no H resonance */
                if (is_protonated_terminal_oxygen(&seed, dst_term)) continue; /* avoid X=OH contributors */
                if (seed.bonds[dst].order >= seed.bonds[src].order) continue;
                if ((uint8_t)(seed.bonds[dst].order + shift) > 3) continue;
                if (seed.lone_pairs[dst_term] < shift) continue;

                LewisStructure cand;
                memcpy(&cand, &seed, sizeof(cand));

                cand.bonds[src].order = 1;
                cand.lone_pairs[src_term] += shift;

                cand.bonds[dst].order += shift;
                cand.lone_pairs[dst_term] -= shift;

                recompute_formal_charges(&cand);
                if (formal_charge_sum(&cand) != mol.charge) continue;

                bool valid = true;
                for (uint8_t i = 0; i < mol.num_atoms; i++) {
                    int electrons = electrons_on_atom(&cand, i);
                    if (!shell_satisfied(i, electrons, i == mol.central)) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;

                if (resonance_exists(&cand)) continue;

                memcpy(&mol.res[mol.num_res], &cand, sizeof(cand));
                mol.num_res++;
            }
        }
    }
}


/* ── Draw Lewis structure screen ── */
/* Layout helper: render path-like molecules in a straight horizontal line. */
static bool layout_linear_chain(const LewisStructure *ls, int ax[MAX_ATOMS], int ay[MAX_ATOMS])
{
    if (mol.num_atoms < 3) return false;
    if (ls->num_bonds != mol.num_atoms - 1) return false;

    uint8_t deg[MAX_ATOMS];
    int8_t neigh[MAX_ATOMS][2];
    memset(deg, 0, sizeof(deg));
    memset(neigh, -1, sizeof(neigh));

    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        uint8_t a = ls->bonds[b].a;
        uint8_t c = ls->bonds[b].b;
        if (a >= mol.num_atoms || c >= mol.num_atoms) return false;
        if (deg[a] >= 2 || deg[c] >= 2) return false;
        neigh[a][deg[a]++] = (int8_t)c;
        neigh[c][deg[c]++] = (int8_t)a;
    }

    int endpoints = 0;
    int start = -1;
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
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
    for (uint8_t k = 0; k < mol.num_atoms; k++) {
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
        if (k + 1 < mol.num_atoms && cur < 0) return false;
    }

    int step = (SCR_W - 80) / (mol.num_atoms - 1);
    if (step > BOND_LEN) step = BOND_LEN;
    if (step < 22) step = 22;
    int total_w = step * (mol.num_atoms - 1);
    int x0 = LEWIS_CENTER_X - total_w / 2;

    for (uint8_t k = 0; k < mol.num_atoms; k++) {
        uint8_t idx = order[k];
        ax[idx] = x0 + k * step;
        ay[idx] = LEWIS_CENTER_Y;
    }
    return true;
}

/* Layout helper: place atoms by graph distance from central atom. */
static bool layout_tree_from_central(const LewisStructure *ls, int ax[MAX_ATOMS], int ay[MAX_ATOMS])
{
    if (mol.num_atoms == 0) return false;

    int8_t dist[MAX_ATOMS];
    int8_t parent[MAX_ATOMS];
    uint8_t q[MAX_ATOMS];
    uint8_t qh = 0, qt = 0;

    for (uint8_t i = 0; i < MAX_ATOMS; i++) {
        dist[i] = -1;
        parent[i] = -1;
    }

    dist[mol.central] = 0;
    q[qt++] = mol.central;

    while (qh < qt) {
        uint8_t u = q[qh++];
        for (uint8_t b = 0; b < ls->num_bonds; b++) {
            int v = -1;
            if (ls->bonds[b].a == u) v = ls->bonds[b].b;
            else if (ls->bonds[b].b == u) v = ls->bonds[b].a;
            if (v < 0 || v >= mol.num_atoms) continue;
            if (dist[(uint8_t)v] != -1) continue;
            dist[(uint8_t)v] = dist[u] + 1;
            parent[(uint8_t)v] = (int8_t)u;
            q[qt++] = (uint8_t)v;
        }
    }

    int max_dist = 0;
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        if (dist[i] < 0) return false; /* disconnected */
        if (dist[i] > max_dist) max_dist = dist[i];
    }

    ax[mol.central] = LEWIS_CENTER_X;
    ay[mol.central] = LEWIS_CENTER_Y;

    /* First shell around central */
    uint8_t first[MAX_ATOMS];
    uint8_t n_first = 0;
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
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
        for (uint8_t i = 0; i < mol.num_atoms; i++) {
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
            for (uint8_t j = 0; j < mol.num_atoms; j++) {
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

static void draw_lewis(void)
{
    gfx_FillScreen(COL_WHITE);

    if (mol.num_res == 0 || mol.num_atoms == 0) {
        gfx_SetTextFGColor(COL_BLACK);
        gfx_SetTextBGColor(COL_WHITE);
        safe_print("No valid structure", 80, 120);
        if ((mol.total_ve & 1) != 0) {
            safe_print("Odd electron count", 84, 132);
            safe_print("(radicals unsupported)", 64, 142);
        } else if (mol.charge != 0) {
            char cbuf[24] = "Current charge: ";
            append_int(cbuf, sizeof(cbuf), mol.charge);
            safe_print(cbuf, 76, 132);
            safe_print("[alpha] to change", 88, 142);
        }
        safe_print("[clear] to go back", 80, 156);
        return;
    }

    LewisStructure *ls = &mol.res[mol.cur_res];

    /* ── Header: formula and VE count ── */
    gfx_SetColor(COL_BLUE);
    gfx_FillRectangle(0, 0, SCR_W, 24);
    gfx_SetTextFGColor(COL_WHITE);
    gfx_SetTextBGColor(COL_BLUE);

    /* Build formula string */
    char formula[48] = "";
    /* Count occurrences of each element */
    uint8_t counts[NUM_ELEMENTS];
    memset(counts, 0, sizeof(counts));
    for (uint8_t i = 0; i < mol.num_atoms; i++)
        counts[mol.atoms[i].elem]++;

    for (uint8_t e = 0; e < NUM_ELEMENTS; e++) {
        if (counts[e] == 0) continue;
        append_str(formula, sizeof(formula), elements[e].symbol);
        if (counts[e] > 1) {
            append_int(formula, sizeof(formula), counts[e]);
        }
    }

    if (mol.charge != 0) {
        if (mol.charge > 0) append_str(formula, sizeof(formula), " +");
        else append_str(formula, sizeof(formula), " ");
        append_int(formula, sizeof(formula), mol.charge);
    }

    safe_print(formula, 4, 4);

    /* VE count */
    char vebuf[20] = "VE: ";
    append_int(vebuf, sizeof(vebuf), mol.total_ve);
    safe_print(vebuf, 200, 4);

    /* Formal charge sum */
    {
        int fc_sum = 0;
        for (uint8_t i = 0; i < mol.num_atoms; i++)
            fc_sum += ls->formal_charge[i];
        char fcbuf[24] = "FC: ";
        if (fc_sum > 0) append_str(fcbuf, sizeof(fcbuf), "+");
        append_int(fcbuf, sizeof(fcbuf), fc_sum);
        safe_print(fcbuf, 260, 4);
    }

    /* Resonance indicator */
    if (mol.num_res > 1) {
        char rbuf[20] = "Res: ";
        append_int(rbuf, sizeof(rbuf), mol.cur_res + 1);
        append_str(rbuf, sizeof(rbuf), "/");
        append_int(rbuf, sizeof(rbuf), mol.num_res);
        safe_print(rbuf, 4, 14);
    }

    gfx_SetTextBGColor(COL_WHITE);

    /* ── Compute atom positions ── */
    int ax[MAX_ATOMS], ay[MAX_ATOMS];

    if (mol.num_atoms == 1) {
        ax[0] = LEWIS_CENTER_X;
        ay[0] = LEWIS_CENTER_Y;
    } else if (mol.num_atoms == 2) {
        ax[0] = LEWIS_CENTER_X - BOND_LEN / 2;
        ay[0] = LEWIS_CENTER_Y;
        ax[1] = LEWIS_CENTER_X + BOND_LEN / 2;
        ay[1] = LEWIS_CENTER_Y;
    } else {
        bool has_multiple = false;
        for (uint8_t b = 0; b < ls->num_bonds; b++) {
            if (ls->bonds[b].order > 1) {
                has_multiple = true;
                break;
            }
        }

        if (has_multiple && layout_linear_chain(ls, ax, ay)) {
            /* positions already assigned */
        } else if (layout_tree_from_central(ls, ax, ay)) {
            /* connectivity-aware tree layout */
        } else {
            /* Central atom at center */
            ax[mol.central] = LEWIS_CENTER_X;
            ay[mol.central] = LEWIS_CENTER_Y;

            /* Terminal atoms arranged radially */
            int n_term = mol.num_atoms - 1;
            int term_idx = 0;
            for (uint8_t i = 0; i < mol.num_atoms; i++) {
                if (i == mol.central) continue;
                int angle_idx;
                if (n_term <= 12) {
                    angle_idx = (term_idx * 12) / n_term;
                } else {
                    angle_idx = term_idx % 12;
                }
                ax[i] = LEWIS_CENTER_X + (int)(cos_tbl[angle_idx] * BOND_LEN / 256);
                ay[i] = LEWIS_CENTER_Y + (int)(sin_tbl[angle_idx] * BOND_LEN / 256);
                term_idx++;
            }
        }
    }

    /* ── Draw bonds ── */
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        uint8_t a = ls->bonds[b].a;
        uint8_t bb = ls->bonds[b].b;
        int x1 = ax[a], y1 = ay[a];
        int x2 = ax[bb], y2 = ay[bb];

        gfx_SetColor(COL_BLACK);

        if (ls->bonds[b].order == 1) {
            gfx_Line(x1, y1, x2, y2);
        } else if (ls->bonds[b].order == 2) {
            /* Two parallel lines offset by 2px perpendicular */
            int dx = x2 - x1;
            int dy = y2 - y1;
            /* Perpendicular normalized (approx) */
            int px = -dy;
            int py = dx;
            /* Normalize: find length approx */
            int len = 1;
            if (abs(px) > abs(py)) len = abs(px); else len = abs(py);
            if (len == 0) len = 1;
            int ox = px * 3 / len;
            int oy = py * 3 / len;
            gfx_Line(x1 + ox, y1 + oy, x2 + ox, y2 + oy);
            gfx_Line(x1 - ox, y1 - oy, x2 - ox, y2 - oy);
        } else if (ls->bonds[b].order == 3) {
            int dx = x2 - x1;
            int dy = y2 - y1;
            int px = -dy;
            int py = dx;
            int len = 1;
            if (abs(px) > abs(py)) len = abs(px); else len = abs(py);
            if (len == 0) len = 1;
            int ox = px * 4 / len;
            int oy = py * 4 / len;
            gfx_Line(x1, y1, x2, y2); /* center line */
            gfx_Line(x1 + ox, y1 + oy, x2 + ox, y2 + oy);
            gfx_Line(x1 - ox, y1 - oy, x2 - ox, y2 - oy);
        }
    }

    /* ── Draw atoms (symbols) and lone pairs ── */
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        const Element *e = &elements[mol.atoms[i].elem];
        int sx = ax[i] - (int)strlen(e->symbol) * 4; /* center text */
        int sy = ay[i] - 4;

        /* Background rectangle behind symbol to cover bond lines */
        int tw = (int)strlen(e->symbol) * 8 + 2;
        gfx_SetColor(COL_WHITE);
        gfx_FillRectangle(sx - 1, sy - 1, tw, 10);

        /* Symbol */
        gfx_SetTextFGColor(e->color);
        if (sx >= 0 && sx < SCR_W && sy >= 0 && sy < SCR_H)
            safe_print(e->symbol, sx, sy);

        /* Lone pairs as dots — placed AWAY from bond directions.
         * 4 cardinal slots: 0=top, 1=bottom, 2=left, 3=right.
         * Mark slots occupied by bonds, then fill lone pairs into free slots. */
        if (ls->lone_pairs[i] > 0) {
            gfx_SetColor(e->color);

            /* Determine which cardinal directions have bonds */
            bool slot_used[4] = {false, false, false, false};
            for (uint8_t b = 0; b < ls->num_bonds; b++) {
                int other = -1;
                if (ls->bonds[b].a == i) other = ls->bonds[b].b;
                else if (ls->bonds[b].b == i) other = ls->bonds[b].a;
                else continue;
                int bdx = ax[other] - ax[i];
                int bdy = ay[other] - ay[i];
                /* Pick the cardinal direction closest to bond vector */
                if (abs(bdy) >= abs(bdx)) {
                    /* More vertical */
                    if (bdy < 0) slot_used[0] = true; /* top */
                    else         slot_used[1] = true; /* bottom */
                } else {
                    /* More horizontal */
                    if (bdx < 0) slot_used[2] = true; /* left */
                    else         slot_used[3] = true; /* right */
                }
            }

            /* Collect free slots in priority order: top, bottom, left, right */
            uint8_t free_slots[4];
            uint8_t n_free = 0;
            for (uint8_t s = 0; s < 4; s++) {
                if (!slot_used[s]) free_slots[n_free++] = s;
            }
            /* If not enough free slots, reuse occupied ones */
            for (uint8_t s = 0; s < 4 && n_free < 4; s++) {
                if (slot_used[s]) free_slots[n_free++] = s;
            }

            /* Slot positions and draw directions */
            int slot_x[4] = {ax[i], ax[i], ax[i] - DOT_DIST, ax[i] + DOT_DIST};
            int slot_y[4] = {ay[i] - DOT_DIST, ay[i] + DOT_DIST, ay[i], ay[i]};

            for (uint8_t lp = 0; lp < ls->lone_pairs[i] && lp < 4; lp++) {
                uint8_t s = free_slots[lp];
                int px = slot_x[s];
                int py = slot_y[s];

                if (s < 2) {
                    /* Top or bottom — horizontal dot pair */
                    if (px - 3 >= 0 && px + 3 < SCR_W && py >= 0 && py < SCR_H) {
                        gfx_FillCircle(px - 3, py, DOT_R);
                        gfx_FillCircle(px + 3, py, DOT_R);
                    }
                } else {
                    /* Left or right — vertical dot pair */
                    if (px >= 0 && px < SCR_W && py - 3 >= 0 && py + 3 < SCR_H) {
                        gfx_FillCircle(px, py - 3, DOT_R);
                        gfx_FillCircle(px, py + 3, DOT_R);
                    }
                }
            }
        }

        /* Formal charge — drawn as red text on white background */
        if (ls->formal_charge[i] != 0) {
            char fcbuf[6] = "";
            if (ls->formal_charge[i] > 0) {
                fcbuf[0] = '+';
                int_to_str(ls->formal_charge[i], fcbuf + 1);
            } else {
                int_to_str(ls->formal_charge[i], fcbuf);
            }
            gfx_SetTextFGColor(COL_RED);
            gfx_SetTextBGColor(COL_WHITE);
            int fcx = ax[i] + (int)strlen(e->symbol) * 4 + 2;
            int fcy = ay[i] - 12;
            if (fcx >= 0 && fcx < SCR_W - 16 && fcy >= 0 && fcy < SCR_H)
                safe_print(fcbuf, fcx, fcy);
        }
    }

    /* ── Footer ── */
    gfx_SetTextFGColor(COL_DKGRAY);
    gfx_SetTextBGColor(COL_WHITE);
    if (mol.num_res > 1)
        safe_print("[L/R] resonance  [clear] back", 30, SCR_H - 10);
    else
        safe_print("[clear] back to periodic table", 50, SCR_H - 10);
}

/* ── Main ── */
int main(void)
{
    gfx_Begin();
    gfx_SetDrawBuffer();

    init_pt_grid();

    /* Initialize cursor to Carbon (period 2, group 14 -> row 1, col 13) */
    cur_row = 1;
    cur_col = 13;

    memset(&mol, 0, sizeof(mol));

    bool running = true;
    bool show_lewis = false;
    bool warning = false;
    uint8_t warning_timer = 0;

    /* Key repeat delay */
    uint8_t key_delay = 0;

    /* Set up timer for frame rate */
    timer_Control = TIMER1_ENABLE | TIMER1_32K | TIMER1_0INT | TIMER1_DOWN;
    timer_1_ReloadValue = FRAME_TICKS;
    timer_1_Counter = FRAME_TICKS;

    while (running) {
        /* Wait for frame tick */
        while (timer_1_Counter > 0) {}
        timer_1_Counter = FRAME_TICKS;

        kb_Scan();

        if (show_lewis) {
            /* ── Lewis structure view ── */
            if (kb_Data[6] & kb_Clear) {
                show_lewis = false;
                continue;
            }
            /* Resonance navigation */
            if (mol.num_res > 1) {
                if (kb_Data[7] & kb_Right) {
                    if (key_delay == 0) {
                        mol.cur_res = (mol.cur_res + 1) % mol.num_res;
                        key_delay = 8;
                    }
                }
                if (kb_Data[7] & kb_Left) {
                    if (key_delay == 0) {
                        mol.cur_res = (mol.cur_res == 0) ? mol.num_res - 1 : mol.cur_res - 1;
                        key_delay = 8;
                    }
                }
            }
            if (key_delay > 0) key_delay--;

            draw_lewis();
        } else {
            /* ── Periodic table selector ── */

            /* Quit */
            if (kb_Data[1] & kb_Mode) {
                running = false;
                continue;
            }

            /* Arrow keys with repeat delay */
            if (key_delay == 0) {
                if (kb_Data[7] & kb_Up)    { move_cursor(-1,  0); key_delay = 6; }
                if (kb_Data[7] & kb_Down)  { move_cursor( 1,  0); key_delay = 6; }
                if (kb_Data[7] & kb_Left)  { move_cursor( 0, -1); key_delay = 6; }
                if (kb_Data[7] & kb_Right) { move_cursor( 0,  1); key_delay = 6; }

                /* Enter — add atom */
                if (kb_Data[6] & kb_Enter) {
                    uint8_t ei = pt_grid[cur_row][cur_col];
                    if (ei != ELEM_NONE && mol.num_atoms < MAX_ATOMS) {
                        /* Check heavy atom limit */
                        int heavy = 0;
                        for (uint8_t i = 0; i < mol.num_atoms; i++)
                            if (mol.atoms[i].elem != 0) heavy++;
                        if (ei != 0 && heavy >= MAX_HEAVY) {
                            warning = true;
                            warning_timer = 40;
                        } else {
                            mol.atoms[mol.num_atoms].elem = ei;
                            mol.num_atoms++;
                        }
                    }
                    key_delay = 8;
                }

                /* Del — remove last atom */
                if (kb_Data[1] & kb_Del) {
                    if (mol.num_atoms > 0) {
                        mol.num_atoms--;
                        if (mol.num_atoms == 0)
                            mol.charge = 0;
                    }
                    key_delay = 8;
                }

                /* Alpha — cycle charge */
                if (kb_Data[2] & kb_Alpha) {
                    /* Cycle: 0 -> +1 -> +2 -> -1 -> -2 -> 0 */
                    if (mol.charge == 0) mol.charge = 1;
                    else if (mol.charge == 1) mol.charge = 2;
                    else if (mol.charge == 2) mol.charge = -1;
                    else if (mol.charge == -1) mol.charge = -2;
                    else mol.charge = 0;
                    key_delay = 8;
                }

                /* 2nd — generate Lewis structure */
                if (kb_Data[1] & kb_2nd) {
                    if (mol.num_atoms >= 1) {
                        generate_resonance();
                        show_lewis = true;
                    }
                    key_delay = 10;
                }
            }
            if (key_delay > 0) key_delay--;

            draw_periodic_table();

            /* Warning overlay */
            if (warning && warning_timer > 0) {
                gfx_SetColor(COL_RED);
                gfx_FillRectangle(40, 100, 240, 30);
                gfx_SetColor(COL_BLACK);
                gfx_Rectangle(40, 100, 240, 30);
                gfx_SetTextScale(1, 1);
                gfx_SetTextFGColor(COL_BLACK);
                gfx_SetTextBGColor(COL_RED);
                safe_print("Max 6 heavy atoms!", 72, 110);
                warning_timer--;
                if (warning_timer == 0) warning = false;
            }
        }

        gfx_SwapDraw();
    }

    gfx_End();
    return 0;
}
