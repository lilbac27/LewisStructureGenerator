#include "ui_vsepr.h"

#include <graphx.h>
#include <stdlib.h>
#include <string.h>

#include "lewis_engine.h"
#include "ui_theme.h"
#include "ui_text.h"

#define VSEPR_CARD_X 196
#define VSEPR_CARD_Y 28
#define VSEPR_CARD_W 120
#define VSEPR_CARD_H 126
#define VSEPR_HIDE_OVERLAP_SCORE 1000

typedef struct {
    int x;
    int y;
    int w;
    int h;
} Rect;

static int rect_intersection_area(const Rect *a, const Rect *b)
{
    if (a->w <= 0 || a->h <= 0 || b->w <= 0 || b->h <= 0) {
        return 0;
    }

    int left = (a->x > b->x) ? a->x : b->x;
    int right = ((a->x + a->w) < (b->x + b->w)) ? (a->x + a->w) : (b->x + b->w);
    int top = (a->y > b->y) ? a->y : b->y;
    int bottom = ((a->y + a->h) < (b->y + b->h)) ? (a->y + a->h) : (b->y + b->h);

    int w = right - left;
    int h = bottom - top;
    if (w <= 0 || h <= 0) {
        return 0;
    }
    return w * h;
}

static Rect line_bounds(int x1, int y1, int x2, int y2, int pad)
{
    int min_x = (x1 < x2) ? x1 : x2;
    int max_x = (x1 > x2) ? x1 : x2;
    int min_y = (y1 < y2) ? y1 : y2;
    int max_y = (y1 > y2) ? y1 : y2;

    Rect r;
    r.x = min_x - pad;
    r.y = min_y - pad;
    r.w = (max_x - min_x) + 1 + (pad * 2);
    r.h = (max_y - min_y) + 1 + (pad * 2);
    return r;
}

static int lone_pair_overlap_score(const LewisStructure *ls, const int ax[MAX_ATOMS], const int ay[MAX_ATOMS], uint8_t atom_idx, const Rect *panel)
{
    if (ls->lone_pairs[atom_idx] == 0) {
        return 0;
    }

    bool slot_used[4] = { false, false, false, false };
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        int other = -1;
        if (ls->bonds[b].a == atom_idx) other = ls->bonds[b].b;
        else if (ls->bonds[b].b == atom_idx) other = ls->bonds[b].a;
        else continue;

        int bdx = ax[(uint8_t)other] - ax[atom_idx];
        int bdy = ay[(uint8_t)other] - ay[atom_idx];
        if (abs(bdy) >= abs(bdx)) {
            if (bdy < 0) slot_used[0] = true;
            else         slot_used[1] = true;
        } else {
            if (bdx < 0) slot_used[2] = true;
            else         slot_used[3] = true;
        }
    }

    uint8_t free_slots[4];
    uint8_t n_free = 0;
    for (uint8_t s = 0; s < 4; s++) {
        if (!slot_used[s]) free_slots[n_free++] = s;
    }
    for (uint8_t s = 0; s < 4 && n_free < 4; s++) {
        if (slot_used[s]) free_slots[n_free++] = s;
    }

    int slot_x[4] = { ax[atom_idx], ax[atom_idx], ax[atom_idx] - DOT_DIST, ax[atom_idx] + DOT_DIST };
    int slot_y[4] = { ay[atom_idx] - DOT_DIST, ay[atom_idx] + DOT_DIST, ay[atom_idx], ay[atom_idx] };

    int overlap = 0;
    for (uint8_t lp = 0; lp < ls->lone_pairs[atom_idx] && lp < 4; lp++) {
        uint8_t s = free_slots[lp];
        int px = slot_x[s];
        int py = slot_y[s];

        if (s < 2) {
            Rect dot_a = { px - 3 - DOT_R, py - DOT_R, (DOT_R * 2) + 1, (DOT_R * 2) + 1 };
            Rect dot_b = { px + 3 - DOT_R, py - DOT_R, (DOT_R * 2) + 1, (DOT_R * 2) + 1 };
            overlap += rect_intersection_area(&dot_a, panel);
            overlap += rect_intersection_area(&dot_b, panel);
        } else {
            Rect dot_a = { px - DOT_R, py - 3 - DOT_R, (DOT_R * 2) + 1, (DOT_R * 2) + 1 };
            Rect dot_b = { px - DOT_R, py + 3 - DOT_R, (DOT_R * 2) + 1, (DOT_R * 2) + 1 };
            overlap += rect_intersection_area(&dot_a, panel);
            overlap += rect_intersection_area(&dot_b, panel);
        }
    }

    return overlap;
}

static int card_overlap_score(const Molecule *mol, const LewisStructure *ls, const int ax[MAX_ATOMS], const int ay[MAX_ATOMS], const Rect *panel)
{
    int overlap = 0;

    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if (ls->bonds[b].a == mol->central || ls->bonds[b].b == mol->central) {
            uint8_t a = ls->bonds[b].a;
            uint8_t c = ls->bonds[b].b;
            Rect bond = line_bounds(ax[a], ay[a], ax[c], ay[c], 2);
            overlap += rect_intersection_area(&bond, panel);
        }
        if (overlap >= VSEPR_HIDE_OVERLAP_SCORE) return overlap;
    }

    for (uint8_t i = 0; i < mol->num_atoms; i++) {
        const Element *e = &elements[mol->atoms[i].elem];
        int sym_w = (int)strlen(e->symbol) * 8 + 2;
        int sx = ax[i] - (int)strlen(e->symbol) * 4;
        int sy = ay[i] - 4;

        Rect symbol = { sx - 1, sy - 1, sym_w, 10 };
        overlap += rect_intersection_area(&symbol, panel);
        if (overlap >= VSEPR_HIDE_OVERLAP_SCORE) return overlap;

        overlap += lone_pair_overlap_score(ls, ax, ay, i, panel);
        if (overlap >= VSEPR_HIDE_OVERLAP_SCORE) return overlap;

        if (ls->formal_charge[i] != 0) {
            char fcbuf[6] = "";
            if (ls->formal_charge[i] > 0) {
                fcbuf[0] = '+';
                int_to_str(ls->formal_charge[i], fcbuf + 1);
            } else {
                int_to_str(ls->formal_charge[i], fcbuf);
            }

            int fcx = ax[i] + (int)strlen(e->symbol) * 4 + 2;
            int fcy = ay[i] - 12;
            Rect fc = { fcx, fcy, (int)strlen(fcbuf) * 8, 8 };
            overlap += rect_intersection_area(&fc, panel);
            if (overlap >= VSEPR_HIDE_OVERLAP_SCORE) return overlap;
        }
    }

    return overlap;
}

bool draw_vsepr_info_card(const Molecule *mol, const LewisStructure *ls, const int ax[MAX_ATOMS], const int ay[MAX_ATOMS], bool force_visible)
{
    if (mol == NULL || ls == NULL || ax == NULL || ay == NULL) {
        return false;
    }
    if (mol->num_atoms == 0 || mol->central >= mol->num_atoms) {
        return false;
    }

    Rect panel = { VSEPR_CARD_X, VSEPR_CARD_Y, VSEPR_CARD_W, VSEPR_CARD_H };
    if (!force_visible && card_overlap_score(mol, ls, ax, ay, &panel) >= VSEPR_HIDE_OVERLAP_SCORE) {
        return false;
    }

    VseprInfo info;
    bool has_row = lewis_get_vsepr_info(mol, ls, &info);
    const char *ep_geometry = (has_row && info.ep_geometry != NULL) ? info.ep_geometry : "N/A";
    const char *shape = (has_row && info.shape != NULL) ? info.shape : "N/A";
    const char *hybrid = (has_row && info.hybridization != NULL) ? info.hybridization : "N/A";
    const char *bond_angle = (has_row && info.bond_angle != NULL) ? info.bond_angle : "N/A";

    uint8_t sigma_bonds = ls->num_bonds;
    uint8_t pi_bonds = 0;
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        if (ls->bonds[b].order > 1) {
            pi_bonds = (uint8_t)(pi_bonds + (ls->bonds[b].order - 1));
        }
    }

    gfx_SetColor(UI_BORDER);
    gfx_FillRectangle(panel.x, panel.y, panel.w, panel.h);
    gfx_SetColor(UI_SURFACE);
    gfx_FillRectangle(panel.x + 1, panel.y + 1, panel.w - 2, panel.h - 2);
    gfx_SetColor(UI_SELECTED_BG);
    gfx_FillRectangle(panel.x + 1, panel.y + 1, panel.w - 2, 12);
    gfx_SetColor(UI_BORDER);
    gfx_Rectangle(panel.x, panel.y, panel.w, panel.h);

    gfx_SetTextFGColor(UI_SELECTED_TEXT);
    gfx_SetTextBGColor(UI_SELECTED_BG);
    safe_print("VSEPR", panel.x + 4, panel.y + 3);

    gfx_SetTextFGColor(UI_TEXT);
    gfx_SetTextBGColor(UI_SURFACE);

    char pair_buf[28] = "EP:";
    append_int(pair_buf, sizeof(pair_buf), info.valence_pairs);
    append_str(pair_buf, sizeof(pair_buf), " BP:");
    append_int(pair_buf, sizeof(pair_buf), info.bond_pairs);
    append_str(pair_buf, sizeof(pair_buf), " LP:");
    append_int(pair_buf, sizeof(pair_buf), info.lone_pairs);
    safe_print(pair_buf, panel.x + 4, panel.y + 16);

    char bond_buf[28] = "Sig:";
    append_int(bond_buf, sizeof(bond_buf), sigma_bonds);
    append_str(bond_buf, sizeof(bond_buf), " Pi:");
    append_int(bond_buf, sizeof(bond_buf), pi_bonds);
    safe_print(bond_buf, panel.x + 4, panel.y + 26);

    if (!has_row) {
        gfx_SetColor(UI_ALERT_BG);
        gfx_FillRectangle(panel.x + 3, panel.y + 37, panel.w - 6, 10);
        gfx_SetTextFGColor(UI_ALERT_TEXT);
        gfx_SetTextBGColor(UI_ALERT_BG);
        safe_print("No table match", panel.x + 6, panel.y + 38);

        gfx_SetTextFGColor(UI_TEXT);
        gfx_SetTextBGColor(UI_SURFACE);

        safe_print("E-Geom: N/A", panel.x + 4, panel.y + 52);
        safe_print("Shape: N/A", panel.x + 4, panel.y + 68);
        safe_print("Hyb: N/A", panel.x + 4, panel.y + 84);
        safe_print("Angle:", panel.x + 4, panel.y + 98);
        safe_print_wrapped("N/A", panel.x + 4, panel.y + 106, panel.w - 8, 2);
        return true;
    }

    gfx_SetTextFGColor(UI_TEXT);
    safe_print("E-Geom:", panel.x + 4, panel.y + 36);
    safe_print_wrapped(ep_geometry, panel.x + 4, panel.y + 46, panel.w - 8, 2);

    safe_print("Shape:", panel.x + 4, panel.y + 60);
    safe_print_wrapped(shape, panel.x + 4, panel.y + 70, panel.w - 8, 2);

    char hyb_buf[24] = "Hyb: ";
    append_str(hyb_buf, sizeof(hyb_buf), hybrid);
    safe_print(hyb_buf, panel.x + 4, panel.y + 86);

    safe_print("Angle:", panel.x + 4, panel.y + 98);
    safe_print_wrapped(bond_angle, panel.x + 4, panel.y + 106, panel.w - 8, 2);

    return true;
}
