#include "ui_periodic.h"

#include <graphx.h>
#include <stdlib.h>
#include <string.h>

#include "ui_text.h"

void move_cursor(uint8_t *cur_row, uint8_t *cur_col, int dr, int dc)
{
    int nr = (int)(*cur_row) + dr;
    int nc = (int)(*cur_col) + dc;

    if (nr < 0) nr = 0;
    if (nr >= PT_ROWS) nr = PT_ROWS - 1;
    if (nc < 0) nc = 0;
    if (nc >= PT_COLS) nc = PT_COLS - 1;

    if (pt_grid[nr][nc] != ELEM_NONE) {
        *cur_row = (uint8_t)nr;
        *cur_col = (uint8_t)nc;
        return;
    }

    for (int step = 1; step < PT_COLS; step++) {
        int tr = nr + dr * step;
        int tc = nc + dc * step;
        if (tr < 0 || tr >= PT_ROWS || tc < 0 || tc >= PT_COLS) break;
        if (pt_grid[tr][tc] != ELEM_NONE) {
            *cur_row = (uint8_t)tr;
            *cur_col = (uint8_t)tc;
            return;
        }
    }

    if (dr != 0) {
        int best = -1;
        int best_dist = 999;
        for (int c = 0; c < PT_COLS; c++) {
            if (pt_grid[nr][c] != ELEM_NONE) {
                int d = abs(c - (int)(*cur_col));
                if (d < best_dist) {
                    best_dist = d;
                    best = c;
                }
            }
        }
        if (best >= 0) {
            *cur_row = (uint8_t)nr;
            *cur_col = (uint8_t)best;
        }
    } else {
        int best = -1;
        int best_dist = 999;
        for (int r = 0; r < PT_ROWS; r++) {
            if (pt_grid[r][nc] != ELEM_NONE) {
                int d = abs(r - (int)(*cur_row));
                if (d < best_dist) {
                    best_dist = d;
                    best = r;
                }
            }
        }
        if (best >= 0) {
            *cur_row = (uint8_t)best;
            *cur_col = (uint8_t)nc;
        }
    }
}

void draw_periodic_table(const Molecule *mol, uint8_t cur_row, uint8_t cur_col)
{
    gfx_FillScreen(COL_WHITE);

    uint8_t sel_elem = pt_grid[cur_row][cur_col];

    /* Top info bar */
    gfx_SetColor(COL_BLUE);
    gfx_FillRectangle(0, INFO_Y, SCR_W, INFO_H);
    gfx_SetTextFGColor(COL_WHITE);
    gfx_SetTextBGColor(COL_BLUE);
    safe_print("Lewis Dot Structure Generator", 36, 4);

    if (sel_elem != ELEM_NONE) {
        const Element *e = &elements[sel_elem];
        char buf[32] = "";
        append_str(buf, sizeof(buf), "Val e-: ");
        append_int(buf, sizeof(buf), e->valence);
        append_str(buf, sizeof(buf), "  Bonds: ");
        append_int(buf, sizeof(buf), e->bond_cap);
        safe_print(buf, 60, 18);
    }

    /* Selected atoms bar */
    gfx_SetColor(COL_GRAY);
    gfx_FillRectangle(0, SEL_Y, SCR_W, SEL_H);
    gfx_SetTextFGColor(COL_BLACK);
    gfx_SetTextBGColor(COL_GRAY);

    if (mol->num_atoms > 0) {
        int x = 4;
        gfx_SetTextScale(1, 1);
        for (uint8_t i = 0; i < mol->num_atoms && x < SCR_W - 24; i++) {
            const Element *e = &elements[mol->atoms[i].elem];
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

        gfx_SetTextFGColor(COL_BLACK);
        gfx_SetTextBGColor(COL_GRAY);

        int total_ve = 0;
        for (uint8_t i = 0; i < mol->num_atoms; i++) {
            total_ve += elements[mol->atoms[i].elem].valence;
        }
        total_ve -= mol->charge;

        char buf[32] = "VE: ";
        append_int(buf, sizeof(buf), total_ve);
        safe_print(buf, 4, SEL_Y + 22);

        if (mol->charge != 0) {
            char cbuf[32] = "Charge: ";
            if (mol->charge > 0) append_str(cbuf, sizeof(cbuf), "+");
            append_int(cbuf, sizeof(cbuf), mol->charge);
            safe_print(cbuf, 80, SEL_Y + 22);
        }

        char abuf[24] = "Atoms: ";
        append_int(abuf, sizeof(abuf), mol->num_atoms);
        safe_print(abuf, 180, SEL_Y + 22);
    } else {
        safe_print("Press [enter] to add atoms", 4, SEL_Y + 6);
        safe_print("[2nd] generate  [mode] quit", 4, SEL_Y + 22);
    }

    /* Charge indicator */
    {
        char cbuf[12] = "Chg:0";
        if (mol->charge > 0) {
            cbuf[4] = '+';
            int_to_str(mol->charge, cbuf + 5);
        } else if (mol->charge < 0) {
            int_to_str(mol->charge, cbuf + 4);
        }
        gfx_SetTextFGColor(COL_BLACK);
        gfx_SetTextBGColor(COL_GRAY);
        safe_print(cbuf, SCR_W - 56, SEL_Y + 22);
    }

    /* Periodic table grid */
    int pt_total_w = PT_COLS * PT_CELL_W;
    int pt_x0 = (SCR_W - pt_total_w) / 2;

    for (int r = 0; r < PT_ROWS; r++) {
        for (int c = 0; c < PT_COLS; c++) {
            uint8_t ei = pt_grid[r][c];
            if (ei == ELEM_NONE) continue;

            int cx = pt_x0 + c * PT_CELL_W;
            int cy = PT_Y + r * PT_CELL_H;

            gfx_SetColor(elements[ei].color);
            gfx_FillRectangle(cx + 1, cy + 1, PT_CELL_W - 2, PT_CELL_H - 2);

            if (r == cur_row && c == cur_col) {
                gfx_SetColor(COL_WHITE);
                gfx_Rectangle(cx, cy, PT_CELL_W, PT_CELL_H);
                gfx_Rectangle(cx + 1, cy + 1, PT_CELL_W - 2, PT_CELL_H - 2);
            } else {
                gfx_SetColor(COL_BLACK);
                gfx_Rectangle(cx, cy, PT_CELL_W, PT_CELL_H);
            }

            gfx_SetTextFGColor(COL_BLACK);
            gfx_SetTextBGColor(elements[ei].color);
            int tx = cx + (PT_CELL_W - (int)strlen(elements[ei].symbol) * 8) / 2;
            int ty = cy + (PT_CELL_H - 8) / 2;
            if (tx >= 0 && tx < SCR_W && ty >= 0 && ty < SCR_H) {
                safe_print(elements[ei].symbol, tx, ty);
            }
        }
    }

    /* Element info card */
    if (sel_elem != ELEM_NONE) {
        const Element *e = &elements[sel_elem];

        gfx_SetColor(e->color);
        gfx_FillRectangle(CARD_X, CARD_Y, CARD_W, CARD_H);
        gfx_SetColor(COL_WHITE);
        gfx_FillRectangle(CARD_X + 2, CARD_Y + 2, CARD_W - 4, CARD_H - 4);
        gfx_SetColor(COL_BLACK);
        gfx_Rectangle(CARD_X, CARD_Y, CARD_W, CARD_H);

        gfx_SetTextFGColor(COL_DKGRAY);
        gfx_SetTextBGColor(COL_WHITE);
        {
            char abuf[4];
            int_to_str(e->atomic_num, abuf);
            safe_print(abuf, CARD_X + 5, CARD_Y + 5);
        }

        gfx_SetTextScale(3, 3);
        {
            int sym_w = (int)strlen(e->symbol) * 24;
            int sx = CARD_X + (CARD_W - sym_w) / 2;
            int sy = CARD_Y + 16;
            gfx_SetTextFGColor(e->color);
            gfx_SetTextBGColor(COL_WHITE);
            safe_print(e->symbol, sx, sy);
        }
        gfx_SetTextScale(1, 1);

        {
            int name_w = (int)strlen(e->name) * 8;
            int nx = CARD_X + (CARD_W - name_w) / 2;
            gfx_SetTextFGColor(COL_BLACK);
            gfx_SetTextBGColor(COL_WHITE);
            safe_print(e->name, nx, CARD_Y + CARD_H - 22);
        }

        {
            char vbuf[12] = "e-: ";
            append_int(vbuf, sizeof(vbuf), e->valence);
            int vw = (int)strlen(vbuf) * 8;
            gfx_SetTextFGColor(COL_DKGRAY);
            gfx_SetTextBGColor(COL_WHITE);
            safe_print(vbuf, CARD_X + (CARD_W - vw) / 2, CARD_Y + CARD_H - 11);
        }
    }

    gfx_SetTextFGColor(COL_DKGRAY);
    gfx_SetTextBGColor(COL_WHITE);
    safe_print("[enter]add [del]undo [alpha]chg [2nd]go", 4, SCR_H - 10);
}
