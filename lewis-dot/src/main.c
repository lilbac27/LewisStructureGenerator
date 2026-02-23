/*
 * Lewis Dot Structure Generator for TI-84 Plus CE
 * Generates and displays Lewis dot structures for small molecules.
 * Supports resonance structures with [left]/[right] navigation.
 *
 * Controls (Periodic Table screen):
 *   Arrow keys  - move cursor on periodic table
 *   [enter]     - add highlighted element to molecule
 *   [del]       - remove last added atom
 *   [alpha]     - cycle charge (0, +1, +2, -1, -2)
 *   [2nd]       - generate Lewis structure
 *   [mode]      - quit
 *
 * Controls (Lewis Structure screen):
 *   [left]/[right] - cycle resonance structures
 *   [alpha]        - cycle charge and regenerate
 *   [2nd]          - hide/show VSEPR card (force-show when auto-hidden)
 *   [clear]        - return to periodic table
 */

#include <graphx.h>
#include <keypadc.h>
#include <sys/timers.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "layout.h"
#include "lewis_engine.h"
#include "lewis_model.h"
#include "ui_periodic.h"
#include "ui_theme.h"
#include "ui_text.h"
#include "ui_vsepr.h"

static Molecule mol;
static uint8_t cur_row = 0;
static uint8_t cur_col = 0;

static void cycle_charge(Molecule *m)
{
    /* Cycle: 0 -> +1 -> +2 -> -1 -> -2 -> 0 */
    if (m->charge == 0) m->charge = 1;
    else if (m->charge == 1) m->charge = 2;
    else if (m->charge == 2) m->charge = -1;
    else if (m->charge == -1) m->charge = -2;
    else m->charge = 0;
}

static bool draw_lewis(bool vsepr_force_visible, bool vsepr_card_enabled)
{
    gfx_FillScreen(UI_BG);

    if (mol.num_res == 0 || mol.num_atoms == 0) {
        gfx_SetTextFGColor(UI_TEXT);
        gfx_SetTextBGColor(UI_BG);
        safe_print("No valid structure", 80, 114);
        safe_print(invalid_reason_message(mol.invalid_reason), 24, 128);
        safe_print("[alpha] change charge  [clear] back", 28, 148);
        return false;
    }

    LewisStructure *ls = &mol.res[mol.cur_res];

    gfx_SetColor(UI_SELECTED_BG);
    gfx_FillRectangle(0, 0, SCR_W, 24);
    gfx_SetTextFGColor(UI_SELECTED_TEXT);
    gfx_SetTextBGColor(UI_SELECTED_BG);

    /* Build formula string */
    char formula[48] = "";
    uint8_t counts[NUM_ELEMENTS];
    memset(counts, 0, sizeof(counts));
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        counts[mol.atoms[i].elem]++;
    }

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

    char vebuf[20] = "VE: ";
    append_int(vebuf, sizeof(vebuf), mol.total_ve);
    safe_print(vebuf, 200, 4);

    /* Formal charge sum */
    {
        int fc_sum = 0;
        for (uint8_t i = 0; i < mol.num_atoms; i++) {
            fc_sum += ls->formal_charge[i];
        }
        char fcbuf[24] = "FC: ";
        if (fc_sum > 0) append_str(fcbuf, sizeof(fcbuf), "+");
        append_int(fcbuf, sizeof(fcbuf), fc_sum);
        safe_print(fcbuf, 260, 4);
    }

    if (mol.num_res > 1) {
        char rbuf[20] = "Res: ";
        append_int(rbuf, sizeof(rbuf), mol.cur_res + 1);
        append_str(rbuf, sizeof(rbuf), "/");
        append_int(rbuf, sizeof(rbuf), mol.num_res);
        safe_print(rbuf, 4, 14);
    }

    gfx_SetTextBGColor(UI_BG);

    int ax[MAX_ATOMS];
    int ay[MAX_ATOMS];

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

        if (has_multiple && layout_linear_chain(&mol, ls, ax, ay)) {
            /* Assigned in helper */
        } else if (layout_tree_from_central(&mol, ls, ax, ay)) {
            /* Assigned in helper */
        } else {
            /* Fallback: simple radial arrangement around central */
            static const int16_t cos_tbl[12] = {
                 256,  222,  128,    0, -128, -222,
                -256, -222, -128,    0,  128,  222
            };
            static const int16_t sin_tbl[12] = {
                   0,  128,  222,  256,  222,  128,
                   0, -128, -222, -256, -222, -128
            };

            ax[mol.central] = LEWIS_CENTER_X;
            ay[mol.central] = LEWIS_CENTER_Y;

            int n_term = mol.num_atoms - 1;
            int term_idx = 0;
            for (uint8_t i = 0; i < mol.num_atoms; i++) {
                if (i == mol.central) continue;
                int angle_idx = (n_term <= 12) ? ((term_idx * 12) / n_term) : (term_idx % 12);
                ax[i] = LEWIS_CENTER_X + (int)(cos_tbl[angle_idx] * BOND_LEN / 256);
                ay[i] = LEWIS_CENTER_Y + (int)(sin_tbl[angle_idx] * BOND_LEN / 256);
                term_idx++;
            }
        }
    }

    /* Draw bonds */
    for (uint8_t b = 0; b < ls->num_bonds; b++) {
        uint8_t a = ls->bonds[b].a;
        uint8_t bb = ls->bonds[b].b;
        int x1 = ax[a];
        int y1 = ay[a];
        int x2 = ax[bb];
        int y2 = ay[bb];

        gfx_SetColor(COL_BLACK);

        if (ls->bonds[b].order == 1) {
            gfx_Line(x1, y1, x2, y2);
        } else if (ls->bonds[b].order == 2) {
            int dx = x2 - x1;
            int dy = y2 - y1;
            int px = -dy;
            int py = dx;
            int len = (abs(px) > abs(py)) ? abs(px) : abs(py);
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
            int len = (abs(px) > abs(py)) ? abs(px) : abs(py);
            if (len == 0) len = 1;
            int ox = px * 4 / len;
            int oy = py * 4 / len;
            gfx_Line(x1, y1, x2, y2);
            gfx_Line(x1 + ox, y1 + oy, x2 + ox, y2 + oy);
            gfx_Line(x1 - ox, y1 - oy, x2 - ox, y2 - oy);
        }
    }

    /* Draw atoms (symbols), lone pairs, and formal charges */
    for (uint8_t i = 0; i < mol.num_atoms; i++) {
        const Element *e = &elements[mol.atoms[i].elem];
        int sx = ax[i] - (int)strlen(e->symbol) * 4;
        int sy = ay[i] - 4;

        int tw = (int)strlen(e->symbol) * 8 + 2;
        gfx_SetColor(UI_BG);
        gfx_FillRectangle(sx - 1, sy - 1, tw, 10);

        gfx_SetTextFGColor(UI_TEXT);
        if (sx >= 0 && sx < SCR_W && sy >= 0 && sy < SCR_H) {
            safe_print(e->symbol, sx, sy);
        }

        if (ls->lone_pairs[i] > 0) {
            gfx_SetColor(UI_TEXT);

            bool slot_used[4] = { false, false, false, false };
            for (uint8_t b = 0; b < ls->num_bonds; b++) {
                int other = -1;
                if (ls->bonds[b].a == i) other = ls->bonds[b].b;
                else if (ls->bonds[b].b == i) other = ls->bonds[b].a;
                else continue;
                int bdx = ax[other] - ax[i];
                int bdy = ay[other] - ay[i];
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

            int slot_x[4] = { ax[i], ax[i], ax[i] - DOT_DIST, ax[i] + DOT_DIST };
            int slot_y[4] = { ay[i] - DOT_DIST, ay[i] + DOT_DIST, ay[i], ay[i] };

            for (uint8_t lp = 0; lp < ls->lone_pairs[i] && lp < 4; lp++) {
                uint8_t s = free_slots[lp];
                int px = slot_x[s];
                int py = slot_y[s];

                if (s < 2) {
                    if (px - 3 >= 0 && px + 3 < SCR_W && py >= 0 && py < SCR_H) {
                        gfx_FillCircle(px - 3, py, DOT_R);
                        gfx_FillCircle(px + 3, py, DOT_R);
                    }
                } else {
                    if (px >= 0 && px < SCR_W && py - 3 >= 0 && py + 3 < SCR_H) {
                        gfx_FillCircle(px, py - 3, DOT_R);
                        gfx_FillCircle(px, py + 3, DOT_R);
                    }
                }
            }
        }

        if (ls->formal_charge[i] != 0) {
            char fcbuf[6] = "";
            if (ls->formal_charge[i] > 0) {
                fcbuf[0] = '+';
                int_to_str(ls->formal_charge[i], fcbuf + 1);
            } else {
                int_to_str(ls->formal_charge[i], fcbuf);
            }
            gfx_SetTextFGColor(UI_TEXT);
            gfx_SetTextBGColor(UI_BG);
            int fcx = ax[i] + (int)strlen(e->symbol) * 4 + 2;
            int fcy = ay[i] - 12;
            if (fcx >= 0 && fcx < SCR_W - 16 && fcy >= 0 && fcy < SCR_H) {
                safe_print(fcbuf, fcx, fcy);
            }
        }
    }

    bool card_drawn = false;
    if (vsepr_card_enabled) {
        card_drawn = draw_vsepr_info_card(&mol, ls, ax, ay, vsepr_force_visible);
    }

    gfx_SetTextFGColor(UI_TEXT);
    gfx_SetTextBGColor(UI_BG);
    if (mol.num_res > 1) {
        if (card_drawn) {
            safe_print("[L/R] [alpha]chg [2nd]hide [clear]back", 0, SCR_H - 10);
        } else {
            safe_print("[L/R] [alpha]chg [2nd]show [clear]back", 0, SCR_H - 10);
        }
    } else {
        if (card_drawn) {
            safe_print("[alpha]chg [2nd]hide [clear]periodic", 0, SCR_H - 10);
        } else {
            safe_print("[alpha]chg [2nd]show [clear]periodic", 0, SCR_H - 10);
        }
    }

    /*
     * Keep the VSEPR panel on the topmost layer when enabled so it does not
     * get visually cut by any Lewis-structure render pass.
     */
    if (card_drawn && vsepr_card_enabled) {
        draw_vsepr_info_card(&mol, ls, ax, ay, vsepr_force_visible);
    }

    return card_drawn;
}

int main(void)
{
    gfx_Begin();
    gfx_SetDrawBuffer();

    init_pt_grid();

    /* Initialize cursor to Carbon (period 2, group 14 -> row 1, col 13) */
    cur_row = 1;
    cur_col = 13;

    molecule_reset(&mol);

    bool running = true;
    bool show_lewis = false;
    bool vsepr_force_visible = false;
    bool vsepr_card_enabled = true;
    bool last_card_drawn = false;
    bool warning = false;
    uint8_t warning_timer = 0;
    uint8_t key_delay = 0;

    timer_Control = TIMER1_ENABLE | TIMER1_32K | TIMER1_0INT | TIMER1_DOWN;
    timer_1_ReloadValue = FRAME_TICKS;
    timer_1_Counter = FRAME_TICKS;

    while (running) {
        while (timer_1_Counter > 0) {}
        timer_1_Counter = FRAME_TICKS;

        kb_Scan();

        if (show_lewis) {
            if (kb_Data[6] & kb_Clear) {
                show_lewis = false;
                vsepr_force_visible = false;
                vsepr_card_enabled = true;
                last_card_drawn = false;
                continue;
            }

            if (key_delay == 0 && (kb_Data[1] & kb_2nd)) {
                if (!vsepr_card_enabled) {
                    vsepr_card_enabled = true;
                    vsepr_force_visible = false;
                } else if (!last_card_drawn && !vsepr_force_visible) {
                    vsepr_force_visible = true;
                } else {
                    vsepr_card_enabled = false;
                    vsepr_force_visible = false;
                }
                key_delay = 8;
            }

            if (key_delay == 0 && (kb_Data[2] & kb_Alpha)) {
                cycle_charge(&mol);
                generate_resonance(&mol);
                key_delay = 8;
            }

            if (mol.num_res > 1) {
                if ((kb_Data[7] & kb_Right) && key_delay == 0) {
                    mol.cur_res = (mol.cur_res + 1) % mol.num_res;
                    key_delay = 8;
                }
                if ((kb_Data[7] & kb_Left) && key_delay == 0) {
                    mol.cur_res = (mol.cur_res == 0) ? mol.num_res - 1 : mol.cur_res - 1;
                    key_delay = 8;
                }
            }
            if (key_delay > 0) key_delay--;

            last_card_drawn = draw_lewis(vsepr_force_visible, vsepr_card_enabled);
        } else {
            if (kb_Data[1] & kb_Mode) {
                running = false;
                continue;
            }

            if (key_delay == 0) {
                if (kb_Data[7] & kb_Up)    { move_cursor(&cur_row, &cur_col, -1,  0); key_delay = 6; }
                if (kb_Data[7] & kb_Down)  { move_cursor(&cur_row, &cur_col,  1,  0); key_delay = 6; }
                if (kb_Data[7] & kb_Left)  { move_cursor(&cur_row, &cur_col,  0, -1); key_delay = 6; }
                if (kb_Data[7] & kb_Right) { move_cursor(&cur_row, &cur_col,  0,  1); key_delay = 6; }

                if (kb_Data[6] & kb_Enter) {
                    uint8_t ei = pt_grid[cur_row][cur_col];
                    if (ei != ELEM_NONE && mol.num_atoms < MAX_ATOMS) {
                        int heavy = 0;
                        for (uint8_t i = 0; i < mol.num_atoms; i++) {
                            if (mol.atoms[i].elem != ELEM_H) heavy++;
                        }

                        if (ei != ELEM_H && heavy >= MAX_HEAVY) {
                            warning = true;
                            warning_timer = 40;
                        } else {
                            mol.atoms[mol.num_atoms].elem = ei;
                            mol.num_atoms++;
                        }
                    }
                    key_delay = 8;
                }

                if (kb_Data[1] & kb_Del) {
                    if (mol.num_atoms > 0) {
                        mol.num_atoms--;
                        if (mol.num_atoms == 0) {
                            mol.charge = 0;
                            mol.invalid_reason = INVALID_NONE;
                        }
                    }
                    key_delay = 8;
                }

                if (kb_Data[2] & kb_Alpha) {
                    cycle_charge(&mol);
                    key_delay = 8;
                }

                if (kb_Data[1] & kb_2nd) {
                    if (mol.num_atoms >= 1) {
                        generate_resonance(&mol);
                        show_lewis = true;
                        vsepr_force_visible = false;
                        vsepr_card_enabled = true;
                        last_card_drawn = false;
                    }
                    key_delay = 10;
                }
            }
            if (key_delay > 0) key_delay--;

            draw_periodic_table(&mol, cur_row, cur_col);

            if (warning && warning_timer > 0) {
                gfx_SetColor(UI_ALERT_BG);
                gfx_FillRectangle(40, 100, 240, 30);
                gfx_SetColor(UI_ALERT_TEXT);
                gfx_Rectangle(40, 100, 240, 30);
                gfx_SetTextScale(1, 1);
                gfx_SetTextFGColor(UI_ALERT_TEXT);
                gfx_SetTextBGColor(UI_ALERT_BG);
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
