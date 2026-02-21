#ifndef UI_PERIODIC_H
#define UI_PERIODIC_H

#include <stdint.h>

#include "lewis_model.h"

void move_cursor(uint8_t *cur_row, uint8_t *cur_col, int dr, int dc);
void draw_periodic_table(const Molecule *mol, uint8_t cur_row, uint8_t cur_col);

#endif
