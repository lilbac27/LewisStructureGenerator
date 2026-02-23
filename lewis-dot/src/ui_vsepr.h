#ifndef UI_VSEPR_H
#define UI_VSEPR_H

#include <stdbool.h>

#include "lewis_model.h"

bool draw_vsepr_info_card(const Molecule *mol, const LewisStructure *ls, const int ax[MAX_ATOMS], const int ay[MAX_ATOMS], bool force_visible);

#endif
