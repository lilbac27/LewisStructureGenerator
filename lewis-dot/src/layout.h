#ifndef LAYOUT_H
#define LAYOUT_H

#include <stdbool.h>

#include "lewis_model.h"

bool layout_linear_chain(const Molecule *mol, const LewisStructure *ls, int ax[MAX_ATOMS], int ay[MAX_ATOMS]);
bool layout_tree_from_central(const Molecule *mol, const LewisStructure *ls, int ax[MAX_ATOMS], int ay[MAX_ATOMS]);

#endif
