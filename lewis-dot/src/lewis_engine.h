#ifndef LEWIS_ENGINE_H
#define LEWIS_ENGINE_H

#include "lewis_model.h"

typedef struct {
    uint8_t valence_pairs;
    uint8_t bond_pairs;
    uint8_t lone_pairs;
    const char *ep_geometry;
    const char *shape;
    const char *hybridization;
    const char *bond_angle;
} VseprInfo;

void generate_resonance(Molecule *mol);
bool lewis_get_vsepr_info(const Molecule *mol, const LewisStructure *ls, VseprInfo *out);
const char *invalid_reason_message(InvalidReason reason);

#endif
