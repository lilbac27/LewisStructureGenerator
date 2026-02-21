#ifndef VSEPR_DATA_H
#define VSEPR_DATA_H

typedef struct {
    const char *valence_pairs;
    const char *ep_geometry;
    const char *bond_pairs;
    const char *lone_pairs;
    const char *shape;
    const char *hybridization;
} vsepr_row_t;

#define VSEPR_NUM_ROWS 13

const vsepr_row_t vsepr_data[VSEPR_NUM_ROWS] = {
    {"2", "Linear", "2", "0", "Linear", "sp"},
    {"3", "Trigonal Planar", "3", "0", "Trigonal Planar", "sp2"},
    {"3", "Trigonal Planar", "2", "1", "Bent", "sp2"},
    {"4", "Tetrahedral", "4", "0", "Tetrahedral", "sp3"},
    {"4", "Tetrahedral", "3", "1", "Trigonal Pyramidal", "sp3"},
    {"4", "Tetrahedral", "2", "2", "Bent", "sp3"},
    {"5", "Trigonal Bipyramidal", "5", "0", "Trigonal Bipyramidal", "sp3d"},
    {"5", "Trigonal Bipyramidal", "4", "1", "Seesaw", "sp3d"},
    {"5", "Trigonal Bipyramidal", "3", "2", "T-shaped", "sp3d"},
    {"5", "Trigonal Bipyramidal", "2", "3", "Linear", "sp3d"},
    {"6", "Octahedral", "6", "0", "Octahedral", "sp3d2"},
    {"6", "Octahedral", "5", "1", "Square Pyramidal", "sp3d2"},
    {"6", "Octahedral", "4", "2", "Square Planar", "sp3d2"},
};

const char *vsepr_headers[6] = {
    "valence electron pairs",
    "electron pair geometry",
    "bond pairs",
    "lone pairs",
    "shape",
    "hybridization",
};

#endif
