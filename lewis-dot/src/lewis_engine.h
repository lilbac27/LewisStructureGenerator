#ifndef LEWIS_ENGINE_H
#define LEWIS_ENGINE_H

#include "lewis_model.h"

void generate_resonance(Molecule *mol);
const char *invalid_reason_message(InvalidReason reason);

#endif
