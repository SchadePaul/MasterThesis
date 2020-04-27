#include "node.h"

#ifndef NJST_H
#define NJST_H

void njstFromFile(struct node **root, const char *filename, int norm, int weighted, int miniNJ, int branchLength);

#endif
