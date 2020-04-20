#include "node.h"

#ifndef NJST_H
#define NJST_H

void njstFromFile(struct node **root, const char *filename, int branchLength, int minNJst, int normDistance, double quartil);

#endif
