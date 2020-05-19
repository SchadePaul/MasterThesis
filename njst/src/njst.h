#include "node.h"

#ifndef NJST_H
#define NJST_H

void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root, int miniPairs);

#endif
