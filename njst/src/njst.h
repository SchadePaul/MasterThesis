#include "node.h"

#ifndef NJST_H
#define NJST_H

void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames);

void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root, double miniPairs, double quartil, int closeFriends);
void test(struct node **speciesTree, const char *filename);

#endif
