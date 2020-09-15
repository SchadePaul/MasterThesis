#include "node.h"

#ifndef PARSE_H
#define PARSE_H

void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames);
void readFileToArray(const char *filename, char ***newickTree, int *numberOfTrees);
void newickTreeToTree(char *newickTree, struct node **tree, char ***allLeafNames, int *numberOfLeafNames);
void printTree(struct node *tree, int length);
void saveTree(struct node *tree, const char *name);
int compNumberOfLeaves(struct node *current);

#endif

