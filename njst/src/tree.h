#ifndef TREE_H
#define TREE_H

void leafToLeafDistance(struct node *root, double **dist, int size, char **name, int normDistance, int branchLength);
void freeTree(struct node *tree);

#endif
