#ifndef TREE_H
#define TREE_H

void freeTree(struct node *tree);
void leafToLeafDistance(struct node *root, double **dist, char **name, int branchLength);

#endif
