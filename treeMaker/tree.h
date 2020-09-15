#ifndef TREE_H
#define TREE_H

void removeRoot(struct node **root);
void freeTree(struct node *tree);
void leafToLeafDistance(struct node *root, double **dist, char **name, char branchLength);

#endif

