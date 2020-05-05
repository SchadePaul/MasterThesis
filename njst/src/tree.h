#ifndef TREE_H
#define TREE_H

void leafToLeafDistance(struct node *root, double **dist, char **name, int normDistance, int branchLength);
void freeTree(struct node *tree);
void tagAndRoot(struct node **tree);
int scoreAndTag(struct node *tree);

#endif
