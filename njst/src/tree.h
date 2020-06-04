#ifndef TREE_H
#define TREE_H

void leafToLeafDistance(struct node *root, double **dist, char **name, int branchLength);
void freeTree(struct node *tree);
void tagAndRoot(struct node *tree);
int scoreAndTag(struct node *tree, char **names);
int scoreAndTag2(struct node *tree);
void deletedTaggedDistance(struct node *current, double **dist, int *index);

#endif
