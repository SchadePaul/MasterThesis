#ifndef PARSE_H
#define PARSE_H

void readFileToArray(const char *filename, char ***newickTree, int *numberOfTrees);
void newickTreeToTree(char *newickTree, struct node **tree, char ***allNames, int *numberOfNames);
void printTree(struct node *tree);
void saveTree(struct node *tree, const char *name);

#endif
