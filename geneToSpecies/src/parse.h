#ifndef PARSE_H
#define PARSE_H

void readFileToTrees(const char *filename, struct node ***trees_ptr, char ***taxa_ptr, int *numberOfTrees, int *numberOfTaxa);
void printTree(struct node *tree, int length);
void saveTree(struct node *tree, const char *name);
int compNumberOfLeaves(struct node *current);

#endif
