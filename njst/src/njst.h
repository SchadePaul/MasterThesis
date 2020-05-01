#include "node.h"

#ifndef NJST_H
#define NJST_H

void taggedNJFromFile(struct node **root, const char *filename, int norm, int rooted);
void miniNJFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength);
void njstFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength);

#endif
