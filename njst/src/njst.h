#include "node.h"

#ifndef NJST_H
#define NJST_H

void taggedNJFromFile(struct node **root, const char *filename, int norm, int rooted, int weighted, int square);
void miniNJFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength);
void njstFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength, int square);
void ustar(struct node **root, const char *filename);

#endif
