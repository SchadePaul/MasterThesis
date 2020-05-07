#ifndef NODE_H
#define NODE_H

int static const maxNameLength = 4;
int static const maxValueLength = 50;
char static const placeholderName = '*';

struct node {
    char name[4];
    int numberOfLeaves;
    double distToParent;
    struct node *parent;
    struct node *firstChild;
    struct node *nextSibling;
    int tag;    // tag == 0: Speciation, tag == 1: Duplication
    int idNo;
    int score;
};

#endif
