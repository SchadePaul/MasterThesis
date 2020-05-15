#ifndef NODE_H
#define NODE_H

static int const maxNameLength = 40;
static int const maxValueLength = 50;
static char const placeholderName = '*';

struct node {
    char name[40];
    int numberOfLeaves;
    double distToParent;
    struct node *parent;
    struct node *firstChild;
    struct node *nextSibling;
    int tag;    // tag == 0: Speciation, tag == 1: Duplication
    int idNo;   // given at first read to identify while rooting
    int score;
};

#endif
