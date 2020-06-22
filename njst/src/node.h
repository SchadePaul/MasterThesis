#ifndef NODE_H
#define NODE_H

static int const maxNameLength = 4;
static int const maxValueLength = 50;
static char const placeholderName = '*';

struct node {
    char name[4];
    
    int numberOfLeaves;
    int numberOfChildren;
    int numberOfNodes;
    
    double distToParent;
    
    struct node *parent;
    struct node *firstChild;
    struct node *nextSibling;
    
    int tag;
    int *tag2;    // tag == 0: Speciation, tag == 1: Duplication
    int idNo;   // given at first read to identify while rooting
    int score;
    
};

#endif
