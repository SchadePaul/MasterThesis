#ifndef NODE_H
#define NODE_H

static int const maxNameLength = 40;
static int const maxValueLength = 20;
static char const placeholderName = '*';

struct node {
    char name[40];
    
    int numberOfLeaves;
    int numberOfChildren;
    int numberOfNodes;
    
    double distToParent;
    
    struct node *parent;
    struct node *firstChild;
    struct node *nextSibling;
    
    /*
        tag == 0: Speciation,
        tag == 1: Duplication
    */
    int tag;
    
    int idNo;   // given at first read to identify while rooting
    int score;
    
};

#endif
