#ifndef NODE_H
#define NODE_H

int static const maxNameLength = 5;
int static const maxValueLength = 50;

struct node {
    char name[maxNameLength];
    int numberOfTerminalNodes;
    double distToParent;
    struct node *parent;
    struct node *firstChild;
    struct node *nextSibling;
};

#endif
