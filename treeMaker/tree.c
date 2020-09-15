#include "node.h"
#include "parse.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

void freeTree(struct node *tree);
void leafToLeafDistance(struct node *root, double **dist, char **name, char branchLength);
void subDist(struct node *current, double **dist, char **name, char branchLength, int size, int index);
void removeRoot(struct node **root);

void removeRoot(struct node **root) {
    struct node *children = (*root)->firstChild;
    struct node *lastSiblingOfChild = 0;
    while (children->numberOfChildren == 0) {
        children = children->nextSibling;
    }
    lastSiblingOfChild = children->firstChild;
    while (lastSiblingOfChild->nextSibling != 0) {
        lastSiblingOfChild = lastSiblingOfChild->nextSibling;
    }
    children = (*root)->firstChild;
    
    
    while (children != 0) {
        if (children != lastSiblingOfChild->parent) {
            lastSiblingOfChild->nextSibling = children;
            children->parent = lastSiblingOfChild->parent;
            lastSiblingOfChild = lastSiblingOfChild->nextSibling;
        } else {
            children = children->nextSibling;
        }
        if (children->nextSibling == lastSiblingOfChild->parent) {
            children->nextSibling = children->nextSibling->nextSibling;
        }
        children = children->nextSibling;
    }
    
    lastSiblingOfChild->parent->parent = 0;
    free(*root);
    (*root) = lastSiblingOfChild->parent;
    compNumberOfLeaves(*root);
}


void subDist(struct node *current, double **dist, char **name, char branchLength, int size, int index) {
    if (current->numberOfChildren == 0) {
        strcpy(name[index], current->name);
    } else {
        struct node *child = current->firstChild;
        for (int childNo = 0; childNo < current->numberOfChildren; childNo++) {
            if (child->parent->tag == 0) {
                for (int i = 0; i < index; i++) {
                    for (int j = index; j < index + child->numberOfLeaves; j++) {
                        dist[i][j - 1 - i] += (branchLength != 1) ? 1 : child->distToParent;
                    }
                }
                for (int i = index; i < index + child->numberOfLeaves; i++) {
                    for (int j = index + child->numberOfLeaves; j < size; j++) {
                        dist[i][j - 1 - i] += (branchLength != 1) ? 1 : child->distToParent;
                    }
                }
            } else {
                for (int i = index; i < index + child->numberOfLeaves; i++) {
                    if (child->nextSibling != 0) {
                        for (int j = index + child->numberOfLeaves; j < index + child->numberOfLeaves + child->nextSibling->numberOfLeaves; j++) {
                            dist[i][j - 1 - i] -= (DBL_MAX + 1);
                        }
                    }
                }
            }
            subDist(child, dist, name, branchLength, size, index);
            index = index + child->numberOfLeaves;
            child = child->nextSibling;
        }
    }
}

void leafToLeafDistance(struct node *root, double **dist, char **name, char branchLength) {
    // LCA is counted double therefore all distances are lowered by one
    if (branchLength != 1) {
        for (int i = 0; i < root->numberOfLeaves - 1; i++) {
            for (int j = 0; j < root->numberOfLeaves - 1 - i; j++) {
                dist[i][j] = -1;
            }
        }
    }
    subDist(root, dist, name, branchLength, root->numberOfLeaves, 0);
}

void freeTree(struct node *tree) {
    
    struct node *before = tree;
    struct node *current = tree;
    int direction = 0;
    while (before != NULL) {
        
        // go depth first
        while (current->firstChild != NULL) {
            direction = 1;
            before = current;
            current = current->firstChild;
        }
        
        // go to next sibling and restart while loop
        if (current->nextSibling != NULL) {
            direction = 2;
            before = current;
            current = current->nextSibling;
            continue;
        }
        
        // free current node
        free(current);
        
        // delete pointer to node
        if (direction == 1) {
            before->firstChild = NULL;
        } else if (direction == 2) {
            before->nextSibling = NULL;
        } else {
            before = NULL;
        }
        
        // restart from the top
        current = tree;
        direction = 0;
    }
}
