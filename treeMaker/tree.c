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
void leafToLeafDistance(struct node *root, double **dist, char **name, char branchLength, char astralTag, char notCountTag);
void subDist(struct node *current, double **dist, char **name, char branchLength, int size, int index, char astralTag, char notCountTag);
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


void subDist(struct node *current, double **dist, char **name, char branchLength, int size, int index, char astralTag, char notCountTag) {
    if (current->numberOfChildren == 0) {
        strcpy(name[index], current->name);
    } else {
        struct node *child = current->firstChild;
        for (int childNo = 0; childNo < current->numberOfChildren; childNo++) {
            if (astralTag == 0 || (astralTag == 1 && child->parent->tag == 0)) {
                for (int i = 0; i < index; i++) {
                    for (int j = index; j < index + child->numberOfLeaves; j++) {
                        dist[i][j - 1 - i] += (branchLength != 1) ? 1 : child->distToParent;
                        if (notCountTag == 1 && astralTag == 1 && child->tag == 1) {
                            dist[i][j - 1 - i] -= (branchLength != 1) ? 1 : child->distToParent;
                        }
                    }
                }
                for (int i = index; i < index + child->numberOfLeaves; i++) {
                    for (int j = index + child->numberOfLeaves; j < size; j++) {
                        dist[i][j - 1 - i] += (branchLength != 1) ? 1 : child->distToParent;
                        if (notCountTag == 1 && astralTag == 1 && child->tag == 1) {
                            dist[i][j - 1 - i] -= (branchLength != 1) ? 1 : child->distToParent;
                        }
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
            subDist(child, dist, name, branchLength, size, index, astralTag, notCountTag);
            index = index + child->numberOfLeaves;
            child = child->nextSibling;
        }
    }
}

void leafToLeafDistance(struct node *root, double **dist, char **name, char branchLength, char astralTag, char notCountTag) {
    // LCA is counted double therefore all distances are lowered by one
    if (branchLength != 1) {
        for (int i = 0; i < root->numberOfLeaves - 1; i++) {
            for (int j = 0; j < root->numberOfLeaves - 1 - i; j++) {
                dist[i][j] = -1;
            }
        }
    }
    subDist(root, dist, name, branchLength, root->numberOfLeaves, 0, astralTag, notCountTag);
}

void freeTree(struct node *tree) {
    if (tree->firstChild != 0) {
        freeTree(tree->firstChild);
    }
    if (tree->nextSibling != 0) {
        freeTree(tree->nextSibling);
    }
    free(tree);
}
