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
static void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size);
static void goingDown(int index, int numberOfLeaves, double dist, double **allDist);
void leafToLeafDistance(struct node *root, double **dist, char **name, int branchLength);

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


// Distance functions

static void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size) {
    for (int i = index - numberOfLeaves + 1; i <= index; i++) {
        for (int j = index + 1; j < size; j++) {
            allDist[i][j - i] += dist;
        }
    }
}

static void goingDown(int index, int numberOfLeaves, double dist, double **allDist) {
    for (int i = 0; i < index; i++) {
        for (int j = index; j < index + numberOfLeaves; j++) {
            allDist[i][j - i] += dist;
        }
    }
}

void leafToLeafDistance(struct node *root, double **dist, char **name, int branchLength) {
    int size = root->numberOfLeaves;
    int index = 0;
    struct node *current = root;
    
    // go to first leave
    while (current->firstChild != NULL) {
        current = current->firstChild;
    }
    
    while (1) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            if (!branchLength) {
                goingDown(index, current->numberOfLeaves, 1.0, dist);
            } else {
                goingDown(index, current->numberOfLeaves, current->distToParent, dist);
            }
        }
        
        strcpy(name[index], current->name);
        
        while (current->nextSibling == 0) {
            if (index < size - 1) {
                if (!branchLength) {
                    goingUp(index, current->numberOfLeaves, 1.0, dist, size);
                } else {
                    goingUp(index, current->numberOfLeaves, current->distToParent, dist, size);
                }
            }
            
            current = current->parent;
            if (current == root) {
                break;
            }
        }
        if (current == root) {
            break;
        }
        if (!branchLength) {
            goingUp(index, current->numberOfLeaves, 1.0, dist, size);
        } else {
            goingUp(index, current->numberOfLeaves, current->distToParent, dist, size);
        }
        
        index++;
        current = current->nextSibling;
        if (!branchLength) {
            goingDown(index, current->numberOfLeaves, 1.0, dist);
        } else {
            goingDown(index, current->numberOfLeaves, current->distToParent, dist);
        }
    }
}
