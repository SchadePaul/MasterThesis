#include "node.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void freeTree(struct node *tree) {
    struct node *before = tree;
    struct node *current = tree;
    int didDelete = 1;
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

void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size, double norm) {
    for (int i = index - numberOfLeaves + 1; i <= index; i++) {
        for (int j = index + 1; j < size; j++) {
            allDist[i][j - i] += dist / norm ;
        }
    }
}

void goingDown(int index, int numberOfLeaves, double dist, double **allDist, double norm) {
    for (int i = 0; i < index; i++) {
        for (int j = index; j < index + numberOfLeaves; j++) {
            allDist[i][j - i] += dist / norm;
        }
    }
}


void leafToLeafDistance(struct node *root, double **dist, int size, char **name, int normDistance, int branchLength) {
    int index = 0;
    struct node *current = root;
    double norm = 1;
    if (normDistance) {
        norm = (double) size;
    }
    
    // go to first leave
    while (current->firstChild != NULL) {
        current = current->firstChild;
    }
    
    while (1) {
        while (current->firstChild != NULL) {
            current = current->firstChild;
            if (!branchLength) {
                goingDown(index, current->numberOfLeaves, 1.0, dist, norm);
            } else {
                goingDown(index, current->numberOfLeaves, current->distToParent, dist, norm);
            }
        }
        
        strcpy(name[index], current->name);
        
        while (current->nextSibling == NULL) {
            if (index < size - 1) {
                if (!branchLength) {
                    goingUp(index, current->numberOfLeaves, 1.0, dist, size, norm);
                } else {
                    goingUp(index, current->numberOfLeaves, current->distToParent, dist, size, norm);
                }
            }
            
            struct node *toFree = current;
            current = current->parent;
            free(toFree);
            if (current == root) {
                break;
            }
        }
        if (current == root) {
            free(current);
            break;
        }
        if (!branchLength) {
            goingUp(index, current->numberOfLeaves, 1.0, dist, size, norm);
        } else {
            goingUp(index, current->numberOfLeaves, current->distToParent, dist, size, norm);
        }
        index++;
        struct node *toFree = current;
        current = current->nextSibling;
        free(toFree);
        if (!branchLength) {
            goingDown(index, current->numberOfLeaves, 1.0, dist, norm);
        } else {
            goingDown(index, current->numberOfLeaves, current->distToParent, dist, norm);
        }
    }
}
