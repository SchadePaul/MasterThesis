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

void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size) {
    for (int i = index - numberOfLeaves + 1; i <= index; i++) {
        for (int j = index + 1; j < size; j++) {
            allDist[i][j - i] += 1;
        }
    }
}

void goingDown(int index, int numberOfLeaves, double dist, double **allDist) {
    for (int i = 0; i < index; i++) {
        for (int j = index; j < index + numberOfLeaves; j++) {
            allDist[i][j - i] += 1;
        }
    }
}


void leafToLeafDistance(struct node *root, double **dist, int size, char **name) {
    int index = 0;
    struct node *current = root;
    
    // go to first leave
    while (current->firstChild != NULL) {
        current = current->firstChild;
    }
    
    while (1) {
        while (current->firstChild != NULL) {
            current = current->firstChild;
            goingDown(index, current->numberOfLeaves, current->distToParent, dist);
        }
        
        strcpy(name[index], current->name);
        
        while (current->nextSibling == NULL) {
            if (index < size - 1) {
                goingUp(index, current->numberOfLeaves, current->distToParent, dist, size);
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
        
        goingUp(index, current->numberOfLeaves, current->distToParent, dist, size);
        index++;
        struct node *toFree = current;
        current = current->nextSibling;
        free(toFree);
        goingDown(index, current->numberOfLeaves, current->distToParent, dist);
        
    }
}
