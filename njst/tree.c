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

void editDistanceForUp(int index, struct node *current, double **dist, int size) {
    // Edit distances for going up
    for (int i = 0; i < index - current->numberOfLeaves + 1; i++) {
        for (int j = index + 1; j < size; j++) {
            dist[i][j - i] -= current->distToParent;
        }
    }
    for (int i = index - current->numberOfLeaves + 1; i <= index; i++) {
        for (int j = index + 1; j < size; j++) {
            dist[i][j - i] += current->distToParent;
        }
    }
}
void editDistanceForDown(int index, struct node *current, double **dist, int size) {
    // Edit distances for going down
    for (int i = 0; i < index; i++) {
        for (int j = index; j < size; j++) {
            dist[i][j - i] += current->distToParent;
        }
    }
}



// TODO: Last node might be not done right, only works if no internal nodes appear
void leafToLeafDistance(struct node *root, double **dist, int size, char **name) {
    struct node *current = root;
    int index = 0;
    // Iterate through all terminal nodes ahead
    while (index < size - 1) {
        // Go depth first
        while (current->firstChild != NULL) {
            current = current->firstChild;
            editDistanceForDown(index, current, dist, size);
        }
        // Get name of terminal node
        strcpy(name[index], current->name);
        
        // Go up until next sibling is found
        while (current->nextSibling == NULL) {
            editDistanceForUp(index, current, dist, size);
            current = current->parent;
        }
        
        // Go to next sibling
        editDistanceForUp(index, current, dist, size);
        current = current->nextSibling;
        index++;
        
        // Get name of last terminal node
        if (index == size - 1) {
            strcpy(name[index], current->name);
        }
        editDistanceForDown(index, current, dist, size);
    }
    
}
