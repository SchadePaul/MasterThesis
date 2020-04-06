#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"

void calcQ(double**q, double **distance, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            q[i][j] = (double) (size - 2) * distance[i][j];
            for (int k = 0; k < size; k++) {
                if (k != i) {
                    q[i][j] -= distance[i][k];
                }
                if (k != j) {
                    q[i][j] -= distance[j][k];
                }
            }
        }
    }
}

void findMinOfQ(double **q, int size, int *i, int *j) {
    //TODO find max int
    double min = 999999;
    for (int ii = 0; ii < size; ii++) {
        for (int jj = ii + 1; jj < size; jj++) {
            if(q[ii][jj] < min) {
                *i = ii;
                *j = jj;
                min = q[ii][jj];
            }
        }
    }
}

void join(int i, int j, struct node *root, int size, char
          name, double **distance) {
    struct node *newNode = (struct node*) calloc(sizeof(struct node), 1);
    struct node *current = root->firstChild;
    root->firstChild = newNode;
    root->firstChild->name[0] = name;
    for (int currentpos = 0; currentpos < size; currentpos++) {
        if (currentpos == i) {
            newNode->firstChild = current;
            current->parent = newNode;
            double dist = distance[i][j] / (double) 2;
            double sum = 0;
            for (int ii = 0; ii < size; ii++) {
                sum += (distance[i][ii] - distance[j][ii]);
            }
            dist += sum / (double) (2 * (size - 2));
            current->distToParent = dist;
            if (current->nextSibling != 0) {
                current = current->nextSibling;
            }
        } else if (currentpos == j) {
            newNode->firstChild->nextSibling = current;
            current->parent = newNode;
            double dist = distance[i][j] - newNode->firstChild->distToParent;
            current->distToParent = dist;
            if (current->nextSibling != 0) {
                current = current->nextSibling;
            }
            newNode->firstChild->nextSibling->nextSibling = 0;
        } else {
            struct node *nextSiblingInsert = root->firstChild;
            while (nextSiblingInsert->nextSibling != 0) {
                nextSiblingInsert = nextSiblingInsert->nextSibling;
            }
            nextSiblingInsert->nextSibling = current;
            if (current->nextSibling != 0) {
                current = current->nextSibling;
            }
            nextSiblingInsert->nextSibling->nextSibling = 0;
            
        }
    }
    
}

void calcD(double **distance, double **newDistance, int size, int i, int j) {
    for (int ii = 0; ii < size; ii++) {
        for (int jj = 0; jj < size; jj++) {
            int adjustI = -1 + (ii > i) + (ii >= j);
            int adjustJ = -1 + (jj > i) + (jj >= j);
            if (ii == jj) {
                newDistance[ii][jj] = 0;
            } else if (ii == 0) {
                newDistance[ii][jj] = (distance[i][jj + adjustJ] + distance[j][jj + adjustJ] - distance[i][j]) / 2.0;
            } else if (jj == 0) {
                newDistance[ii][jj] = (distance[i][ii + adjustI] + distance[j][ii + adjustI] - distance[i][j]) / 2.0;
            } else {
                newDistance[ii][jj] = distance[ii + adjustI][jj + adjustJ];
            }
        }
    }
}

void removeRoot(struct node **root, double dist) {
    (*root)->firstChild->nextSibling->distToParent = dist;
    (*root)->firstChild->nextSibling->parent = (*root)->firstChild;
    (*root)->firstChild->firstChild->nextSibling->nextSibling = (*root)->firstChild->nextSibling;
    (*root) = (*root)->firstChild;
    (*root)->parent = NULL;
    (*root)->nextSibling = NULL;
}

void makeTreeFromDistanceMatrix(double **distance, int size, struct node **root, char **names) {
    if (size == 0) {
        return;
    }
    char name = placeholderName;
    (*root)->firstChild = (struct node *) calloc(sizeof(struct node), 1);
    struct node *current = (*root)->firstChild;
    double **currentDistance = (double**) calloc(sizeof(double*), size);
    double **newDistance = (double**) calloc(sizeof(double*), size);
    double **q = (double**) calloc(sizeof(double*), size);
    for (int i = 0; i < size; i++) {
        currentDistance[i] = (double*) calloc(sizeof(double), size);
        newDistance[i] = (double*) calloc(sizeof(double), size);
        q[i] = (double*) calloc(sizeof(double), size);
        current->parent = *root;
        strcpy(current->name, names[i]);
        if (i < size - 1) {
            current->nextSibling = (struct node *) calloc(sizeof(struct node), 1);
            current = current->nextSibling;
        }
    }
    int *i = calloc(sizeof(int), 1);
    int *j = calloc(sizeof(int), 1);
    int currentSize = size;

    while (currentSize > 2) {
        calcQ(q, distance, currentSize);
        findMinOfQ(q, currentSize, i ,j);
        join(*i, *j, *root, currentSize, name, distance);
        currentSize--;
        calcD(distance, newDistance,currentSize,*i,*j);
        double **tmp = distance;
        distance = newDistance;
        newDistance = tmp;
    }
    
    if (size > 2) {
        removeRoot(root, distance[0][1]);
    }
    
    for (int i = 0; i < size; i++) {
        free(currentDistance[i]);
        free(distance[i]);
        free(q[i]);
    }
    free(currentDistance);
    free(distance);
    free(q);
    free(i);
    free(j);
}
