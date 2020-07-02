#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"
#include "parse.h"
#include <limits.h>

static void calcQ(double**q, double **distance, int size);
static void findMin(double **q, int size, int *i, int *j);
static void join(int i, int j, struct node *root, int size, char name, double **distance);
static void calcD(double **distance, double **newDistance, int size, int i, int j);
static void calcD2(double **distance, double **newDistance, int size, int i, int j, struct node *root);
static void calcD3(double **distance, double **newDistance, int size, int i, int j, struct node *root, int *weights);
void makeTreeFromDistanceArray(double **distance, int size, struct node **root, char **names);
void makeTree2(double **distance, int size, struct node **root, char **names);
void makeTree3(double **distance, int size, struct node **root, char **names);


//UPGMA TODO
void makeTree3(double **distance, int size, struct node **root, char **names) {
    if (size == 0) {
        return;
    }
    
    int *weights = (int *) calloc(sizeof(int), size);
    for (int i = 0; i < size; i++) {
        weights[i] = 1;
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
        findMin(distance, currentSize, i ,j);
        join(*i, *j, *root, currentSize, name, distance);
        currentSize--;
        calcD3(distance, newDistance,currentSize,*i,*j, (*root), weights);
        double **tmp = distance;
        distance = newDistance;
        newDistance = tmp;
    }
    
    (*root)->firstChild->distToParent = distance[0][1];
    (*root)->firstChild->parent = (*root);
    (*root)->firstChild->nextSibling->parent = (*root);
    (*root)->name[0] = placeholderName;
    
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


// WPGMA
void makeTree2(double **distance, int size, struct node **root, char **names) {
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
        findMin(distance, currentSize, i ,j);
        join(*i, *j, *root, currentSize, name, distance);
        currentSize--;
        calcD2(distance, newDistance,currentSize,*i,*j, (*root));
        double **tmp = distance;
        distance = newDistance;
        newDistance = tmp;
    }
    
    (*root)->firstChild->distToParent = distance[0][1];
    (*root)->firstChild->parent = (*root);
    (*root)->firstChild->nextSibling->parent = (*root);
    (*root)->name[0] = placeholderName;
    
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

static void calcQ(double**q, double **distance, int size) {
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

static void findMin(double **q, int size, int *i, int *j) {
    double min = (double) INT_MAX;
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

static void join(int i, int j, struct node *root, int size, char name, double **distance) {
//    printf("join %d\t%d\n", i, j);
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

static void calcD3(double **distance, double **newDistance, int size, int i, int j, struct node *root, int *weights) {
    for (int ii = 0; ii < size; ii++) {
        for (int jj = 0; jj < size; jj++) {
            int adjustI = -1 + (ii > i) + (ii >= j);
            int adjustJ = -1 + (jj > i) + (jj >= j);
            if (ii == jj) {
                newDistance[ii][jj] = 0;
            } else if (ii == 0) {
//                if (distance[i][jj + adjustJ] > distance[j][jj + adjustJ]) {
//                    newDistance[ii][jj] = distance[i][jj + adjustJ];
//                } else {
//                    newDistance[ii][jj] = distance[j][jj + adjustJ];
//                }
                newDistance[ii][jj] = (distance[i][jj + adjustJ] * weights[i] + distance[j][jj + adjustJ] * weights[j]) / (weights[i] + weights[j]);
            } else if (jj == 0) {
//                if (distance[i][ii + adjustI] > distance[j][ii + adjustI]) {
//                    newDistance[ii][jj] = distance[i][ii + adjustI];
//                } else {
//                    newDistance[ii][jj] = distance[j][ii + adjustI];
//                }
                newDistance[ii][jj] = (distance[i][ii + adjustI] * weights[i] + distance[j][ii + adjustI] * weights[j]) / (weights[i] + weights[j]);
            } else {
                newDistance[ii][jj] = distance[ii + adjustI][jj + adjustJ];
            }
        }
    }
    int *tmp = calloc(sizeof(int), size);
    tmp[0] = weights[i] + weights[j];
    for (int ii = 1; ii < size; ii++) {
        int adjustI = -1 + (ii > i) + (ii >= j);
        tmp[ii] = weights[i + adjustI];
    }
    for (int ii = 0; ii < size; ii++) {
        weights[ii] = tmp[ii];
    }
    free(tmp);
//    printf("\t");
//    struct node *current = root->firstChild;
//    for (int i = 0; i < size; i++) {
//        printf("%s\t", current->name);
//        current = current->nextSibling;
//    }
//    printf("\n\n");
//
//    for (int i = 0; i < size; i++) {
//        printf("\t");
//        double min = 3020203;
//        for (int j = 0; j < size; j++) {
//            if (newDistance[i][j] < min && newDistance[i][j] != 0) {
//                min = newDistance[i][j];
//            }
//        }
//        for (int j = 0; j < size; j++) {
//            if (newDistance[i][j] == min) {
//                printf("\033[1m\033[31m");
//            }
//            printf("%.2f\t", newDistance[i][j]);
//            printf("\033[0m");
//        }
//        printf("\n");
//    }
//    printf("\n");
}

static void calcD2(double **distance, double **newDistance, int size, int i, int j, struct node *root) {
    for (int ii = 0; ii < size; ii++) {
        for (int jj = 0; jj < size; jj++) {
            int adjustI = -1 + (ii > i) + (ii >= j);
            int adjustJ = -1 + (jj > i) + (jj >= j);
            if (ii == jj) {
                newDistance[ii][jj] = 0;
            } else if (ii == 0) {
//                if (distance[i][jj + adjustJ] > distance[j][jj + adjustJ]) {
//                    newDistance[ii][jj] = distance[i][jj + adjustJ];
//                } else {
//                    newDistance[ii][jj] = distance[j][jj + adjustJ];
//                }
                newDistance[ii][jj] = (distance[i][jj + adjustJ] + distance[j][jj + adjustJ]) / 2.0;
            } else if (jj == 0) {
//                if (distance[i][ii + adjustI] > distance[j][ii + adjustI]) {
//                    newDistance[ii][jj] = distance[i][ii + adjustI];
//                } else {
//                    newDistance[ii][jj] = distance[j][ii + adjustI];
//                }
                newDistance[ii][jj] = (distance[i][ii + adjustI] + distance[j][ii + adjustI]) / 2.0;
            } else {
                newDistance[ii][jj] = distance[ii + adjustI][jj + adjustJ];
            }
        }
    }
//    printf("\t");
//    struct node *current = root->firstChild;
//    for (int i = 0; i < size; i++) {
//        printf("%s\t", current->name);
//        current = current->nextSibling;
//    }
//    printf("\n\n");
//    
//    for (int i = 0; i < size; i++) {
//        printf("\t");
//        double min = 3020203;
//        for (int j = 0; j < size; j++) {
//            if (newDistance[i][j] < min && newDistance[i][j] != 0) {
//                min = newDistance[i][j];
//            }
//        }
//        for (int j = 0; j < size; j++) {
//            if (newDistance[i][j] == min) {
//                printf("\033[1m\033[31m");
//            }
//            printf("%.2f\t", newDistance[i][j]);
//            printf("\033[0m");
//        }
//        printf("\n");
//    }
//    printf("\n");
}


static void calcD(double **distance, double **newDistance, int size, int i, int j) {
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


void makeTreeFromDistanceArray(double **distance, int size, struct node **root, char **names) {
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
        findMin(q, currentSize, i ,j);
        join(*i, *j, *root, currentSize, name, distance);
        currentSize--;
        calcD(distance, newDistance,currentSize,*i,*j);
        double **tmp = distance;
        distance = newDistance;
        newDistance = tmp;
    }
    
    (*root)->firstChild->distToParent = distance[0][1];
    (*root)->firstChild->parent = (*root);
    (*root)->firstChild->nextSibling->parent = (*root);
    (*root)->name[0] = placeholderName;
    
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
