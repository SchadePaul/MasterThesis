#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"
#include "parse.h"
#include <limits.h>

static void calcQ(double**q, double **distance, int size);
static void findMinOfQ(double **q, int size, int *i, int *j);
static void join(int i, int j, struct node *root, int size, char name, double **distance);
static void calcD(double **distance, double **newDistance, int size, int i, int j);
void makeTreeFromDistanceArray(double **distance, int size, struct node **root, char **names);
void makeTree2(double **distance, int size, struct node **root, char **names);
static void findMin(double **distance, int size, int *index_1, int *index_2);
static void calcDistance(double **oldDistance, double **newDistance, int size, int index_1, int index_2);
static void joiny(struct node *root, int size, int index_1, int index_2);

static void joiny(struct node *root, int size, int index_1, int index_2) {
    printf("join at size\t%d\t%d\t%d\n", size, index_1, index_2);
    struct node *new = (struct node *) calloc(sizeof(struct node), 1);
    new->name[0] = placeholderName;
    new->distToParent = 1;
    struct node *a;
    struct node *b;

    struct node *current = root->firstChild;

    for (int i = 0; i < size; i++) {
        struct node *next = current->nextSibling;
        if (i == 0 && index_1 == 0) {
            a = current;
            if (index_2 == 1) {
                b = next;
                root->firstChild = new;
                new->nextSibling = b->nextSibling;
                break;
            } else {
                root->firstChild = new;
                new->nextSibling = next;
                current = next;
                continue;
            }
        }
        
        if (i == index_1 - 1) {
            a = next;
            next = next->nextSibling;
            current->nextSibling = new;
            current = new;
            i++;
            if (i == index_2 - 1) {
                b = next;
                next = next->nextSibling;
                i++;
            }
        } else if (i == index_2 - 1){
            b = next;
            next = next->nextSibling;
            i++;
        }
        current->nextSibling = next;
        current = next;
    }
    
    new->firstChild = a;
    a->nextSibling = b;
    b->nextSibling = 0;
    a->parent = new;
    b->parent = new;
    
    
}

static void calcDistance(double **oldDistance, double **newDistance, int size, int index_1, int index_2) {
    for (int i = 0; i < size; i++) {
        int ii = i;
        if (i == index_2) {
            ii = index_1;
        } else if (i > index_2) {
            ii--;
        }
        for (int j = 0; j < size; j++) {
            int jj = j;
            if (j == index_2) {
                jj = index_1;
            } else if (j > index_2) {
                jj--;
            }
            if (newDistance[ii][jj] == 0 || oldDistance[i][j] < newDistance[ii][jj]) {
                newDistance[ii][jj] = oldDistance[i][j];
            }
            oldDistance[i][j] = 0;
        }
    }
}

static void findMin(double **distance, int size, int *index_1, int *index_2) {
    int min = INT_MAX;
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (distance[i][j] < min) {
                min = distance[i][j];
                *index_1 = i;
                *index_2 = j;
            }
        }
    }
}

void makeTree2(double **distance, int size, struct node **root, char **names) {
    if (size == 0) {
        return;
    } else {
        printf("size: %d\n", size);
    }
    double **distance_1 = (double **) calloc(sizeof(double *), size);
    double **distance_2 = (double **) calloc(sizeof(double *), size);
    for (int i = 0; i < size; i++) {
        distance_1[i] = (double *) calloc(sizeof(double), size);
        distance_2[i] = (double *) calloc(sizeof(double), size);
        for (int j = 0; j < size; j++) {
            distance_1[i][j] = distance[i][j];
        }
    }
    struct node *current = (struct node *) calloc(sizeof(struct node), 1);
    strcpy(current->name, names[0]);
    (*root)->name[0] = placeholderName;
    (*root)->firstChild = current;
    current->parent = *root;
    current->distToParent = 1;
    for (int i = 1; i < size; i++) {
        struct node *next = (struct node *) calloc(sizeof(struct node), 1);
        strcpy(next->name, names[i]);
        current->nextSibling = next;
        next->parent = *root;
        current = next;
        current->distToParent = 1;
    }
    
    for (int i = 0; i < size - 2; i++) {
        int index_1;
        int index_2;
        findMin(distance_1, size - i, &index_1, &index_2);
        calcDistance(distance_1, distance_2, size - i, index_1, index_2);
        double **tmp = distance_1;
        distance_1 = distance_2;
        distance_2 = tmp;
        joiny(*root, size - i, index_1, index_2);
//        compNumberOfLeaves(*root);
//        printTree(*root, 4);
    }
    
    compNumberOfLeaves(*root);
    
    for (int i = 0; i < size; i++) {
        free(distance_1[i]);
        free(distance_2[i]);
    }
    free(distance_1);
    free(distance_2);
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

static void findMinOfQ(double **q, int size, int *i, int *j) {
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
        findMinOfQ(q, currentSize, i ,j);
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
