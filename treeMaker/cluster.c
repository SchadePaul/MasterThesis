#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include "node.h"
#include "parse.h"

void UWPGMA(struct node **root, char **names, double **dist, int size, char wpgma);
static void findMin(double **matrix, int size, int *i, int *j);
static void calcDistUWPGMA(double dist, struct node *nodeI, struct node *nodeJ, double *distI, double*distJ);
static void updatedDistWPGMA(double ***matrix, int joinI, int joinJ, int size);
static void updatedDistUPGMA(double ***matrix, struct node **nodes, int joinI, int joinJ, int size);
static void join(struct node ***nodes, int size, int i, int j, double distI, double distJ);
void NJ(struct node **root, char **names, double **dist, int size);
static void calcQ(double **dist, double **q, int size);
static void calcDistNJ(double **dist, int i, int j, int size, double *distI, double *distJ);
static void updatedDistNJ(double ***matrix, int joinI, int joinJ, int size);

void NJ(struct node **root, char **names, double **dist, int size) {
    double **clusterDist = (double **) calloc((size_t) size, sizeof(double *));
    for (int i = 0; i < size; i++) {
        clusterDist[i] = (double *) calloc((size_t) size, sizeof(double *));
        for (int j = 0; j < size; j++) {
            clusterDist[i][j] = dist[i][j];
        }
    }
    struct node **nodes = (struct node **) calloc((size_t) size, sizeof(struct node *));
    for (int i = 0; i < size; i++) {
        nodes[i] = (struct node *) calloc(1, sizeof(struct node));
        strcpy(nodes[i]->name, names[i]);
        nodes[i]->numberOfLeaves = 1;
    }
    while (size > 2) {
//        for (int i = 0; i < size; i++) {
//            for (int j = 0; j < size; j++) {
//                printf("%.3f\t", clusterDist[i][j]);
//            }
//            printf("\n");
//        }
//        printf("\n");
        double **q = (double **) calloc((size_t) size, sizeof(double *));
        for (int i = 0; i < size; i++) {
            q[i] = (double *) calloc((size_t) size, sizeof(double));
        }
        calcQ(clusterDist, q, size);
        int i = 0;
        int j = 0;
        findMin(q, size, &i, &j);
        double distI = 0.0;
        double distJ = 0.0;
        calcDistNJ(clusterDist, i, j, size, &distI, &distJ);
        updatedDistNJ(&clusterDist, i, j, size);
        join(&nodes, size, i, j, distI, distJ);
        for (int i = 0; i < size; i++) {
            free(q[i]);
        }
        free(q);
        size--;
        compNumberOfLeaves(nodes[0]);
        if (size == 2) {
            nodes[1]->parent = nodes[0];
            nodes[1]->distToParent = clusterDist[0][1];
            nodes[0]->firstChild->nextSibling->nextSibling = nodes[1];
        }
    }
    (*root) = nodes[0];
}

static void calcQ(double **dist, double **q, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double qij = (size - 2) * dist[i][j];
            for (int k = 0; k < size; k++) {
                qij -= dist[i][k];
                qij -= dist[j][k];
            }
            q[i][j] = qij;
        }
    }
}

static void calcDistNJ(double **dist, int i, int j, int size, double *distI, double *distJ) {
    *distI = (1 / 2.0) * dist[i][j];
    double add = 0.0;
    for (int k = 0; k < size; k++) {
        add += dist[i][k];
        add -= dist[j][k];
    }
    add *= (1 / (2 * (size - 2)));
    *distI += add;
    *distJ = dist[i][j] - *distI;
}

static void updatedDistNJ(double ***matrix, int joinI, int joinJ, int size) {
    double **newMatrix = (double **) calloc((size_t) size - 1, sizeof(double *));
    for (int i = 0; i < size - 1; i++) {
        newMatrix[i] = (double *) calloc((size_t) size - 1, sizeof(double *));
    }
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - 1; j++) {
            if (i == j) {
                newMatrix[i][j] = 0.0;
            } else if (i < j) {
                if (i == 0) {
                    int indJ = j - (j <= joinI) + (j >= joinJ);
                    newMatrix[0][j] = (1 / 2.0) * ((*matrix)[joinI][indJ] + (*matrix)[joinJ][indJ] - (*matrix)[joinI][joinJ]);
                } else {
                    int indI = i - (i <= joinI) + (i >= joinJ);
                    int indJ = j - (j <= joinI) + (j >= joinJ);
                    newMatrix[i][j] = (*matrix)[indI][indJ];
                }
                
            } else {
                newMatrix[i][j] = newMatrix[j][i];
            }
        }
    }
    for (int i = 0; i < size; i++) {
        free((*matrix)[i]);
    }
    free(*matrix);
    (*matrix) = newMatrix;
}

void UWPGMA(struct node **root, char **names, double **dist, int size, char wpgma) {
    double **clusterDist = (double **) calloc((size_t) size, sizeof(double *));
    for (int i = 0; i < size; i++) {
        clusterDist[i] = (double *) calloc((size_t) size, sizeof(double *));
        for (int j = 0; j < size; j++) {
            clusterDist[i][j] = dist[i][j];
        }
    }
    struct node **nodes = (struct node **) calloc((size_t) size, sizeof(struct node *));
    for (int i = 0; i < size; i++) {
        nodes[i] = (struct node *) calloc(1, sizeof(struct node));
        strcpy(nodes[i]->name, names[i]);
        nodes[i]->numberOfLeaves = 1;
    }
    while (size > 1) {
//        for (int i = 0; i < size; i++) {
//            for (int j = 0; j < size; j++) {
//                printf("%.3f\t", clusterDist[i][j]);
//            }
//            printf("\n");
//        }
//        printf("\n");
        int i = 0;
        int j = 0;
        findMin(clusterDist, size, &i, &j);
        double distI = 0.0;
        double distJ = 0.0;
        calcDistUWPGMA(clusterDist[i][j], nodes[i], nodes[j], &distI, &distJ);
        if (wpgma == 1) {
            updatedDistWPGMA(&clusterDist, i, j, size);
        } else {
            updatedDistUPGMA(&clusterDist, nodes, i ,j ,size);
        }
        join(&nodes, size, i, j, distI, distJ);
        size--;
        compNumberOfLeaves(nodes[0]);
    }
    (*root) = nodes[0];
}

static void updatedDistWPGMA(double ***matrix, int joinI, int joinJ, int size) {
    double **newMatrix = (double **) calloc((size_t) size - 1, sizeof(double *));
    for (int i = 0; i < size - 1; i++) {
        newMatrix[i] = (double *) calloc((size_t) size - 1, sizeof(double *));
    }
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - 1; j++) {
            if (i == j) {
                newMatrix[i][j] = 0.0;
            } else if (i < j) {
                if (i == 0) {
                    int indJ = j - (j <= joinI) + (j >= joinJ);
                    newMatrix[0][j] = ((*matrix)[joinI][indJ] + (*matrix)[joinJ][indJ]) / 2.0;
                } else {
                    int indI = i - (i <= joinI) + (i >= joinJ);
                    int indJ = j - (j <= joinI) + (j >= joinJ);
                    newMatrix[i][j] = (*matrix)[indI][indJ];
                }
            } else {
                newMatrix[i][j] = newMatrix[j][i];
            }
        }
    }
    for (int i = 0; i < size; i++) {
        free((*matrix)[i]);
    }
    free(*matrix);
    (*matrix) = newMatrix;
}

static void updatedDistUPGMA(double ***matrix, struct node **nodes, int joinI, int joinJ, int size) {
    double **newMatrix = (double **) calloc((size_t) size - 1, sizeof(double *));
    for (int i = 0; i < size - 1; i++) {
        newMatrix[i] = (double *) calloc((size_t) size - 1, sizeof(double *));
    }
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - 1; j++) {
            if (i == j) {
                newMatrix[i][j] = 0.0;
            } else if (i < j) {
                if (i == 0) {
                    int indJ = j - (j <= joinI) + (j >= joinJ);
                    int numberI = nodes[joinI]->numberOfLeaves;
                    int numberJ = nodes[joinJ]->numberOfLeaves;
                    newMatrix[0][j] = ((*matrix)[joinI][indJ] * numberI + (*matrix)[joinJ][indJ] * numberJ) / (numberI + numberJ);
                } else {
                    int intI = i - (i <= joinI) + (i >= joinJ);
                    int intJ = j - (j <= joinI) + (j >= joinJ);
                    newMatrix[i][j] = (*matrix)[intI][intJ];
                }
            } else {
                newMatrix[i][j] = newMatrix[j][i];
            }
        }
    }
    for (int i = 0; i < size; i++) {
        free((*matrix)[i]);
    }
    free(*matrix);
    (*matrix) = newMatrix;
}

static void calcDistUWPGMA(double dist, struct node *nodeI, struct node *nodeJ, double *distI, double*distJ) {
    double totalI = 0.0;
    double totalJ = 0.0;
    struct node *current = nodeI;
    while (current->firstChild != 0) {
        current = current->firstChild;
        totalI += current->distToParent;
    }
    current = nodeJ;
    while (current->firstChild != 0) {
        current = current->firstChild;
        totalJ += current->distToParent;
    }
    *distI = (dist / 2.0) - totalI;
    *distJ = (dist / 2.0) - totalJ;
}

static void findMin(double **matrix, int size, int *i, int *j) {
    if (size > 0) {
        double min = DBL_MAX;
        for (int row = 0; row < size; row++) {
            for (int col = 0; col < size; col++) {
                if (row != col && matrix[row][col] < min) {
                    min = matrix[row][col];
                    *i = row;
                    *j = col;
                }
            }
        }
    }
}

static void join(struct node ***nodes, int size, int i, int j, double distI, double distJ) {
    struct node **newNodes = (struct node **) calloc((size_t) size - 1, sizeof(struct node *));
    newNodes[0] = (struct node *) calloc(1, sizeof(struct node));
    newNodes[0]->name[0] = placeholderName;
    newNodes[0]->firstChild = (*nodes)[i];
    (*nodes)[i]->distToParent = distI;
    (*nodes)[i]->parent = newNodes[0];
    (*nodes)[i]->nextSibling = (*nodes)[j];
    (*nodes)[j]->distToParent = distJ;
    (*nodes)[j]->parent = newNodes[0];
    for (int node = 1; node < size - 1; node++) {
        int intI = node - (node <= i) + (node >= j);
        newNodes[node] = (*nodes)[intI];
    }
    free(*nodes);
    *nodes = newNodes;
}
