#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "node.h"
#include "tree.h"
#include "parse.h"

static int check(int index, int index2, int index3, int index4, char **names);
static void score(struct node *thisNode, int index, char **names);
void astralTag(struct node *root);
void astralRoot(struct node **root);
static void reRoot(struct node *thisNode, int *bestScore, struct node **bestRoot);
static void traverse(struct node *newTree, struct node *oldTree, struct node *avoid);
void mad(struct node **root);
static void nodeToNodeDist(struct node *current, double **dist, int index, int size, int *leaveIndex, int *leaveToId);
static void relDev(struct node *current, double **relDeviation, double **dist, int index, int *leaveIndices, int leaveIndex);
static void calcRho(double **dist, int from, int to, int indexJ, int size, double *rho, int *leafIndices);
static void calcScore(double **nodeDist, double **relDev, int from, int to, double disIJ, int size, double *rho, double *score, int *bestId, int *leafIndices, int indexNode);
static void reRoot_MAD(struct node **root, int bestId, double rho);

static void nodeToNodeDist(struct node *current, double **dist, int index, int size, int *leaveIndex, int *leaveToId) {
    if (current->numberOfChildren != 0) {
        struct node *child = current->firstChild;
        index++;
        for (int childNo = 0; childNo < current->numberOfChildren; childNo++) {
            for (int i = 0; i < index; i++) {
                for (int j = index; j < index + child->numberOfNodes; j++) {
                    dist[i][j - 1 - i] += child->distToParent;
                }
            }
            for (int i = index; i < index + child->numberOfNodes; i++) {
                for (int j = index + child->numberOfNodes; j < size; j++) {
                    dist[i][j - 1 - i] += child->distToParent;
                }
            }
            nodeToNodeDist(child, dist, index, size, leaveIndex, leaveToId);
            index += child->numberOfNodes;
            child = child->nextSibling;
        }
    } else {
        leaveToId[*leaveIndex] = index;
        *leaveIndex += 1;
    }
}

static void relDev(struct node *current, double **relDeviation, double **dist, int index, int *leafIndices, int leafIndex) {
    if (current->numberOfChildren != 0) {
        struct node *child = current->firstChild;
        int childIndex = index + 1;
        while (child != 0) {
            struct node *next = child->nextSibling;
            int nextIndex = leafIndex + child->numberOfLeaves;
            while (next != 0) {
                for (int i = leafIndex; i < leafIndex + child->numberOfLeaves; i++) {
                    for (int j = nextIndex; j < nextIndex + next->numberOfLeaves; j++) {
                        relDeviation[i][j - 1 - i] = fabs(((2.0 * dist[index][leafIndices[i] - 1 - index]) / dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]]) - 1.0);
//                        printf("relDev for:\t%d\t%d\t%d\t\t%.5f\t%.5f\t%.5f\n", index, leafIndices[i], leafIndices[j], dist[index][leafIndices[i] - 1 - index], dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]], relDeviation[i][j - 1 - i]);
//                        printf("%d\t%d\t%f\n", i, j, fabs(((2.0 * dist[index][leafIndices[i] - 1 - index]) / dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]]) - 1.0));
                    }
                }
                nextIndex += next->numberOfLeaves;
                next = next->nextSibling;
            }
            relDev(child, relDeviation, dist, childIndex, leafIndices, leafIndex);
            childIndex += child->numberOfNodes;
            leafIndex += child->numberOfLeaves;
            child = child->nextSibling;
            
        }
    }
}

static void calcRho(double **dist, int from, int to, int indexJ, int size, double *rho, int *leafIndices) {
    // from == i and is child node of j
    double divider = 0.0;
    *rho = 0.0;
    
    for (int i = from; i < to; i++) {
        for (int j = 0; j < from; j++) {
            divider += pow(dist[leafIndices[j]][leafIndices[i] - 1 - leafIndices[j]], -2.0);
        }
        for (int j = to; j < size; j++) {
            divider += pow(dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]], -2.0);
        }
    }
    divider *= 2.0 * dist[indexJ][leafIndices[from] - 1 - indexJ];
    divider = 1.0 / divider;
    
    for (int i = 0; i < from; i++) {
        for (int j = from; j < to; j++) {
            double thisDist = (j > from) ? dist[from][j - 1 - from] : 0.0;
            *rho += ((dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]] - 2.0 * thisDist) * pow(dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]], -2.0) * divider);
        }
    }
    for (int i = from; i < to; i++) {
        double thisDist = (i > from) ? dist[from][i - 1 - from] : 0.0;
        for (int j = to; j < size; j++) {
            *rho += ((dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]] - 2.0 * thisDist) * pow(dist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]], -2.0) * divider);
        }
    }
    
    if (*rho < 0.0) {
        *rho = 0.0;
    } else if (*rho > 1.0) {
        *rho = 1.0;
    }
}

static void calcScore(double **nodeDist, double **relDev, int from, int to, double disIJ, int size, double *rho, double *score, int *bestId, int *leafIndices, int indexNode) {
    double thisRho = 0.0;
    calcRho(nodeDist, from, to, 0, size, &thisRho, leafIndices);
//    printf("%d:%d:\trho: %f\t\t", from, to, thisRho);
    double thisScore = 0.0;
    
    for (int i = 0; i < from; i++) {
        for (int j = i + 1; j < from; j++) {
            thisScore += pow(relDev[i][j - 1 - i], 2.0);
        }
        for (int j = to; j < size; j++) {
            thisScore += pow(relDev[i][j - 1 - i], 2.0);
        }
    }
    for (int i = from; i < to; i++) {
        for (int j = i + 1; j < to; j++) {
            thisScore += pow(relDev[i][j - 1 - i], 2.0);
        }
    }
    for (int i = to; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            thisScore += pow(relDev[i][j - 1 - i], 2.0);
        }
    }
    
//    printf("score: ext: %f\t\t", thisScore);
    double disIRho = thisRho * disIJ;
    for (int i = from; i < to; i++) {
        double thisDis = (leafIndices[i] > indexNode) ? nodeDist[indexNode][leafIndices[i] - 1 - indexNode] : 0.0;
        thisDis += disIRho;
        thisDis *= 2;
        for (int j = 0; j < from; j++) {
            thisScore += pow((thisDis / nodeDist[leafIndices[j]][leafIndices[i] - 1 - leafIndices[j]]) - 1.0, 2.0);
        }
        for (int j = to; j < size; j++) {
            thisScore += pow((thisDis / nodeDist[leafIndices[i]][leafIndices[j] - 1 - leafIndices[i]]) - 1.0, 2.0);
        }
    }
    if (thisScore < *score) {
        *score = thisScore;
        *bestId = indexNode;
        *rho = thisRho;
    }
//    printf("%f\n", thisScore);
    
}

static void reRoot_MAD(struct node **root, int bestId, double rho) {
    struct node *current = (*root);
    int id = 0;
    while (1 == 1) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            id += 1;
            if (id == bestId) {
                break;
            }
        }
        if (id == bestId) {
            break;
        }
        while (current->nextSibling == 0) {
            current = current->parent;
        }
        current = current->nextSibling;
        id += 1;
        if (id == bestId) {
            break;
        }
    }
    
    struct node *newRoot = (struct node *) calloc(1, sizeof(struct node));
    newRoot->name[0] = placeholderName;
    
    newRoot->firstChild = (struct node *) calloc(1, sizeof(struct node));
    strcpy(newRoot->firstChild->name, current->name);
    newRoot->firstChild->parent = newRoot;
    newRoot->firstChild->distToParent = current->distToParent * rho;
    traverse(newRoot->firstChild, current, current->parent);
    
    newRoot->firstChild->nextSibling = (struct node *) calloc(1, sizeof(struct node));
    strcpy(newRoot->firstChild->nextSibling->name, current->parent->name);
    newRoot->firstChild->nextSibling->parent = newRoot;
    newRoot->firstChild->nextSibling->distToParent = current->distToParent * (1 - rho);
    traverse(newRoot->firstChild->nextSibling, current->parent, current);
    
    compNumberOfLeaves(newRoot);
    astralTag(newRoot);
    
    freeTree(*root);
    (*root) = newRoot;
    
    
}

void mad(struct node **root) {
    int size = (*root)->numberOfNodes;
    int *leafIndices = (int *) calloc((size_t) (*root)->numberOfLeaves, sizeof(int));
    int index = 0;
    double **nodeToNodeDistances = (double **) calloc((size_t) (size - 1), sizeof(double *));
    for (int i = 0; i < size - 1; i++) {
        nodeToNodeDistances[i] = (double *) calloc((size_t) (size - 1 - i), sizeof(double));
    }
    
    nodeToNodeDist((*root), nodeToNodeDistances, 0, size, &index, leafIndices);
    
    size = (*root)->numberOfLeaves;
    
    double **relDeviation = (double **) calloc((size_t) (size - 1), sizeof(double *));
    for (int i = 0; i < size - 1; i++) {
        relDeviation[i] = (double *) calloc((size_t) (size - 1 - i), sizeof(double));
    }
    
    relDev((*root), relDeviation, nodeToNodeDistances, 0, leafIndices, 0);
    
//    for (int i = 0; i < size - 1; i++) {
//        for (int j = 0; j < size - 1 - i; j++) {
//            printf("%.3f  ", relDeviation[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
    
    double rho = 0.0;
    double score = DBL_MAX;
    int bestId = 0;

    struct node *current = (*root);
    int indexNode = 0;
    int indexLeaf = 0;
    while (1 == 1) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            indexNode += 1;
            calcScore(nodeToNodeDistances, relDeviation, indexLeaf, indexLeaf + current->numberOfLeaves, current->distToParent, size, &rho, &score, &bestId, leafIndices, indexNode);
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current == (*root)) {
                break;
            }
        }
        if (current == (*root)) {
            break;
        }
        indexLeaf += 1;
        indexNode += 1;
        current = current->nextSibling;
        calcScore(nodeToNodeDistances, relDeviation, indexLeaf, indexLeaf + current->numberOfLeaves, current->distToParent, size, &rho, &score, &bestId, leafIndices, indexNode);
    }
    
//    printf("WINNER:\n%d\n%f\n", bestId, score);
    
    reRoot_MAD(root, bestId, rho);
    
//    printTree((*root), 5);
    
    for (int i = 0; i < size - 1; i++) {
        free(nodeToNodeDistances[i]);
        free(relDeviation[i]);
    }
    if (relDeviation != 0) {
        free(relDeviation);
    }
    if (nodeToNodeDistances != 0) {
        free(nodeToNodeDistances);
    }
}

void astralRoot(struct node **root) {
    //go to every node
    struct node *current = (*root);
    struct node *bestRoot = 0;
    int bestScore = INT_MAX;
    do {
        while (current->firstChild != 0) {
            current = current->firstChild;
            reRoot(current, &bestScore, &bestRoot);
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current == (*root)) {
                break;
            }
        }
        if (current == (*root)) {
            break;
        }
        current = current->nextSibling;
        reRoot(current, &bestScore, &bestRoot);
    } while (current != (*root));
    
    freeTree(*root);
    (*root) = bestRoot;
}

static void reRoot(struct node *thisNode, int *bestScore, struct node **bestRoot) {
    struct node *newRoot = (struct node *) calloc(1, sizeof(struct node));
    newRoot->name[0] = placeholderName;
    
    newRoot->firstChild = (struct node *) calloc(1, sizeof(struct node));
    strcpy(newRoot->firstChild->name, thisNode->name);
    newRoot->firstChild->parent = newRoot;
    newRoot->firstChild->distToParent = thisNode->distToParent / 2.0;
    traverse(newRoot->firstChild, thisNode, thisNode->parent);
    
    newRoot->firstChild->nextSibling = (struct node *) calloc(1, sizeof(struct node));
    strcpy(newRoot->firstChild->nextSibling->name, thisNode->parent->name);
    newRoot->firstChild->nextSibling->parent = newRoot;
    newRoot->firstChild->nextSibling->distToParent = thisNode->distToParent / 2.0;
    traverse(newRoot->firstChild->nextSibling, thisNode->parent, thisNode);
    
    compNumberOfLeaves(newRoot);
    astralTag(newRoot);
    if (newRoot->score < *bestScore) {
        *bestScore = newRoot->score;
        if (*bestRoot != 0) {
            freeTree(*bestRoot);
        }
        (*bestRoot) = newRoot;
//        printTree(newRoot, 6);
    } else {
        freeTree(newRoot);
    }
}

static void traverse(struct node *newTree, struct node *oldTree, struct node *avoid) {
    if (oldTree->parent != 0 && oldTree->parent != avoid) {
        newTree->firstChild = (struct node *) calloc(1, sizeof(struct node));
        newTree->firstChild->parent = newTree;
        strcpy(newTree->firstChild->name, oldTree->parent->name);
        newTree->firstChild->distToParent = oldTree->distToParent;
        traverse(newTree->firstChild, oldTree->parent, oldTree);
    }

    struct node *child = oldTree->firstChild;
    struct node *before = newTree->firstChild;
    while (child != 0) {
        if (child != avoid) {
            struct node *new = (struct node *) calloc(1, sizeof(struct node));
            new->parent = newTree;
            strcpy(new->name, child->name);
            new->distToParent = child->distToParent;
            if (before == 0) {
                newTree->firstChild = new;
            } else {
                before->nextSibling = new;
            }
            traverse(new, child, oldTree);
            before = new;
        }
        child = child->nextSibling;
    }
}

void astralTag(struct node *root) {
    int size = root->numberOfLeaves;
    char **names = (char **) calloc(size, sizeof(char *));
    for (int i = 0; i < size; i++) {
        names[i] = (char *) calloc(maxNameLength + 1, sizeof(char));
    }
    score(root, 0, names);
    for (int i = 0; i < size; i++) {
        free(names[i]);
    }
    if (names != 0) {
        free(names);
    }
}

static void score(struct node *thisNode, int index, char **names) {
    thisNode->score = 0;
    thisNode->tag = 0;
    if (thisNode->firstChild != 0) {
        struct node *next = thisNode->firstChild;
        int nextIndex = index;
        while (next != 0) {
            score(next, nextIndex, names);
            nextIndex += next->numberOfLeaves;
            next = next->nextSibling;
        }
    } else {
        strcpy(names[index], thisNode->name);
    }
    
    struct node *child = thisNode->firstChild;
    int index1 = index;
    
    while (child != 0) {
        int index2 = index1 + child->numberOfLeaves;
        struct node *next = child->nextSibling;
        int index3 = index2;
        while (next != 0) {
            int index4 = index3 + next->numberOfLeaves;
            int score = check(index1, index2, index3, index4, names);
            if (score > 0) {
                child->parent->tag = 1;
                child->parent->score += score;
            }
            next = next->nextSibling;
            index3 = index4;
        }
        child->parent->score += child->score;
        child = child->nextSibling;
        index1 = index2;
    }
}

static int check(int index, int index2, int index3, int index4, char **names) {
    int score = 0;
    int overlap = 0;
    int aNotInB = 0;
    int bNotInA = 0;
    
    for (int i = index; i < index2; i++) {
        int hit = 0;
        for (int j = index3; j < index4; j++) {
            if (strcmp(names[i], names[j]) == 0) {
                overlap++;
                hit = 1;
                break;
            }
        }
        if (hit == 0) {
            aNotInB = 1;
            if (overlap > 0) {
                break;
            }
        }
    }
    for (int i = index3; i < index4; i++) {
        int hit = 0;
        for (int j = index; j < index2; j++) {
            if (strcmp(names[i], names[j]) == 0) {
                overlap++;
                hit = 1;
                break;
            }
        }
        if (hit == 0) {
            bNotInA = 1;
            if (overlap > 0) {
                break;
            }
        }
    }
    if (aNotInB == 0 && bNotInA == 0) {
        score = 1;
    } else if (aNotInB == 0 || bNotInA == 0) {
        score = 2;
    } else if (overlap > 0) {
        score = 3;
    }
    return score;
}
