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
//static void copySubtree(struct node *from, struct node *at);
static int check(char **names, int index1, int index2, int index3, int *aInB, int *bInA);
static int subScoreAndTag(struct node *current, int *index, char **names, int *aInB, int *bInA);
int scoreAndTag(struct node *tree, char **names);
//static void reRoot(struct node *atNode, struct node **best, int *bestScore);
void tagAndRoot(struct node *tree);
static void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size);
static void goingDown(int index, int numberOfLeaves, double dist, double **allDist);
void leafToLeafDistance(struct node *root, double **dist, char **name, int branchLength);
void deletedTaggedDistance(struct node *current, double **dist, int *index);
int compNumberOfNodes(struct node *tree);
static void setScoreTagLeavesZero(struct node *tree);
int scoreAndTag2(struct node *tree);
static int subScoreAndTag2(struct node *current, char **names, int *index);
static int check2(char **names, int index0, int index1, int index2, int index3);
static void goingUp2(double **dist, double distance, int index, int size, int maxIndex);
static void goingDown2(double **dist, double distance, int index, int size);
static void nodeToNodeDistance(struct node *root, double **dist);
static void leafToLeafRelativeDeviation(struct node *root, double **dist, double **relDev, int *indexOfId);
static void calcRelDiv(struct node *current, double **dist, double **relDev, int *indexOfId);
static double rms(struct node *root, struct node *node, double **dist, double **relDev, int *indexOfId, double *rho);
static void findNodeByID(struct node *root, int id, struct node **nodeWithId);
static void removeSibling(struct node *sibling, struct node *parent);
static void madRoot(struct node **root, int topId, double rho);
void mad(struct node **root);

int compNumberOfNodes(struct node *tree) {
    int number = 1;
    if (tree->firstChild != 0) {
        number += compNumberOfNodes(tree->firstChild);
    }
    if (tree->nextSibling != 0) {
        number += compNumberOfNodes(tree->nextSibling);
    }
    return number;
}

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

static int check(char **names, int index1, int index2, int index3, int *aInB, int *bInA) {
    
    // Score according to Astral-Pro tag procedure
    int score = 0;
    int aInBCount = 0;
    int bInACount = 0;
    
    for (int i = index1; i < index2; i++) {
        for (int j = index2; j <= index3; j++) {
            if (aInB[i - index1] == 0 || bInA[j - index2] == 0) {
                if (strcmp(names[i], names[j]) == 0) {
                    aInB[i - index1] = 1;
                    bInA[j - index2] = 1;
                }
            }
        }
    }
    for (int i = 0; i < (index2 - index1); i++) {
        aInBCount += aInB[i];
        aInB[i] = 0;
    }
    for (int i = 0; i < (index3 - index2 + 1); i++) {
        bInACount += bInA[i];
        bInA[i] = 0;
    }
    if (aInBCount == (index2 - index1)) {
        if (bInACount == (index3 - index2 + 1)) {
            score = 1;
        } else {
            score = 2;
        }
    } else if (bInACount == (index3 - index2 + 1)) {
        score = 2;
    } else if (aInBCount == 0 && bInACount == 0){
        score = 0;
    } else {
        score = 3;
    }
    return score;
}


static int subScoreAndTag(struct node *current, int *index, char **names, int *aInB, int *bInA) {
    int score = 0;
    if (current->firstChild != 0) {
        struct node *this = current->firstChild;
        int index1 = *index;
        score += subScoreAndTag(this, index, names, aInB, bInA);
        struct node *next = this->nextSibling;
        int a = 0;
        while (1) {
            if (next == 0) {
                break;
            }
            *index += 1;
            int index2 = *index;
            score += subScoreAndTag(next, index, names, aInB, bInA);
            int index3 = *index;
            if (current->score == 0) {
                int checker = check(names, index1, index2, index3, aInB, bInA);
                if (checker != 0) {
                    score += checker;
                    current->tag = 1;
                }
            } else {
                score = current->score;
                if (score == -1) {
                    score = 0;
                }
            }
            this = next;
            next = next->nextSibling;
            index1 = index2;
            a++;
        }
        current->score = score;
        if (current->score == 0) {
            current->score = -1;
        }
    } else {
        strcpy(names[*index], current->name);
    }
    return score;
}

int scoreAndTag(struct node *tree, char **names) {
    // Tag procedre from Astral-Pro
    int index = 0;
    int score = 0;
    int *aInB = (int *) calloc(sizeof(int), tree->numberOfLeaves);
    int *bInA = (int *) calloc(sizeof(int), tree->numberOfLeaves);
    score = subScoreAndTag(tree, &index, names, aInB, bInA);
    free(aInB);
    free(bInA);
    return score;
}

int check2(char **names, int index0, int index1, int index2, int index3) {
    int score = 0;
    int aInB = 0;
    int bInA = 0;
    for (int i = index0; i < index1; i++) {
        for (int j = index2; j < index3; j++) {
            if (strcmp(names[i], names[j]) == 0) {
                aInB += 1;
                break;
            }
        }
    }
    
    for (int i = index2; i < index3; i++) {
        for (int j = index0; j < index1; j++) {
            if (strcmp(names[i], names[j]) == 0) {
                bInA += 1;
                break;
            }
        }
    }
    
    if (aInB == (index1 - index0) || bInA == (index3 - index2)) {
        if (aInB == (index1 - index0) && bInA == (index3 - index2)) {
            score = 1;
        } else {
            score = 2;
        }
    } else if (aInB != 0) {
        score = 3;
    }
    return score;
}

int subScoreAndTag2(struct node *current, char **names, int *index) {
    int score = 0;
    if (current->firstChild != 0) {
        int children = current->numberOfChildren;
        int *indices = calloc(sizeof(int), children + 1);
        int i = 0;
        indices[0] = *index;
        struct node *workOn = current->firstChild;
        while (workOn != 0) {
            i += 1;
            subScoreAndTag2(workOn, names, index);
            score += workOn->score;
            workOn = workOn->nextSibling;
            indices[i] = *index + 1;
            if (workOn != 0) {
                *index += 1;
            }
        }
        for (i = 0; i < children - 1; i++) {
            for (int j = i + 1; j < children; j++) {
                int checker = check2(names, indices[i], indices[i + 1], indices[j], indices[j + 1]);
                if (checker != 0) {
                    (current->tag2)[i * children + j - (((i + 1) * (i + 2)) / 2)] = 1;
                } else {
                    (current->tag2)[i * children + j - (((i + 1) * (i + 2)) / 2)] = 0;
                }
                score += checker;
            }
        }
        
    } else {
        strcpy(names[*index], current->name);
    }
    if (score != 0) {
        current->score = score;
    }
    return score;
}

int scoreAndTag2(struct node *tree) {
    int size = tree->numberOfLeaves;
    char **names = calloc(sizeof(char *), size);
    for (int i = 0; i < size; i++) {
        names[i] = calloc(sizeof(char), maxNameLength);
    }
    int *index = calloc(sizeof(int), 1);
    int score = subScoreAndTag2(tree, names, index);
    for (int i = 0; i < size; i++) {
        free(names[i]);
    }
    free(names);
    return score;
}




void setScoreTagLeavesZero(struct node *tree) {
    struct node *current = tree;
    current->tag = 0;
    current->score = 0;
    current->numberOfLeaves = 0;
    while (current->firstChild != 0) {
        current = current->firstChild;
        current->tag = 0;
        current->score = 0;
        current->numberOfLeaves = 0;
    }
    while (current->nextSibling == 0) {
        current = current->parent;
        if (current == tree) {
            return;
        }
    }
    current = current->nextSibling;
    current->tag = 0;
    current->score = 0;
    current->numberOfLeaves = 0;
}

static void reRoot2(struct node *tree, struct node *root) {
    if (tree->parent != root) {
        struct node *current = tree;
        struct node *newRoot = (struct node *) calloc(sizeof(struct node), 1);
        struct node *toParent = newRoot;
        toParent->firstChild = current;
        
        struct node *currentSibling = current->parent->firstChild;
        if (currentSibling == current) {
            if (currentSibling->nextSibling != 0) {
                current->parent->firstChild = currentSibling->nextSibling;
            }
        } else {
            while (currentSibling->nextSibling != 0) {
                if (currentSibling->nextSibling == current) {
                    if (currentSibling->nextSibling->nextSibling != 0) {
                        currentSibling->nextSibling = currentSibling->nextSibling->nextSibling;
                        break;
                    } else {
                        currentSibling->nextSibling = 0;
                        break;
                    }
                }
                currentSibling = currentSibling->nextSibling;
            }
        }
        
        struct node *tmpCurrent = currentSibling->parent->parent;
        struct node *tmpParent = currentSibling->parent;
        current->nextSibling = currentSibling->parent;
//        current->nextSibling->score = 0;
//        current->nextSibling->tag = 0;
//        current->nextSibling->numberOfLeaves = 0;
        current->nextSibling->parent = toParent;
        while (1) {
            current = tmpCurrent;
            toParent = tmpParent;
            if (current == root) {
                current = root->firstChild;
                currentSibling = toParent->firstChild;
                while (currentSibling->nextSibling != 0) {
                    currentSibling = currentSibling->nextSibling;
                }
                while (current != 0) {
                    if (current != toParent) {
                        current->parent = toParent;
                        currentSibling->nextSibling = current;
                        currentSibling = currentSibling->nextSibling;
                    }
                    if (current->nextSibling == toParent) {
                        current->nextSibling = current->nextSibling->nextSibling;
                    }
                    current = current->nextSibling;
                    
                }
                toParent->nextSibling = 0;
                break;
            }
            currentSibling = current->firstChild;
            if (currentSibling == toParent) {
                if (currentSibling->nextSibling != 0) {
                    current->firstChild = currentSibling->nextSibling;
                    currentSibling->nextSibling = 0;
                }
            } else {
                while (currentSibling->nextSibling != 0) {
                    if (currentSibling->nextSibling == toParent) {
                        struct node *tmp = currentSibling->nextSibling->nextSibling;
                        if (tmp != 0) {
                            currentSibling->nextSibling->nextSibling = 0;
                            currentSibling->nextSibling = tmp;
                            break;
                        } else {
                            currentSibling->nextSibling = 0;
                            break;
                        }
                    }
                    currentSibling = currentSibling->nextSibling;
                }
            }
            tmpParent = current;
            tmpCurrent = current->parent;
            current->parent = toParent;
//            current->score = 0;
//            current->tag = 0;
//            current->numberOfLeaves = 0;
            currentSibling = toParent->firstChild;
            while (currentSibling->nextSibling != 0) {
                currentSibling = currentSibling->nextSibling;
            }
            currentSibling->nextSibling = current;
            toParent = current;
            current = tmpCurrent;
        }
        root->firstChild = tree;
        root->firstChild->parent = root;
        root->firstChild->nextSibling->parent = root;
//        root->score = 0;
//        root->tag = 0;
//        root->numberOfLeaves = 0;
        free(newRoot);
        setScoreTagLeavesZero(root);
        compNumberOfLeaves(root);
    }
//    printTree(root);
}

static void addRoot(struct node *tree) {
    struct node *newNode = calloc(sizeof(struct node), 1);
    newNode->parent = tree;
    newNode->numberOfChildren = 2;
    newNode->numberOfLeaves = tree->numberOfLeaves - tree->firstChild->numberOfLeaves;
    newNode->firstChild = tree->firstChild->nextSibling;
    newNode->name[0] = placeholderName;
    newNode->tag2 = calloc(sizeof(int), 1);
    tree->firstChild->nextSibling->parent = newNode;
    tree->firstChild->nextSibling->nextSibling->parent = newNode;
    tree->firstChild->nextSibling = newNode;
    tree->numberOfChildren = 2;
}

void tagAndRoot(struct node *tree) {
    struct node *current = tree;
    int bestScore = INT_MAX;;
    int score = INT_MAX;;
    int bestId = 0;
    int id = 1;
    int maxId = compNumberOfNodes(tree);
    int size = tree->numberOfLeaves;
    char **names = (char **) calloc(sizeof(char *), size);
    for (int i = 0; i < size; i++) {
        names[i] = (char *) calloc(sizeof(char), maxNameLength);
    }
    
    if (tree->numberOfChildren == 3) {
        addRoot(tree);
    }
    
    while (id < maxId) {
//        printf("id: %d\n", id);
        while (current->firstChild != 0) {
            current = current->firstChild;
            if (current->idNo == id) {
                reRoot2(current, tree);
                score = scoreAndTag2(tree);
                if (score < bestScore) {
                    bestScore = score;
                    bestId = id;
                }
                id++;
                if (score == 0) {
                    break;
                }
            }
        }
        if (score == 0) {
            break;
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current->parent == 0) {
                break;
            }
        }
        if (current->parent != 0) {
            current = current->nextSibling;
            if (current->idNo == id) {
                reRoot2(current, tree);
                score = scoreAndTag2(tree);
                if (score < bestScore) {
                    bestScore = score;
                    bestId = id;
                }
                id++;
            }
        }
    }
    current = tree;
    while (score != 0) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            if (current->idNo == bestId) {
                reRoot2(current, tree);
                scoreAndTag2(tree);
                break;
            }
        }
        if (tree->firstChild->idNo == bestId) {
            break;
        }
        while (current->nextSibling == 0) {
            current = current->parent;
        }
        current = current->nextSibling;
        if (current->idNo == bestId) {
            reRoot2(current, tree);
            scoreAndTag2(tree);
            break;
        }
    }
    for (int i = 0; i < size; i++) {
        free(names[i]);
    }
    free(names);
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

void deletedTaggedDistance(struct node *current, double **dist, int *index) {
    int children = current->numberOfChildren;
    if (children != 0) {
        struct node *workOn = current->firstChild;
        for (int i = 0; i < children; i++) {
            int index0 = *index;
            int index1 = index0 + workOn->numberOfLeaves;
            int index2 = index1;
            struct node *workOn2 = workOn->nextSibling;
            for (int j = i + 1; j < children; j++) {
                int index3 = index2 + workOn2->numberOfLeaves;
                int index = (i * children + j) - (((i + 1) * (i + 2)) / 2);
                if (current->tag2[index] != 0) {
                    for (int k = index0; k < index1; k++) {
                        for (int l = index2; l < index3; l++) {
                            dist[k][l - k] = 0;
                        }
                    }
                }
                index2 += workOn2->numberOfLeaves;
                workOn2 = workOn2->nextSibling;
            }
            deletedTaggedDistance(workOn, dist, index);
            workOn = workOn->nextSibling;
        }
    } else {
        *index += 1;
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
//            printf("go down\n");
            current = current->firstChild;
            if (!branchLength) {
                goingDown(index, current->numberOfLeaves, 1.0, dist);
            } else {
                goingDown(index, current->numberOfLeaves, current->distToParent, dist);
            }
        }
        
        strcpy(name[index], current->name);
        
        while (current->nextSibling == 0) {
//            printf("go up\n");
            if (index < size - 1) {
                if (!branchLength) {
                    goingUp(index, current->numberOfLeaves, 1.0, dist, size);
                } else {
                    goingUp(index, current->numberOfLeaves, current->distToParent, dist, size);
                }
            }
            
//            struct node *toFree = current;
            current = current->parent;
//            free(toFree);
            if (current == root) {
                break;
            }
        }
        if (current == root) {
//            free(current);
            break;
        }
        if (!branchLength) {
            goingUp(index, current->numberOfLeaves, 1.0, dist, size);
        } else {
            goingUp(index, current->numberOfLeaves, current->distToParent, dist, size);
        }
        
        index++;
//        struct node *toFree = current;
//        printf("go next\n");
        current = current->nextSibling;
//        free(toFree);
        if (!branchLength) {
            goingDown(index, current->numberOfLeaves, 1.0, dist);
        } else {
            goingDown(index, current->numberOfLeaves, current->distToParent, dist);
        }
    }
}

static void goingDown2(double **dist, double distance, int index, int size) {
    for (int i = 0; i < index; i++) {
        for (int j = index; j < index + size; j++) {
            dist[i][j] += distance;
        }
    }
}

static void goingUp2(double **dist, double distance, int index, int size, int maxIndex) {
    for (int i = index; i < index + size; i++) {
        for (int j = index + size; j < maxIndex; j++) {
            dist[i][j] += distance;
        }
    }
}

static void nodeToNodeDistance(struct node *root, double **dist) {
    
    int maxIndex = root->numberOfNodes;
    
    struct node *current = root;
    while (1) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            goingDown2(dist, current->distToParent, current->idNo, current->numberOfNodes);
        }
        while (current->nextSibling == 0) {
            goingUp2(dist, current->distToParent, current->idNo, current->numberOfNodes, maxIndex);
            current = current->parent;
            if (current == root) {
                return;
            }
        }
        goingUp2(dist, current->distToParent, current->idNo, current->numberOfNodes, maxIndex);
        current = current->nextSibling;
        goingDown2(dist, current->distToParent, current->idNo, current->numberOfNodes);
    }
}


static void calcRelDiv(struct node *current, double **dist, double **relDev, int *indexOfId) {
    // check if internal node
    if (current->numberOfChildren == 0) {
        return;
    }
    struct node *first = current->firstChild;
    struct node *second;
    
    //go through all children
    for (int i = 0; i < current->numberOfChildren - 1; i++) {
        second = first->nextSibling;
        for (int j = i + 1; j < current->numberOfChildren; j++) {
            
            for (int k = first->idNo; k < first->idNo + first->numberOfNodes; k++) {
                if (indexOfId[k] == -1) {
                    continue;
                }
                for (int l = second->idNo; l < second->idNo + second->numberOfNodes; l++) {
                    if (indexOfId[l] == -1) {
                        continue;
                    }
                    relDev[indexOfId[k]][indexOfId[l]] = fabs(((2 * dist[current->idNo][k]) / dist[k][l]) - 1);
                }
            }
            second = second->nextSibling;
        }
        first = first->nextSibling;
    }
}

static void leafToLeafRelativeDeviation(struct node *root, double **dist, double **relDev, int *indexOfId) {
    // Traverse Tree and and calc relative deviation for internal nodes
    struct node *current = root;
    while (1) {
        calcRelDiv(current, dist, relDev, indexOfId);
        while (current->firstChild != 0) {
            current = current->firstChild;
            calcRelDiv(current, dist, relDev, indexOfId);
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current == root) {
                return;
            }
        }
        current = current->nextSibling;
    }
}

static double rms(struct node *root, struct node *node, double **dist, double **relDev, int *indexOfId, double *rho) {
    
    double rms = 0.0;
    double distIJ = node->distToParent;
    *rho = 0.0;
    
    int indexI_1 = node->idNo;
    int indexI_2 = node->idNo + node->numberOfNodes;
    int i = node->idNo;
    
    
    
    // Calculate rho
    
    for (int b = indexI_1; b < indexI_2; b++) {
        for (int c = 0; c < root->numberOfNodes; c++) {
            double tmp = (dist[b][c] - 2 * dist[b][i]) / pow(dist[b][c], 2.0);
            if ((c < indexI_1 || c >= indexI_2) && indexOfId[b] != -1 && indexOfId[c] != -1) {
                double tmp2 = 0.0;
                for (int bb = indexI_1; bb < indexI_2; bb++) {
                    for (int cc = 0; cc < root->numberOfNodes; cc++) {
                        if ((cc < indexI_1 || cc >= indexI_2) && indexOfId[bb] != -1 && indexOfId[cc] != -1) {
                            tmp2 += 1 / pow(dist[bb][cc], 2.0);
                        }
                    }
                }
                tmp = tmp / (2 * distIJ * tmp2);
                *rho += tmp;
            }
        }
    }
    if (*rho < 0.0) {
        *rho = 0.0;
    } else if (*rho > 1.0) {
        *rho = 1.0;
    }
    
    for (int b = 0; b < root->numberOfNodes; b++) {
        if (indexOfId[b] == -1) {
            continue;
        }
        for (int c = b + 1; c < root->numberOfNodes; c++) {
            if (indexOfId[c] == -1) {
                continue;
            }
            if ((b >= indexI_1 && b < indexI_2 && c >= indexI_1 && c < indexI_2) || (!(b >= indexI_1 && b < indexI_2) && !(c >= indexI_1 && c < indexI_2))) {
                rms += pow(relDev[indexOfId[b]][indexOfId[c]], 2.0);
//                printf("same side\t%d\t%d\t\t%f\n", b, c, relDev[indexOfId[b]][indexOfId[c]] * relDev[indexOfId[b]][indexOfId[c]]);
            } else {
                if (b < indexI_1 || b >= indexI_2) {
                    rms += pow((((2 * (dist[c][i] + *rho * distIJ)) / dist[b][c]) - 1), 2.0);
//                    printf("c closer\t%d\t%d\t\t%f\n", b, c, pow((((2 * (dist[c][i] + *rho * distIJ)) / dist[b][c]) - 1), 2.0));
                } else {
                    rms += pow((((2 * (dist[b][i] + *rho * distIJ)) / dist[b][c]) - 1), 2.0);
//                    printf("b closer\t%d\t%d\t\t%f\n", b, c, pow((((2 * (dist[b][i] + *rho * distIJ)) / dist[b][c]) - 1), 2.0));
                }
            }
        }
    }
//    printf("%f\t%f\t%d\n", rms, *rho, i);
    return rms;
}

static void findNodeByID(struct node *root, int id, struct node **nodeWithId) {
    *nodeWithId = root;
    while (1) {
        while ((*nodeWithId)->firstChild != 0) {
            (*nodeWithId) = (*nodeWithId)->firstChild;
            if ((*nodeWithId)->idNo == id) {
                return;
            }
        }
        while ((*nodeWithId)->nextSibling == 0) {
            (*nodeWithId) = (*nodeWithId)->parent;
        }
        (*nodeWithId) = (*nodeWithId)->nextSibling;
        if ((*nodeWithId)->idNo == id) {
            return;
        }
    }
}

static void removeSibling(struct node *sibling, struct node *parent) {
    if (parent->firstChild == sibling) {
        parent->firstChild = sibling->nextSibling;
    } else {
        struct node *child = parent->firstChild;
        while (child->nextSibling != sibling) {
            child = child->nextSibling;
        }
        child->nextSibling = sibling->nextSibling;
    }
    sibling->nextSibling = 0;
}

static void madRoot(struct node **root, int topId, double rho) {
    
    // Find node where new root should be inserted
    struct node *current;
    findNodeByID((*root), topId, &current);
    struct node *myParent = current->parent;
    // Create new root node
    struct node *newRoot = (struct node *) calloc(sizeof(struct node), 1);
    struct node *newParent = newRoot;
    
    // Set new parent
    current->parent = newParent;
    newParent->firstChild = current;
    
    // Remove as child
    removeSibling(current, myParent);
    
    // Add as sibling
    current->nextSibling = myParent;
    
    current = myParent->parent;
    myParent->parent = newParent;
    
    newParent = myParent;
    myParent = current->parent;
    
    
    while (1) {
        removeSibling(newParent, current);
        current->parent = newParent;
        newParent = newParent->firstChild;
        while (newParent->nextSibling != 0) {
            newParent = newParent->nextSibling;
        }
        newParent->nextSibling = current;
        if (current == (*root)) {
            break;
        }
        newParent = current;
        current = myParent;
        myParent = myParent->parent;
    }
    (*root) = newRoot;
    compNumberOfLeaves(*root);
}

void mad(struct node **root) {
    
    // Allocate space for index, distance, relative deviation, rho
    int *indexOfId = (int *) calloc(sizeof(int), (*root)->numberOfNodes);
    double **nodeToNodeDist = (double **) calloc(sizeof(double *), (*root)->numberOfNodes);
    double **relDev = (double **) calloc(sizeof(double *), (*root)->numberOfLeaves);
    double *rho = (double *) calloc(sizeof(double), 1);
    
    for (int j = 0; j < (*root)->numberOfNodes; j++) {
        nodeToNodeDist[j] = (double *) calloc(sizeof(double), (*root)->numberOfNodes);
    }
    for (int j = 0; j < (*root)->numberOfLeaves; j++) {
        relDev[j] = (double *) calloc(sizeof(double), (*root)->numberOfLeaves);
    }
    
    // Calc all node to node distance
    nodeToNodeDistance(*root, nodeToNodeDist);
    
    // Fill array
    for (int i = 0; i < (*root)->numberOfNodes; i++) {
        for (int j = 0; j < (*root)->numberOfNodes; j++) {
            if (j < i) {
                nodeToNodeDist[i][j] = nodeToNodeDist[j][i];
            }
        }
    }
    
//    for (int i = 0; i < (*root)->numberOfNodes; i++) {
//        for (int j = 0; j < (*root)->numberOfNodes; j++) {
//            printf("%f\t", nodeToNodeDist[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
    
    // Traverse Tree and save IDs of leaves
    struct node *current = *root;
    int index = 0;
    while (1) {
        if (current->firstChild == 0) {
            indexOfId[current->idNo] = index;
        } else {
            indexOfId[current->idNo] = -1;
        }
        while (current->firstChild != 0) {
            current = current->firstChild;
            if (current->firstChild == 0) {
                indexOfId[current->idNo] = index;
            } else {
                indexOfId[current->idNo] = -1;
            }
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current == *root) {
                break;
            }
        }
        if (current == *root) {
            break;
        }
        current = current->nextSibling;
        index += 1;
    }
    
    // Calc relative Deviation
    leafToLeafRelativeDeviation(*root, nodeToNodeDist, relDev, indexOfId);
    
    for (int i = 0; i < (*root)->numberOfLeaves; i++) {
        for (int j = 0; j < (*root)->numberOfLeaves; j++) {
            if (j < i) {
                relDev[i][j] = relDev[j][i];
            }
        }
    }
    
//    for (int i = 0; i < (*root)->numberOfLeaves; i++) {
//        for (int j = 0; j < (*root)->numberOfLeaves; j++) {
//            printf("%f\t", relDev[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
    
    
    
    // get scores for all possible roots
    double topScore = DBL_MAX;
    double topRho = 2;
    double score = DBL_MAX;
    int topId = 0;
    
    while (1) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            score = rms(*root, current, nodeToNodeDist, relDev, indexOfId, rho);
            if (score < topScore) {
                topScore = score;
                topRho = *rho;
                topId = current->idNo;
            }
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current == *root) {
                break;
            }
        }
        if (current == *root) {
            break;
        }
        current = current->nextSibling;
        score = rms(*root, current, nodeToNodeDist, relDev, indexOfId, rho);
        if (score < topScore) {
            topScore = score;
            topRho = *rho;
            topId = current->idNo;
        }
    }
   
    free(indexOfId);
    free(rho);
    
    for (int j = 0; j < (*root)->numberOfNodes; j++) {
        free(nodeToNodeDist[j]);
    }
    for (int j = 0; j < (*root)->numberOfLeaves; j++) {
        free(relDev[j]);
    }
    free(nodeToNodeDist);
    free(relDev);
    madRoot(root, topId, topRho);
}
