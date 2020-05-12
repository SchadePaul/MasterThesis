#include "node.h"
#include "parse.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

void freeTree(struct node *tree);
//static void copySubtree(struct node *from, struct node *at);
static int check(char **names, int index1, int index2, int index3);
static int subScoreAndTag(struct node *current, int *index, char **names);
int scoreAndTag(struct node *tree, char **names);
//static void reRoot(struct node *atNode, struct node **best, int *bestScore);
void tagAndRoot(struct node *tree);
static void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size, double norm);
static void goingDown(int index, int numberOfLeaves, double dist, double **allDist, double norm);
void leafToLeafDistance(struct node *root, double **dist, char **name, int normDistance, int branchLength);
int compNumberOfNodes(struct node *tree);

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

static int check(char **names, int index1, int index2, int index3) {
    
    // Score according to Astral-Pro tag procedure
    int score = 0;
    int aInB = 0;
    int bInA = 0;
    for (int i = index1; i < index2; i++) {
        for (int j = index2; j <= index3; j++) {
            if (strcmp(names[i], names[j]) == 0) {
                aInB++;
                break;
            }
        }
    }
    for (int i = index2; i <= index3; i++) {
        for (int j = index1; j < index2; j++) {
            if (strcmp(names[i],names[j]) == 0) {
                bInA++;
                break;
            }
        }
    }
    if (aInB == (index2 - index1)) {
        if (bInA == (index3 - index2 + 1)) {
            score = 1;
        } else {
            score = 2;
        }
    } else if (bInA == (index3 - index2 + 1)) {
        score = 2;
    } else if (aInB == 0 && bInA == 0){
        score = 0;
    } else {
        score = 3;
    }
    return score;
}


static int subScoreAndTag(struct node *current, int *index, char **names) {
    int score = 0;
    if (current->firstChild != 0) {
        struct node *this = current->firstChild;
        int index1 = *index;
        score += subScoreAndTag(this, index, names);
        struct node *next = this->nextSibling;
        int a = 0;
        while (1) {
            if (next == 0) {
                break;
            }
            *index += 1;
            int index2 = *index;
            score += subScoreAndTag(next, index, names);
            int index3 = *index;
            if (current->score == 0) {
                int checker = check(names, index1, index2, index3);
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
    score = subScoreAndTag(tree, &index, names);
    return score;
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
        current->nextSibling->score = 0;
        current->nextSibling->tag = 0;
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
            current->score = 0;
            current->tag = 0;
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
        root->score = 0;
        root->tag = 0;
        free(newRoot);
        compNumberOfLeaves(root);
    }
}


void tagAndRoot(struct node *tree) {
    struct node *current = tree;
    int bestScore = 9999999;
    int score = 999999;
    int bestId = 0;
    int id = 1;
    int maxId = compNumberOfNodes(tree);
    int size = tree->numberOfLeaves;
    char **names = (char **) calloc(sizeof(char *), size);
    for (int i = 0; i < size; i++) {
        names[i] = (char *) calloc(sizeof(char), maxNameLength);
    }
    while (id < maxId) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            if (current->idNo == id) {
                reRoot2(current, tree);
                score = scoreAndTag(tree, names);
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
                score = scoreAndTag(tree, names);
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
                scoreAndTag(tree, names);
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
            scoreAndTag(tree, names);
            break;
        }
    }
    compNumberOfLeaves(tree);
    for (int i = 0; i < size; i++) {
        free(names[i]);
    }
    free(names);
}




// Distance functions

static void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size, double norm) {
    for (int i = index - numberOfLeaves + 1; i <= index; i++) {
        for (int j = index + 1; j < size; j++) {
            allDist[i][j - i] += dist / norm ;
        }
    }
}

static void goingDown(int index, int numberOfLeaves, double dist, double **allDist, double norm) {
    for (int i = 0; i < index; i++) {
        for (int j = index; j < index + numberOfLeaves; j++) {
            allDist[i][j - i] += dist / norm;
        }
    }
}


void leafToLeafDistance(struct node *root, double **dist, char **name, int normDistance, int branchLength) {
    int size = root->numberOfLeaves;
    int index = 0;
    struct node *current = root;
    double norm = 1;
    if (normDistance == 1) {
        norm = (double) size;
    } else if (normDistance == 2) {
        norm = log((double) size);
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
        
        // ignore distances between nodes that have LCA marked as duplication, distance will result in 0 or negative and be ignored
        if (current->parent->tag == 1) {
            for (int i = index - current->numberOfLeaves + 1; i <= index; i++) {
                for (int j = index + 1; j <= index + current->nextSibling->numberOfLeaves; j++) {
                    dist[i][j - i] -= size;
                }
            }
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
