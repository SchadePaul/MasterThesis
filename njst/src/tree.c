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
int scoreAndTag(struct node *tree);
//static void reRoot(struct node *atNode, struct node **best, int *bestScore);
void tagAndRoot(struct node *tree);
static void goingUp(int index, int numberOfLeaves, double dist, double **allDist, int size, double norm);
static void goingDown(int index, int numberOfLeaves, double dist, double **allDist, double norm);
void leafToLeafDistance(struct node *root, double **dist, char **name, int normDistance, int branchLength);

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
//
//static void copySubtree(struct node *from, struct node *at) {
//    //Copy subtree of new roots first child
//    struct node *oldCurrent = from;
//    struct node *newCurrent = at;
//    while (1) {
//        while (oldCurrent->firstChild != NULL) {
//            oldCurrent = oldCurrent->firstChild;
//            newCurrent->firstChild = (struct node *) calloc(sizeof(struct node), 1);
//            newCurrent->firstChild->parent = newCurrent;
//            newCurrent = newCurrent->firstChild;
//            strcpy(newCurrent->name, oldCurrent->name);
//            newCurrent->tag = oldCurrent->tag;
//            newCurrent->score = oldCurrent->score;
//            newCurrent->idNo = oldCurrent->idNo;
//        }
//        while (oldCurrent->nextSibling == NULL) {
//            if (oldCurrent == from) {
//                break;
//            }
//            oldCurrent = oldCurrent->parent;
//            newCurrent = newCurrent->parent;
//        }
//        if (oldCurrent == from) {
//            break;
//        }
//        oldCurrent = oldCurrent->nextSibling;
//        newCurrent->nextSibling = (struct node *) calloc(sizeof(struct node), 1);
//        newCurrent->nextSibling->parent = newCurrent->parent;
//        newCurrent = newCurrent->nextSibling;
//        strcpy(newCurrent->name, oldCurrent->name);
//        newCurrent->tag = oldCurrent->tag;
//        newCurrent->score = oldCurrent->score;
//        newCurrent->idNo = oldCurrent->idNo;
//    }
//}

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
    
    // Subroutine of tag procedure from Astral-Pro
    int score = 0;
    if (current->firstChild != NULL) {
        int index1 = *index;
        score = subScoreAndTag(current->firstChild, index, names);
        *index += 1;
        int index2 = *index;
        score += subScoreAndTag(current->firstChild->nextSibling, index, names);
        int index3 = *index;
        if (current->score != 0) {
            if (current->score == -1) {
                score = 0;
            } else {
                score = current->score;
            }
        } else {
            int checker = check(names, index1, index2, index3);
            score += checker;
            if (checker != 0) {
                current->tag = 1;
                current->score = score;
            }
            else {
                current->tag = 0;
                current->score = -1;
            }
        }
    } else {
        strcpy(names[*index], current->name);
    }
    return score;
}

int scoreAndTag(struct node *tree) {
    
    // Tag procedre from Astral-Pro
    int index = 0;
    int score = 0;
    int size = tree->numberOfLeaves;
    char **names = (char **) calloc(sizeof(char *), size);
    for (int i = 0; i < size; i++) {
        names[i] = (char *) calloc(sizeof(char), maxNameLength);
    }
    score = subScoreAndTag(tree, &index, names);
    for (int i = 0; i < size; i++) {
        free(names[i]);
    }
    free(names);
    return score;
}
//
//static void reRoot(struct node *atNode, struct node **best, int *bestScore) {
//
//    struct node *oldCurrent = atNode;
//
//    // No need to reroot of childs of root, because it would result in the same tree
//    if (oldCurrent->parent->parent == NULL) {
//        return;
//    }
//
//    // Set new root
//    struct node *newRoot = (struct node *) calloc(sizeof(struct node), 1);
//    struct node *current = newRoot;
//    current->name[0] = placeholderName;
//    current->firstChild = (struct node *) calloc(sizeof(struct node), 1);
//    current->firstChild->parent = current;
//    current = current->firstChild;
//    strcpy(current->name, oldCurrent->name);
//
//    copySubtree(oldCurrent, current);
//
//    // travers old tree
//    while (oldCurrent->parent != NULL) {
//        current->nextSibling = (struct node *) calloc(sizeof(struct node), 1);
//        current->nextSibling->parent = current->parent;
//        current = current->nextSibling;
//
//        // Decide if going up from first or second child
//        if (oldCurrent->nextSibling != NULL) {
//            // Skip old root
//            if (oldCurrent->parent->parent == NULL) {
//                oldCurrent = oldCurrent->nextSibling;
//            } else {
//                oldCurrent = oldCurrent->parent;
//                strcpy(current->name, oldCurrent->name);
//
//                oldCurrent = oldCurrent->firstChild->nextSibling;
//                current->firstChild = (struct node *) calloc(sizeof(struct node), 1);
//                current->firstChild->parent = current;
//                current = current->firstChild;
//            }
//        } else {
//            // Skip old root
//            if (oldCurrent->parent->parent == NULL) {
//                oldCurrent = oldCurrent->parent->firstChild;
//            } else {
//                oldCurrent = oldCurrent->parent;
//                strcpy(current->name, oldCurrent->name);
//
//
//                oldCurrent = oldCurrent->firstChild;
//                current->firstChild = (struct node *) calloc(sizeof(struct node), 1);
//                current->firstChild->parent = current;
//                current = current->firstChild;
//            }
//        }
//
//        // Copy old subtree to new tree
//        strcpy(current->name, oldCurrent->name);
//
//        copySubtree(oldCurrent, current);
//        oldCurrent = oldCurrent->parent;
//    }
//
//    compNumberOfLeaves(newRoot);
//
//    // Get score and tag new tree
//    int thisScore = scoreAndTag(newRoot);
//    if (thisScore < *bestScore) {
//        *bestScore = thisScore;
//        (*best) = newRoot;
//    } else {
//        freeTree(newRoot);
//    }
//}

static void reRoot2(struct node *tree, struct node *root) {
    if (tree->parent == root) {
        return;
    }
    struct node *newRoot = (struct node *) calloc(sizeof(struct node), 1);
    newRoot->firstChild = tree;
    struct node *toParent = newRoot;
    struct node *current = tree;
    
    if (current->parent->firstChild == current) {
        current->parent->firstChild = current->nextSibling;
    }
    current->nextSibling = current->parent;
    current->parent->score = 0;
    current->parent->tag = 0;
    current->parent = toParent;
    current = current->nextSibling;
    
    while (current->parent->parent != 0) {
        current->firstChild->nextSibling = current->parent;
                current->parent->score = 0;
        current->parent->tag = 0;
        if (current->parent->firstChild == current) {
            current->parent->firstChild = current->nextSibling;
            current->nextSibling = 0;
        }
        current->parent = toParent;
        toParent = current;
        current = current->firstChild->nextSibling;
    }
    
    if (current->parent->firstChild == current) {
        current->firstChild->nextSibling = current->nextSibling;
        current->nextSibling->parent = current;
        current->nextSibling = 0;
        current->parent = toParent;
    } else {
        current->firstChild->nextSibling = current->parent->firstChild;
        current->firstChild->nextSibling->parent = current;
        current->nextSibling = 0;
        current->parent->firstChild->nextSibling = 0;
        current->parent = toParent;
    }
    
    
    root->firstChild = newRoot->firstChild;
    root->firstChild->parent = root;
    root->firstChild->nextSibling->parent = root;
    free(newRoot);
    compNumberOfLeaves(root);
    root->score = 0;
}


void tagAndRoot(struct node *tree) {
    
    struct node *current = tree;
//    printTree(tree);
    int bestScore = scoreAndTag(tree);
    int bestId = 0;
    int id = 1;
    int maxId = tree->numberOfLeaves * 2 - 1;
    while (id < maxId) {
        while (current->firstChild != 0) {
            current = current->firstChild;
            if (current->idNo == id) {
                if (current->parent->parent != 0) {
                    reRoot2(current, tree);
                    int score = scoreAndTag(tree);
                    if (score < bestScore) {
                        bestScore = score;
                        bestId = id;
                    }
                }
                id++;
            }
        }
        while (current->nextSibling == 0) {
            current = current->parent;
            if (current->parent == 0) {
                break;
            }
        }
        if (current->parent == 0) {
            
        } else {
            current = current->nextSibling;
            if (current->idNo == id) {
                if (current->parent->parent != 0) {
                    reRoot2(current, tree);
                    int score = scoreAndTag(tree);
                    if (score < bestScore) {
                        bestScore = score;
                        bestId = id;
                    }
                }
                id++;
            }
        }
    }
    current = tree;
    if (bestId == 0) {
        bestId = 1;
    }
//    printf("bestID: %d\t%d\n", bestId, bestScore);
//    printf("Tag tree of size\t%d\n", tree->numberOfLeaves);
    if (current->score != bestScore) {
        while (1) {
                while (current->firstChild != 0) {
                    current = current->firstChild;
                    if (current->idNo == bestId) {
                        reRoot2(current, tree);
        //                scoreAndTag(tree);
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
        //            scoreAndTag(tree);
                    break;
                }
            }
    }
//    printf("Tag tree of size\t%d\tdone\n", tree->numberOfLeaves);
}







//void tagAndRoot(struct node **tree) {
//    // current position in original tree
//    struct node *current = *tree;
//
//    // pointer to tree with lowest
//    struct node *best = *tree;
//    int bestScore = scoreAndTag(*tree);
//    // if score is zero no need to reroot
//    if (bestScore == 0) {
//        return;
//    }
//    // Travers tree and reRoot at every node
//    while (1) {
//        while (current->firstChild != NULL) {
//            current = current->firstChild;
//            reRoot(current, &best, &bestScore);
//        }
//        while (current->nextSibling == NULL) {
//            current = current->parent;
//            if (current == *tree) {
//                break;
//            }
//        }
//        if (current == *tree) {
//            break;
//        }
//        current = current->nextSibling;
//        reRoot(current, &best, &bestScore);
//    }
//
//    // Set pointer to best tree
//    *tree = best;
//}

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
