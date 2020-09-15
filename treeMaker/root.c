#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
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
static void nodeToNodeDist(struct node *current, double **dist, int index, int size);

static void nodeToNodeDist(struct node *current, double **dist, int index, int size) {
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
            nodeToNodeDist(child, dist, index, size);
            index += child->numberOfNodes;
            child = child->nextSibling;
        }
    }
}

void mad(struct node **root) {
    int size = (*root)->numberOfNodes;
    double **nodeToNodeDistances = (double **) calloc((size_t) (size - 1), sizeof(double));
    for (int i = 0; i < size - 1; i++) {
        nodeToNodeDistances[i] = (double *) calloc((size_t) (size - 1 - i), sizeof(double));
    }
    
    nodeToNodeDist((*root), nodeToNodeDistances, 0, size);
    
    for (int i = 0; i < size - 1; i++) {
        free(nodeToNodeDistances[i]);
    }
    free(nodeToNodeDistances);
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
