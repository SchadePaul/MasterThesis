#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "node.h"
static int check(int index, int index2, int index3, int index4, char **names);
static void score(struct node *thisNode, int index, char **names);
void astralTag(struct node *root);

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
//            printf("check at %d %d\t%d\t%d\t%d\t%d\n", child->idNo, next->idNo, index1, index2, index3, index4);
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
//    printf("A\t");
//    for (int i = index; i < index2; i++) {
//        printf("%s\t", names[i]);
//    }
//    printf("\n");
//    printf("B\t");
//    for (int i = index3; i < index4; i++) {
//        printf("%s\t", names[i]);
//    }
//    printf("\n");
//    printf("\n");
    
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
//    printf("score\t%d\n", score);
//    printf("\n");
    return score;
}
