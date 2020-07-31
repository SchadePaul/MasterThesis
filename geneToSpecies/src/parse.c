#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "node.h"
#include <errno.h>

const int lineWidthIntervall = 50;
const char endOfTree = ';';

int compNumberOfLeaves(struct node *current);
void readFileToTrees(const char *filename, struct node ***trees_ptr, char ***taxa_ptr, int *numberOfTrees, int *numberOfTaxa);
static void addTree(struct node ***trees_ptr, int *numberOfTrees, struct node *root);
static void addValue(char *value, struct node *current);
static void addName(char ***taxa_ptr, int *numberOfTaxa, char *name);
static int isTreeOperator(char c);
static int maxDepth(struct node *current);
static void fillArray(char *array, struct node *current, int width, int offsetX, int offsetY, int length);
void printTree(struct node *tree, int length);
void saveTree(struct node *tree, const char *name);

static void addTree(struct node ***trees_ptr, int *numberOfTrees, struct node *root) {
    struct node **trees = (struct node **) calloc(sizeof(struct node *), *numberOfTrees);
    for (int i = 0; i < *numberOfTrees - 1; i++) {
        trees[i] = trees_ptr[0][i];
    }
    trees[*numberOfTrees - 1] = root;
    if (*numberOfTrees > 1) {
        free(*trees_ptr);
    }
    trees_ptr[0] = trees;
}

static void addValue(char *value, struct node *current) {
    current->distToParent = atof(value);
}

static void addName(char ***taxa_ptr, int *numberOfTaxa, char *name) {
    
    for (int i = 0; i < *numberOfTaxa; i++) {
        if (strcmp(taxa_ptr[0][i], name) == 0) {
            return;
        }
    }
    printf("add Name:\t%s\n", name);
    char **taxa_new = (char **) calloc(sizeof(char *), *numberOfTaxa + 1);
    for (int i = 0; i < *numberOfTaxa; i++) {
        taxa_new[i] = taxa_ptr[0][i];
    }
    
    
    if (*numberOfTaxa > 0) {
        free(taxa_ptr[0]);
    }
    taxa_ptr[0] = taxa_new;
    taxa_new[*numberOfTaxa] = (char *) calloc(sizeof(char), strlen(name) + 1);
    strcpy(taxa_new[*numberOfTaxa], name);
    *numberOfTaxa += 1;
}

void readFileToTrees(const char *filename, struct node ***trees_ptr, char ***taxa_ptr, int *numberOfTrees, int *numberOfTaxa) {
    
    FILE *f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        return;
    }
    
    char c;
    int newTree = 1;
    struct node *current;
    int posInName = 0;
    char *name = (char *) calloc(sizeof(char), maxNameLength);
    int posInValue = 0;
    char *value = (char *) calloc(sizeof(char), maxValueLength);
    int isLeafName = 1;
    
    while (fread(&c, sizeof(c), 1, f) == 1) {
        // skip linebreaks
        if (c == '\n') {
            continue;
        }
        
        if (newTree == 1) {
            // add new root
            *numberOfTrees += 1;
            current = (struct node *) calloc(sizeof(struct node), 1);
            addTree(trees_ptr, numberOfTrees, current);
            newTree = 0;
        }
        
        if (c != endOfTree) {
            if (isTreeOperator(c)) {
                switch (c) {
                    case '(':
                        isLeafName = 1;
                        posInName = 1;
                        posInValue = 0;
                        current->firstChild = (struct node*) calloc(sizeof(struct node), 1);
                        current->firstChild->parent = current;
                        current = current->firstChild;
                        break;
                    case ')':
                        isLeafName = 0;
                        posInName = 1;
                        value[posInValue - 1] = 0;
                        addValue(value, current);
                        posInValue = 0;
                        current = current->parent;
                        break;
                    case ',':
                        isLeafName = 1;
                        posInName = 1;
                        value[posInValue - 1] = 0;
                        addValue(value, current);
                        posInValue = 0;
                        current->nextSibling = (struct node*) calloc(sizeof(struct node), 1);
                        current->nextSibling->parent = current->parent;
                        current = current->nextSibling;
                        break;
                    case ':':
                        name[posInName - 1] = 0;
                        if (isLeafName) {
                            addName(taxa_ptr, numberOfTaxa, name);
                        }
                        posInName = 0;
                        posInValue = 1;
                        break;
                    default:
                        break;
                }
            } else {
                if (posInName > 0) {
                    name[posInName - 1] = c;
                    posInName++;
                } else if (posInValue > 0) {
                    value[posInValue - 1] = c;
                    posInValue++;
                }
            }
        } else {
            newTree = 1;
            compNumberOfLeaves(current);
        }
    }
    free(name);
    free(value);
    fclose(f);
}

int compNumberOfLeaves(struct node *current) {
    
    current->numberOfLeaves = 0;
    current->numberOfChildren = 0;
    current->numberOfNodes = 1;
    
    if (current->firstChild != NULL) {
        struct node *workOn = current->firstChild;
        current->numberOfLeaves += compNumberOfLeaves(workOn);
        current->numberOfChildren += 1;
        current->numberOfNodes += workOn->numberOfNodes;
        while (workOn->nextSibling != NULL) {
            workOn = workOn->nextSibling;
            current->numberOfLeaves += compNumberOfLeaves(workOn);
            current->numberOfChildren += 1;
            current->numberOfNodes += workOn->numberOfNodes;
        }
    } else {
        current->numberOfLeaves = 1;
    }
    return current->numberOfLeaves;
}

static int isTreeOperator(char c) {
    if (c == '(' || c == ')' || c == ',' || c == ':') {
        return 1;
    }
    return 0;
}

static int maxDepth(struct node *current) {
    int max = 1;
    if (current->firstChild != NULL) {
        max = maxDepth(current->firstChild) + 1;
    }
    if (current->nextSibling != NULL) {
        int check = maxDepth(current->nextSibling);
        if (check > max) {
            max = check;
        }
    }
    return max;
}

static void fillArray(char *array, struct node *current, int width, int offsetX, int offsetY, int length) {
    for (int i = 0; i < length; i++) {
        array[offsetY * width * length + offsetX * length + i] = current->name[i];
    }
    
    if (current->firstChild != NULL) {
        array[(offsetY + 1) * width * length + offsetX * length] = '|';
        fillArray(array, current->firstChild, width, offsetX, offsetY + 2, length);
    }
    if (current->nextSibling != NULL) {
        int index = (offsetY - 1) * width * length + (offsetX + current->numberOfLeaves) * length;
        array[index] = '\\';
        fillArray(array, current->nextSibling, width, offsetX + current->numberOfLeaves, offsetY, length);
        index--;
        while (array[index - width * length] == 0) {
            array[index - width * length] = '-';
            if (index % (width * length) == 0) {
                break;
            }
            index--;

        }
    }
}

void printTree(struct node *tree, int length) {
    int depth = maxDepth(tree);
    int terminalNodes = compNumberOfLeaves(tree);
    char *array = calloc(sizeof(char), depth * terminalNodes * length * 2 - 1);
    fillArray(array, tree, terminalNodes, 0, 0, length);
    for (int i = 0; i < depth * 2; i++) {
        for (int j = 0; j < terminalNodes * length; j++) {
            char c = array[i * terminalNodes * length + j];
            if (c != 0) {
                printf("%c", c);
            } else {
                printf("%c", ' ');
            }
        }
        printf("\n");
    }
    free(array);
}




void saveTree(struct node *tree, const char *name) {
    FILE *f = fopen(name, "w");
    
    if (!f) {
        fprintf(stderr, "Error saving file: %s\n", strerror(errno));
        return;
    }
    
    struct node *current = tree;
    int goingUp = 0;
    while (!(current == tree && goingUp == 1)) {
        if (current->firstChild != NULL && goingUp == 0) {
            fprintf(f, "(");
            goingUp = 0;
            current = current->firstChild;
        } else if (current->nextSibling != NULL) {
            if (current->name[0] != placeholderName) {
                fprintf(f, "%s:%f,", current->name, current->distToParent);
            } else {
                fprintf(f, ":%f,", current->distToParent);
            }
            goingUp = 0;
            current = current->nextSibling;
        } else {
            if (current->name[0] != placeholderName) {
                fprintf(f, "%s:%f)", current->name, current->distToParent);
            } else {
                fprintf(f, ":%f)", current->distToParent);
            }
            goingUp = 1;
            current = current->parent;
        }
    }
    if (current->name[0] != placeholderName) {
        fprintf(f, "%s;", current->name);
    } else {
        fprintf(f, ";");
    }
    
    fclose(f);
}
