#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "node.h"

const int lineWidthIntervall = 50;
const char endOfTree = ';';
const char placeholderName = '*';

int compNumberOfTerminalNodes(struct node *current) {
    int num = 0;
    if (current->firstChild != NULL) {
        num = compNumberOfTerminalNodes(current->firstChild);
    } else {
        num = 1;
    }
    current->numberOfTerminalNodes = num;
    if (current->nextSibling != NULL) {
        num += compNumberOfTerminalNodes(current->nextSibling);
    }
    return num;
}

void extendCharArraybyTree(char ***newickTree, int newNumberOfTrees) {
    char **newArray = calloc(sizeof(char *), newNumberOfTrees);
    for (int i = 0; i < newNumberOfTrees - 1; i++) {
        newArray[i] = (*newickTree)[i];
    }
    free(*newickTree);
    *newickTree = newArray;
}

void extenCharSize(char **oldChar, int fromSize, int toSize) {
    char *newChar = calloc(sizeof(char), toSize + 1);
    for (int i = 0; i < fromSize; i++) {
        newChar[i] = (*oldChar)[i];
    }
    free(*oldChar);
    *oldChar = newChar;
}

void readFileToArray(const char *filename, char ***newickTree, int *numberOfTrees) {
    FILE *f = fopen(filename, "rb");
    char c;
    *numberOfTrees = 1;
    int indexTree = 0;
    int posInLine = 0;
    
    // Set char array for first tree;
    *newickTree = (char **) calloc(sizeof(char *), 1);
    
    // Read whole file
    while (fread(&c, sizeof(c), 1, f) == 1) {
        if (c == '\n') {
            continue;
        }
        if (c != endOfTree) {
            if (posInLine % lineWidthIntervall == 0) {
                int fromSize = posInLine;
                int toSize = fromSize + lineWidthIntervall;
                extenCharSize(&(*newickTree)[indexTree], fromSize, toSize);
            }
            (*newickTree)[indexTree][posInLine] = c;
            posInLine++;
        } else {
            indexTree++;
            extendCharArraybyTree(newickTree, indexTree + 1);
            posInLine = 0;
        }
    }
    *numberOfTrees = indexTree;
    fclose(f);
}

void addName(char *name, char ***allNames, int *currentLength) {
    int exists = 0;
    for (int i = 0; i < *currentLength; i++) {
        if (strcmp((*allNames)[i], name) == 0) {
            exists = 1;
            break;
        }
    }
    if (!exists) {
        char **new = (char **) calloc(sizeof(char*), *currentLength + 1);
        for (int i = 0; i < *currentLength; i++) {
            new[i] = (*allNames)[i];
        }
        new[*currentLength] = (char*) calloc(sizeof(char), strlen(name));
        strcpy(new[*currentLength], name);
        free(*allNames);
        *allNames = new;
        *currentLength += 1;
    }
}

int isChar(char c) {
    return (c == '_' || isalpha(c));
}

int isValue(char c) {
    if (isdigit(c) || c == '.') {
        return 1;
    }
    return 0;
}

int isTreeOperator(char c) {
    if (c == '(' || c == ')' || c == ',' || c == ':') {
        return 1;
    }
    return 0;
}

void newickTreeToTree(char *newickTree, struct node **tree, char ***allNames, int *numberOfNames) {

    int posInNewickTree = 0;
    int posInCurrentName = 0;
    int posInCurrentValue = 0;
    (*tree) = (struct node *) calloc(sizeof(struct node), 1);
    struct node *current = (*tree);
    char *value = calloc(sizeof(char), maxValueLength);
    
    while (newickTree[posInNewickTree] != 0) {
        char c = newickTree[posInNewickTree];
        if (isChar(c)) {
            current->name[posInCurrentName] = c;
            posInCurrentName++;
        } else if (isValue(c)) {
            value[posInCurrentValue] = c;
            posInCurrentValue++;
        } else if (isTreeOperator(c)) {
            if (posInCurrentValue > 0) {
                value[posInCurrentValue] = 0;
                current->distToParent = atof(value);
                posInCurrentValue = 0;
            }
            if (posInCurrentName > 0) {
                posInCurrentName = 0;
                addName(current->name, allNames, numberOfNames);
            }
            if (current->name[0] == 0) {
                current->name[0] = placeholderName;
            }
            if (c == '(') {
                current->firstChild = (struct node*) calloc(sizeof(struct node), 1);
                current->firstChild->parent = current;
                current = current->firstChild;
            } else if (c == ',') {
                current->nextSibling = (struct node*) calloc(sizeof(struct node), 1);
                current->nextSibling->parent = current->parent;
                current = current->nextSibling;
            } else if (c == ')') {
                current = current->parent;
            }
        }
        posInNewickTree++;
    }
    free(value);
    compNumberOfTerminalNodes(*tree);
    
}

int maxDepth(struct node *current) {
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

void fillArray(char *array, struct node *current, int width, int offsetX, int offsetY) {
    for (int i = 0; i < maxNameLength; i++) {
        array[offsetY * width * maxNameLength + offsetX * maxNameLength + i] = current->name[i];
    }
    
    if (current->firstChild != NULL) {
        array[(offsetY + 1) * width * maxNameLength + offsetX * maxNameLength] = '|';
        fillArray(array, current->firstChild, width, offsetX, offsetY + 2);
    }
    if (current->nextSibling != NULL) {
        int index = (offsetY - 1) * width * maxNameLength + (offsetX + current->numberOfTerminalNodes) * maxNameLength;
        array[index] = '\\';
        fillArray(array, current->nextSibling, width, offsetX + current->numberOfTerminalNodes, offsetY);
        index--;
        while (array[index - width * maxNameLength] == 0) {
            array[index - width * maxNameLength] = '-';
            if (index % (width * maxNameLength) == 0) {
                break;
            }
            index--;

        }
    }
}

void printTree(struct node *tree) {
    int depth = maxDepth(tree);
    int terminalNodes = compNumberOfTerminalNodes(tree);
    char *array = calloc(sizeof(char), depth * terminalNodes * maxNameLength * 2 - 1);
    fillArray(array, tree, terminalNodes, 0, 0);
    for (int i = 0; i < depth * 2; i++) {
        for (int j = 0; j < terminalNodes * maxNameLength; j++) {
            char c = array[i * terminalNodes * maxNameLength + j];
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
    FILE *f;
    char *rootname = (char *) calloc(sizeof(char), maxNameLength);
    strcpy(rootname, tree->name);
    f = fopen(name, "w");
    int index = 0;
    struct node *current = tree;
    int size = compNumberOfTerminalNodes(current);
    while (index < size - 1) {
        while (current->firstChild != NULL) {
            current = current->firstChild;
            if (current->name[0] != placeholderName) {
                fprintf(f, "(%s:%f", current->name, current->distToParent);
            } else {
                fprintf(f, "(:%f", current->distToParent);
            }
            
        }
        while (current->nextSibling == NULL) {
            current = current->parent;
            fprintf(f, ")");
        }
        current = current->nextSibling;
        if (current->name[0] != placeholderName) {
            fprintf(f, ",%s:%f", current->name, current->distToParent);
        } else {
            fprintf(f, ",:%f", current->distToParent);
        }
        
        index++;
    }
    if (rootname[0] != placeholderName) {
        fprintf(f, ")%s;", rootname);
    } else{
        fprintf(f, ");");
    }
    free(rootname);
    fclose(f);
    
}
