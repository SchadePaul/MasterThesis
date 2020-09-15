#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "node.h"
#include <errno.h>

const int lineWidthIntervall = 1024;
const int readChars = 128;
const char endOfTree = ';';

void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames);
void readFileToArray(const char *filename, char ***newickTree, int *numberOfTrees);
static void extenCharSize(char **oldChar, int fromSize, int toSize);
static void extendCharArraybyTree(char ***newickTree, int newNumberOfTrees);
void newickTreeToTree(char *newickTree, struct node **tree, char ***allLeafNames, int *numberOfLeafNames);
static int isTreeOperator(char c);
static void addName(char *name, char ***allNames, int *currentLength);
int compNumberOfLeaves(struct node *current);
static int maxDepth(struct node *current);
static void fillArray(char *array, struct node *current, int width, int offsetX, int offsetY, int length);
void printTree(struct node *tree, int length);
void saveTree(struct node *tree, const char *name);



void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames) {
    // read file to array of chars
    char **newickTree;
    readFileToArray(filename, &newickTree, numberOfTrees);
    
    free(newickTree[*numberOfTrees]);
    if (errno != 0) {
        return;
    }
    
    // array of chars to trees
    *trees = (struct node **) calloc((size_t) *numberOfTrees, sizeof(struct node *));
    for (int i = 0; i < *numberOfTrees; i++) {
        // parse newick format to tree, add leaf names to allLeafNamesArray
        newickTreeToTree(newickTree[i], &((*trees)[i]), allLeafNames, numberOfLeafNames);
        free(newickTree[i]);
        
    }
    free(newickTree);
}

void readFileToArray(const char *filename, char ***newickTree, int *numberOfTrees) {
    FILE *f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        return;
    }
    
    
    char c[readChars];
    *numberOfTrees = 0;
    int indexTree = 0;
    int posInLine = 0;
    
    // Set char array for first tree;
    *newickTree = (char **) calloc(sizeof(char *), 1);
    // Read whole file
    int read = 0;
    
    do {
        read = fread(&c, sizeof(char), readChars, f);
        for (int i = 0; i < read; i++) {
            // Ignore newline char
            if (c[i] == '\n') {
                continue;
            }

            if (c[i] != endOfTree) {
                if (posInLine % lineWidthIntervall == 0) {
                    int fromSize = posInLine;
                    int toSize = fromSize + lineWidthIntervall;
                    extenCharSize(&((*newickTree)[indexTree]), fromSize, toSize);
                }
                
                // Add char to array
                (*newickTree)[indexTree][posInLine] = c[i];
                posInLine++;
            } else {
                indexTree++;
                extendCharArraybyTree(newickTree, indexTree + 1);
                posInLine = 0;
            }
        }
    } while (read == readChars);
    
//    while (fread(&c, sizeof(char), readChars, f) > 0) {
//
//        // Ignore newline char
//        if (c == '\n') {
//            continue;
//        }
//
//        if (c != endOfTree) {
//            if (posInLine % lineWidthIntervall == 0) {
//                int fromSize = posInLine;
//                int toSize = fromSize + lineWidthIntervall;
//                extenCharSize(&((*newickTree)[indexTree]), fromSize, toSize);
//            }
//
//            // Add char to array
//            (*newickTree)[indexTree][posInLine] = c;
//            posInLine++;
//        } else {
//            indexTree++;
//            extendCharArraybyTree(newickTree, indexTree + 1);
//            posInLine = 0;
//        }
//    }
    *numberOfTrees = indexTree;
    fclose(f);
}

static void extenCharSize(char **oldChar, int fromSize, int toSize) {
    // extend to size + 1 (to have a terminal 0)
    char *newChar = (char *) calloc(sizeof(char), toSize + 1);
    for (int i = 0; i < fromSize; i++) {
        newChar[i] = (*oldChar)[i];
    }
    free(*oldChar);
    *oldChar = newChar;
}

static void extendCharArraybyTree(char ***newickTree, int newNumberOfTrees) {
    char **newArray = (char **) calloc((size_t) newNumberOfTrees, sizeof(char *));
    for (int i = 0; i < newNumberOfTrees - 1; i++) {
        newArray[i] = (*newickTree)[i];
    }
    free(*newickTree);
    *newickTree = newArray;
}

void newickTreeToTree(char *newickTree, struct node **tree, char ***allLeafNames, int *numberOfLeafNames) {

    int posInNewickTree = 0;
    int posInCurrentName = 0;
    int posInCurrentValue = 0;
    int readingValue = 0;
    (*tree) = (struct node *) calloc(sizeof(struct node), 1);
    struct node *current = (*tree);
    int id = 0;
    current->idNo = id;
    id++;
    char *value = (char *) calloc((size_t) maxValueLength +1 , sizeof(char));
    
    // Only add names of leaves
    int dontAddNextName = 0;
    while (newickTree[posInNewickTree] != 0) {
        char c = newickTree[posInNewickTree];
        if (isTreeOperator(c)) {
            readingValue = 0;
            if (posInCurrentValue > 0) {
                value[posInCurrentValue] = 0;
                current->distToParent = atof(value);
                posInCurrentValue = 0;
            }
            if (posInCurrentName > 0) {
                posInCurrentName = 0;
                if (dontAddNextName == 0) {
                    addName(current->name, allLeafNames, numberOfLeafNames);
                } else {
                    dontAddNextName = 0;
                }
            }
            if (current->name[0] == 0) {
                current->name[0] = placeholderName;
                dontAddNextName = 0;
            }
            if (c == '(') {
                current->firstChild = (struct node*) calloc(sizeof(struct node), 1);
                current->firstChild->parent = current;
                current = current->firstChild;
                current->idNo = id;
                id++;
            } else if (c == ',') {
                current->nextSibling = (struct node*) calloc(sizeof(struct node), 1);
                current->nextSibling->parent = current->parent;
                current = current->nextSibling;
                current->idNo = id;
                id++;
                dontAddNextName = 0;
            } else if (c == ')') {
                current = current->parent;
                dontAddNextName = 1;
            } else if (c == ':') {
                readingValue = 1;
            }
        } else if (readingValue == 1) {
            value[posInCurrentValue] = c;
            posInCurrentValue++;
        } else {
            current->name[posInCurrentName] = c;
            posInCurrentName++;
        }
        posInNewickTree++;
    }
    free(value);
    compNumberOfLeaves(*tree);
}

static int isTreeOperator(char c) {
    if (c == '(' || c == ')' || c == ',' || c == ':') {
        return 1;
    }
    return 0;
}

static void addName(char *name, char ***allNames, int *currentLength) {
    int exists = 0;
    for (int i = 0; i < *currentLength; i++) {
        if (strcmp((*allNames)[i], name) == 0) {
            exists = 1;
            break;
        }
    }
    if (!exists) {
        char **new = (char **) calloc((size_t) (*currentLength + 1), sizeof(char *));
        for (int i = 0; i < *currentLength; i++) {
            new[i] = (*allNames)[i];
        }
        new[*currentLength] = (char *) calloc((size_t) strlen(name) + 1, sizeof(char));
        strcpy(new[*currentLength], name);
        if (*currentLength > 0) {
            free(*allNames);
        }
        *allNames = new;
        *currentLength += 1;
    }
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
//    current->tag = 0;
    return current->numberOfLeaves;
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
//    array[offsetY * width * length + offsetX * length + (length - 3)] = (current->score / 10) + '0';
//    array[offsetY * width * length + offsetX * length + (length - 2)] = (current->score % 10) + '0';
    
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
