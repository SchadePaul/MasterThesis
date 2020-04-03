#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "njst.h"
#include "parse.h"
#include "tree.h"
#include "node.h"
#include <unistd.h>

int main(int argc, char **argv) {
    // Read filename
    char *filename = NULL;
    filename = (char *) calloc(sizeof(char), strlen(argv[1]) + 1);
    strncpy(filename, argv[1], strlen(argv[1]));

    // Declare root for species tree
    struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
    njstFromFile(&speciestree, filename);
    printTree(speciestree);
    free(filename);
    filename = (char *) calloc(sizeof(char), strlen(argv[2] + 1));
    strncpy(filename, argv[2], strlen(argv[2]));
    saveTree(speciestree, filename);
    free(filename);
    freeTree(speciestree);
    
    return 0;
}

