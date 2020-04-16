#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "src/njst.h"
#include "src/parse.h"
#include "src/tree.h"
#include "src/node.h"
#include <unistd.h>
#include <getopt.h>
#include <errno.h>

int main(int argc, char **argv) {
    char c;
    char *input;
    char *output;
    while((c = getopt(argc,argv,"i:o:"))!=-1) {
        switch(c) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
        }
    }

    // Declare root for species tree
    struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
    njstFromFile(&speciestree, input);
    if (errno != 0) {
        return errno;
    }
    saveTree(speciestree, output);
    
    freeTree(speciestree);
    
    return 0;
}

