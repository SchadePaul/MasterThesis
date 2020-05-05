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
    int branchLength = 0;
    int mini = 0;
    int norm = 0;
    int weight = 0;
    int tag = 0;
    int root = 0;
    int square = 0;
    int ustar = 0;
    while((c = getopt(argc,argv,"ubmwtrsn:o:i:"))!=-1) {
        switch(c) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'n':
                norm = atoi(optarg);
                break;
            case 's':
                square = 1;
                break;
            case 'r':
                root = 1;
                break;
            case 't':
                tag = 1;
                break;
            case 'w':
                weight = 1;
                break;
            case 'm':
                mini = 1;
                break;
            case 'b':
                branchLength = 1;
                break;
            case 'u':
                ustar = 1;
                break;
        }
    }

    // Declare root for species tree
    struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
    
    inferSpeciesTreeFromGeneTrees(&speciestree, input, mini, ustar, norm, branchLength, weight, square, tag, root);
    
    if (errno != 0) {
        return errno;
    }
    saveTree(speciestree, output);
    
    freeTree(speciestree);
    
    return 0;
}

