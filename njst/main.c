#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "src/njst.h"
#include "src/parse.h"
#include "src/tree.h"
#include "src/node.h"
#include "src/compare.h"
#include <unistd.h>
#include <getopt.h>
#include <errno.h>

int main(int argc, char **argv) {
    char c;
    char *input;
    char *output;
    char *compare = 0;
    int branchLength = 0;
    int mini = 0;
    int norm = 0;
    int weight = 0;
    int tag = 0;
    int root = 0;
    int square = 0;
    int ustar = 0;
    double miniPairs = 0;
    while((c = getopt(argc,argv,"ubmtsc:r:w:p:n:o:i:"))!=-1) {
        switch(c) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'c':
                compare = optarg;
                break;
            case 'n':
                norm = atoi(optarg);
                break;
            case 'p':
                miniPairs = atof(optarg);
            case 's':
                square = 1;
                break;
            case 'r':
                root = atoi(optarg);
                break;
            case 't':
                tag = 1;
                break;
            case 'w':
                weight = atoi(optarg);
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
            default:
                printf("unexpected input: \n%s\n", optarg);
                break;
        }
    }

    if (compare != 0) {
        stepByStep(input, compare);
    } else {
        // Declare root for species tree
        struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
        
        inferSpeciesTreeFromGeneTrees(&speciestree, input, mini, ustar, norm, branchLength, weight, square, tag, root, miniPairs);
        
        if (errno != 0) {
            return errno;
        }
        saveTree(speciestree, output);
        
        freeTree(speciestree);
    }
    
    return 0;
}

