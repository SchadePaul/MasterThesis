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
    int minNJst = 0;
    int normDistance = 0;
    double quartil = 0;
    while((c = getopt(argc,argv,"bmnq:i:o:"))!=-1) {
        switch(c) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'q':
                quartil = atof(optarg);
                break;
            case 'b':
                branchLength = 1;
                break;
            case 'm':
                minNJst = 1;
                break;
            case 'n':
                normDistance = 1;
                break;
        }
    }

    // Declare root for species tree
    struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
    njstFromFile(&speciestree, input, branchLength, minNJst, normDistance, quartil);
    if (errno != 0) {
        return errno;
    }
    saveTree(speciestree, output);
    
    freeTree(speciestree);
    
    return 0;
}

