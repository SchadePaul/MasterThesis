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
    int miniNJ = 0;
    int norm = 0;
    int weighted = 0;
    int tagged = 0;
    int root = 0;
    int square = 0;
    int doustar = 0;
    while((c = getopt(argc,argv,"ubmwtrsn:o:i:"))!=-1) {
        switch(c) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'n':
                norm = (int) atof(optarg);
                break;
            case 's':
                square = 1;
                break;
            case 'r':
                root = 1;
                break;
            case 't':
                tagged = 1;
                break;
            case 'w':
                weighted = 1;
                break;
            case 'm':
                miniNJ = 1;
                break;
            case 'b':
                branchLength = 1;
                break;
            case 'u':
                doustar = 1;
                break;
        }
    }

    // Declare root for species tree
    struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
    
    if (doustar) {
        ustar(&speciestree, input, tagged, root);
    } else {
        if (tagged) {
            taggedNJFromFile(&speciestree, input, norm, root, weighted, square);
        } else {
            if (miniNJ) {
                miniNJFromFile(&speciestree, input, norm, weighted, branchLength);
            } else {
                njstFromFile(&speciestree, input, norm, weighted, branchLength, square);
            }
        }
    }
    
    
    
    if (errno != 0) {
        return errno;
    }
    saveTree(speciestree, output);
    
    freeTree(speciestree);
    
    return 0;
}

