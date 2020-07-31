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
    
    while((c = getopt(argc,argv,"bmw:n:o:i:"))!=-1) {
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
            case 'w':
                weight = atoi(optarg);
                break;
            case 'm':
                mini = 1;
                break;
            case 'b':
                branchLength = 1;
                break;
            default:
                printf("unexpected input: \n%s\n", optarg);
                break;
        }
    }
    
    printf("input:%s\toutput:%s\tbranchLength:%d\tmini:%d\tnorm:%d\tweight:%d\n", input, output, branchLength, mini, norm, weight);
    struct node *speciestree = (struct node*) calloc(sizeof(struct node), 1);
    printf("cha\n");
    calcSpeciesTree(&speciestree, input, branchLength, mini, norm, weight);
    saveTree(speciestree, output);
    freeTree(speciestree);
    return 0;
}

