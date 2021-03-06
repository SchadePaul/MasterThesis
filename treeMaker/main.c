#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include "methods.h"
#include "parse.h"
#include "tree.h"
#include "node.h"

int main(int argc, char **argv) {
    char *input = 0;
    char *output = 0;
    
    char isRooted = 0;
    char hasPoly = 0;
    
    char root = 0;
    char astralTag = 0;
    char notCountTag = 0;
    
    char mini = 0;
    char ustar = 0;
    
    char branchLength = 0;
    
    char norm = 0;
    char weight = 0;
    
    char average = 1;
    char median = 0;
    char mostCommon = 0;
    
    char cluster = 0;
    
    char c;
    while ((c = getopt(argc,argv,"i:o:pRr:tNmubsdn:w:a:c:")) != -1) {
        switch (c) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            
            case 'p':
                hasPoly = 1;
                break;
            case 'R':
                isRooted = 1;
                break;
                
            case 'r':
                // 1 = astral root
                // 2 = mad root
                root = atoi(optarg);
                break;
            case 't':
                astralTag = 1;
                break;
            case 'N':
                notCountTag = 1;
                break;

            case 'm':
                mini = 1;
                break;
            case 'u':
                ustar = 1;
                break;

            case 'b':
                branchLength = 1;
                break;
                
            case 'n':
                norm = atoi(optarg);
                break;
            case 'w':
                weight = atoi(optarg);
                break;
                
            case 'a':
                if (atoi(optarg) == 0) {
                    average = 1;
                } else if (atoi(optarg) == 1) {
                    average = 0;
                    median = 1;
                } else if (atoi(optarg) == 2) {
                    mostCommon = 0;
                    mostCommon = 1;
                }
                break;
            
            case 'c':
                // 0 = Neighbor Joining
                // 1 = UPGMA
                // 2 = WPGMA
                cluster = atoi(optarg);
                break;
                
            default:
                printf("unexpected input: \n%s\n", optarg);
                break;
        }
    }
    
    struct node *finalTree = 0;
    if (input != 0) {    
        makeTree(&finalTree, input, mini, ustar, norm, weight, average, median, mostCommon, cluster, root, astralTag, notCountTag, branchLength);
        if (output != 0) {
            saveTree(finalTree, output);
            freeTree(finalTree);
        }
    } 
}
