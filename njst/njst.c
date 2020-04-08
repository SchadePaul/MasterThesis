#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

void averageDistances(double ***allDistances, double **distance, int numberOfNames, int numberOfTrees) {
    for (int i = 0; i < numberOfNames; i++) {
        for (int j = 0; j < numberOfNames; j++) {
            int counter = 0;
            for (int k = 0; k < numberOfTrees; k++) {
                distance[i][j] += allDistances[i][j][k];
                if (allDistances[i][j][k] != 0) {
                    counter++;
                }
            }
            if (counter != 0) {
                distance[i][j] = distance[i][j] / (double) counter;
            } else {
                distance[i][j] = 0;
            }
        }
    }
}

void njstFromFile(struct node **root, const char *filename) {
    
    // read file to array of chars
    char **newickTree;
    int numberOfTrees;
    readFileToArray(filename, &newickTree, &numberOfTrees);
    
    if (errno != 0) {
        return;
    }
    
    
    // array of chars to trees
    char **allLeafNames;
    int numberOfLeafNames;
    struct node **tree = (struct node **) calloc(sizeof(struct node *), numberOfTrees);
    for (int i = 0; i < numberOfTrees; i++) {
        // parse newick format to tree, add leaf names to allLeafNamesArray
        newickTreeToTree(newickTree[i], &tree[i], &allLeafNames, &numberOfLeafNames);
        free(newickTree[i]);
    }
    free(newickTree);
    
    // allocate space for allLeafDistances array
    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
        for (int j = 0; j < numberOfLeafNames; j++) {
            allLeafDistances[i][j] = (double *) calloc(sizeof(double), numberOfTrees);
        }
    }
    
    // for all trees copmute distances between leaf nodes
    for (int i = 0; i < numberOfTrees; i++) {
        // get number of leafs
        int size = tree[i]->numberOfLeaves;
        double **dist = (double **) calloc(sizeof(double*), size);
        for (int j = 0; j < size; j++) {
            dist[j] = (double *) calloc(sizeof(double), size);
        }
        char **name = (char **) calloc(sizeof(char *), size);
        for (int j = 0; j < size; j++) {
            name[j] = (char *) calloc(sizeof(char), maxNameLength);
        }
        
        // compute leaf to leaf distance for all leafs
        leafToLeafDistance(tree[i], dist, size, name);
        freeTree(tree[i]);
        
        // copy distances to allLeafDistances array
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int ii;
                int jj;
                for (int l = 0; l < numberOfLeafNames; l++) {
                    if (strcmp(name[j], allLeafNames[l]) == 0) {
                        ii = l;
                    }
                    if (strcmp(name[k], allLeafNames[l]) == 0) {
                        jj = l;
                    }
                }
                allLeafDistances[ii][jj][i] = dist[j][k];
            }
        }
        
        for (int j = 0; j < size; j++) {
            free(name[j]);
            free(dist[j]);
        }
        free(dist);
        free(name);
    }
    
    free(tree);
    
    double **distance = (double **) calloc(sizeof(double*), numberOfLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        distance[i] = (double *) calloc(sizeof(double), numberOfLeafNames);
    }
    
    // Average over all trees
    averageDistances(allLeafDistances, distance, numberOfLeafNames, numberOfTrees);
    
    // Free allocate space allDistances
    for (int i = 0; i < numberOfLeafNames; i++) {
        for (int j = 0; j < numberOfLeafNames; j++) {
            free(allLeafDistances[i][j]);
        }
        free(allLeafDistances[i]);
    }
    free(allLeafDistances);
    
    // copmute species tree from distance array
    makeTreeFromDistanceArray(distance, numberOfLeafNames, root, allLeafNames);
    
    // Free allocate space numberOfNames, distance
    for (int i = 0; i < numberOfLeafNames; i++) {
        free(allLeafNames[i]);
    }
    free(allLeafNames);
}
