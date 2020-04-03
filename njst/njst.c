#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    
    // array of chars to trees
    char **allNames;
    int numberOfNames;
    struct node **tree = (struct node **) calloc(sizeof(struct node *), numberOfTrees);
    for (int i = 0; i < numberOfTrees; i++) {
        newickTreeToTree(newickTree[i], &tree[i], &allNames, &numberOfNames);
        free(newickTree[i]);
    }
    free(newickTree);
    double ***allDistances = (double ***) calloc(sizeof(double **), numberOfNames);
    for (int i = 0; i < numberOfNames; i++) {
        allDistances[i] = (double **) calloc(sizeof(double *), numberOfNames);
        for (int j = 0; j < numberOfNames; j++) {
            allDistances[i][j] = (double *) calloc(sizeof(double), numberOfTrees);
        }
    }
    for (int i = 0; i < numberOfTrees; i++) {
        int size = tree[i]->numberOfTerminalNodes;
        double **dist = (double **) calloc(sizeof(double*), size);
        for (int j = 0; j < size; j++) {
            dist[j] = (double *) calloc(sizeof(double), size);
        }
        char **name = (char **) calloc(sizeof(char *), size);
        for (int j = 0; j < size; j++) {
            name[j] = (char *) calloc(sizeof(char), maxNameLength);
        }
        allToAllDistance(tree[i], dist, size, name);
        freeTree(tree[i]);
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                int ii;
                int jj;
                for (int l = 0; l < numberOfNames; l++) {
                    if (strcmp(name[j], allNames[l]) == 0) {
                        ii = l;
                    }
                    if (strcmp(name[k], allNames[l]) == 0) {
                        jj = l;
                    }
                }
                allDistances[ii][jj][i] = dist[j][k];
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
    
    double **distance = (double **) calloc(sizeof(double*), numberOfNames);
    for (int i = 0; i < numberOfNames; i++) {
        distance[i] = (double *) calloc(sizeof(double), numberOfNames);
    }
    
    averageDistances(allDistances, distance, numberOfNames, numberOfTrees);
    
    // Free allocate space allDistances
    for (int i = 0; i < numberOfNames; i++) {
        for (int j = 0; j < numberOfNames; j++) {
            free(allDistances[i][j]);
        }
        free(allDistances[i]);
    }
    free(allDistances);
    
    makeTreeFromDistanceMatrix(distance, numberOfNames, root, allNames);
    
    // Free allocate space numberOfNames, distance
    for (int i = 0; i < numberOfNames; i++) {
        free(allNames[i]);
        free(distance[i]);
    }
    free(distance);
    free(allNames);
    
}
