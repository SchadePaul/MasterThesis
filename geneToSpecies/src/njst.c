#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <float.h>

void calcSpeciesTree(struct node **speciesTree, const char *filename, int branchLength, int mini, int norm, int weight);
static void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames);


void calcSpeciesTree(struct node **speciesTree, const char *filename, int branchLength, int mini, int norm, int weight) {
    
    int numberOfTrees = 0;
    char **taxa;
    int numberOfTaxa;
    struct node **trees;
    
    readFileToTrees(&trees, filename, &numberOfTrees, &taxa, &numberOfTaxa);
    
    double ***taxaDistances = (double ***) calloc(sizeof(double **), (size_t) numberOfTaxa);
    for (int i = 0; i < numberOfTaxa; i++) {
        taxaDistances[i] = (double **) calloc(sizeof(double *), (size_t) numberOfTaxa);
        for (int j = 0; j < numberOfTaxa; j++) {
            taxaDistances[i][j] = (double *) calloc(sizeof(double), (size_t) numberOfTrees * 2);
        }
    }
    
    for (int i = 0; i < numberOfTrees; i++) {
        // Number of Leaves in this tree
        int size = trees[i]->numberOfLeaves;
                        
        // Half-Matrix for distances between leafs
        double **dist = (double **) calloc(sizeof(double *), size);
        for (int j = 0; j < size; j++) {
            dist[j] = (double *) calloc(sizeof(double), size - j);
        }
                
        // Array of taxa in tree
        char **taxaInTree = (char **) calloc(sizeof(char *), size);
        for (int j = 0; j < size; j++) {
            taxaInTree[j] = (char *) calloc(sizeof(char), maxNameLength);
        }
            
        // Compute leaf distances
        leafToLeafDistance(trees[i], dist, taxaInTree, branchLength);
            
        freeTree(trees[i]);
                
        // Find taxon in taxa
        int *indexInTaxa = (int *) calloc(sizeof(int), size);
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < numberOfTaxa; k++) {
                if (strcmp(taxaInTree[j], taxa[k]) == 0) {
                    indexInTaxa[j] = k;
                    break;
                }
            }
        }
        
        double weightFactor = 1;
        double normFactor = 1;
        
        switch (weight) {
            case 1:
                weightFactor = size;
                break;
                
            default:
                break;
        }
        
        switch (norm) {
            case 1:
                normFactor = size;
                break;
                
            default:
                break;
        }
        
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size - j; k++) {
                if (dist[j][k] > 0) {
                    int ii = indexInTaxa[j];
                    int jj = indexInTaxa[k + j];
                    if (jj < ii) {
                        int tmp = ii;
                        ii = jj;
                        jj = tmp;
                    }
                    double currentDistance = (dist[j][k] * weightFactor) * normFactor;
                    if (mini) {
                        if (taxaDistances[ii][jj][2 * i + 1] == 0 || currentDistance < taxaDistances[ii][jj][2 * i + 1]) {
                            taxaDistances[ii][jj][2 * i + 1] = currentDistance;
                            taxaDistances[ii][jj][2 * i] = weightFactor;
                        }
                    } else {
                        taxaDistances[ii][jj][2 * i + 1] += currentDistance;
                        taxaDistances[ii][jj][2 * i] += weightFactor;
                    }
                }
            }
        }
    
        for (int j = 0; j < size; j++) {
            free(taxaInTree[j]);
            free(dist[j]);
        }
        free(indexInTaxa);
        free(taxaInTree);
        free(dist);
    }
    free(trees);
        
    double **distance = (double **) calloc(sizeof(double*), numberOfTaxa);
    for (int i = 0; i < numberOfTaxa; i++) {
        distance[i] = (double *) calloc(sizeof(double), numberOfTaxa);
    }

    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            if (i < j) {
                int count = 0;
                for (int k = 0; k < numberOfTrees; k++) {
                    if (taxaDistances[i][j][2 * k] != 0) {
                        double currentDistance = taxaDistances[i][j][2 * k + 1];
                        distance[i][j] += currentDistance;
                        count += taxaDistances[i][j][2 * k];
                    }
                }
                if (count != 0) {
                    distance[i][j] = distance[i][j] / count;
                }
            } else {
                distance[i][j] = distance[j][i];
            }
        }
    }
    
    makeTreeFromDistanceArray(distance, numberOfTaxa, speciesTree, taxa);
    
        
    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            free(taxaDistances[i][j]);
        }
        free(distance[i]);
        free(taxaDistances[i]);
        free(taxa[i]);
    }
    free(distance);
    free(taxaDistances);
    free(taxa);
    
}




static void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames) {
    // read file to array of chars
    char **newickTree;
    readFileToArray(filename, &newickTree, numberOfTrees);
    
    if (errno != 0) {
        return;
    }
    
    // array of chars to trees
    *trees = (struct node **) calloc(sizeof(struct node *), *numberOfTrees);
    for (int i = 0; i < *numberOfTrees; i++) {
        // parse newick format to tree, add leaf names to allLeafNamesArray
        newickTreeToTree(newickTree[i], &((*trees)[i]), allLeafNames, numberOfLeafNames);
        free(newickTree[i]);
    }
    free(newickTree);
    
}
