#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"
#include "tree.h"
#include "njst.h"

void stepByStep(const char *geneFile, const char *speciesFile);

void stepByStep(const char *geneFile, const char *speciesFile) {
    
    int numberOfTrees = 0;
    char **taxa;
    int numberOfTaxa = 0;
    struct node **trees;
    
    readFileToTrees(&trees, geneFile, &numberOfTrees, &taxa, &numberOfTaxa);
    
    int numberOfSpeciesTree = 0;
    char **speciesTaxa;
    int numberOfSpeciesTaxa = 0;
    struct node **speciesTree;
    
    readFileToTrees(&speciesTree, speciesFile, &numberOfSpeciesTree, &speciesTaxa, &numberOfSpeciesTaxa);
    
    // Array of taxa distances per tree
    double ***taxaDistances = (double ***) calloc(sizeof(double **), (size_t) numberOfTaxa);
    for (int i = 0; i < numberOfTaxa; i++) {
        taxaDistances[i] = (double **) calloc(sizeof(double *), (size_t) numberOfTaxa);
        for (int j = 0; j < numberOfTaxa; j++) {
            taxaDistances[i][j] = (double *) calloc(sizeof(double), (size_t) numberOfTrees * 2);
        }
    }
    
    for (int i = 0; i < numberOfTrees; i++) {
        int size = trees[i]->numberOfLeaves;
        
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
        leafToLeafDistance(trees[i], dist, taxaInTree, 0);
        
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
                    double currentDistance = dist[j][k];
                    taxaDistances[ii][jj][2 * i + 1] += currentDistance;
                    taxaDistances[ii][jj][2 * i] += 1;
                }
            }
        }
        
        for (int j = 0; j < size; j++) {
            free(taxaInTree[j]);
            free(dist[j]);
        }
        free(taxaInTree);
        free(dist);
    }
    
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
                    if (count != 0) {
                        distance[i][j] = distance[i][j] / count;
                    }
                }
            } else {
                distance[i][j] = distance[j][i];
            }
        }
    }
    
    int speciesSize = speciesTree[0]->numberOfLeaves;
    
    double **dist = (double **) calloc(sizeof(double *), speciesSize);
    for (int j = 0; j < speciesSize; j++) {
        dist[j] = (double *) calloc(sizeof(double), speciesSize - j);
    }
    
    // Array of taxa in tree
    char **taxaInSpeciesTree = (char **) calloc(sizeof(char *), speciesSize);
    for (int j = 0; j < speciesSize; j++) {
        taxaInSpeciesTree[j] = (char *) calloc(sizeof(char), maxNameLength);
    }
    
    // Compute leaf distances
    leafToLeafDistance(speciesTree[0], dist, taxaInSpeciesTree, 0);
    
    
    int *indexToShow = (int *) calloc(sizeof(int), numberOfSpeciesTaxa);
    
    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfSpeciesTaxa; j++) {
            if (strcmp(taxa[i], speciesTaxa[j]) == 0) {
                indexToShow[i] = j;
                break;
            }
        }
    }
    
    printf("\nGene Trees Distance Matrix\n");
    
    printf("\t");
    for (int i = 0; i < numberOfTaxa; i++) {
        printf("%s\t\t", taxa[i]);
    }
    
    printf("\n\n");
    for (int i = 0; i < numberOfTaxa; i++) {
        printf("%s\t", taxa[i]);
        for (int j = 0; j < numberOfTaxa; j++) {
            printf("%f\t", distance[i][j]);
        }
        printf("\n");
    }
    
    printf("\nSpeices Tree Distance Matrix\n");
    
    printf("\t");
    for (int i = 0; i < numberOfSpeciesTaxa; i++) {
        printf("%s\t\t", speciesTaxa[indexToShow[i]]);
    }
    
    double maxDist = 0;
    char **muchDist = (char **) calloc(sizeof(char *), 2);
    for (int i = 0; i < 2; i++) {
        muchDist[i] = (char *) calloc(sizeof(char), maxNameLength);
    }
    
    printf("\n\n");
    for (int i = 0; i < numberOfSpeciesTaxa; i++) {
        printf("%s\t", taxa[i]);
        for (int j = 0; j < numberOfSpeciesTaxa; j++) {
//            printf("compare: \t%s\t%s\n", taxa[i], taxa[j]);
            if (indexToShow[i] <= indexToShow[j]) {
                
                double relDis = 0;
                if (i != j) {
                    relDis = dist[indexToShow[i]][indexToShow[j] - indexToShow[i]] / distance[i][j];
                }
                if (relDis < 1 && relDis != 0) {
                    relDis = 1 / relDis;
                }
                if (relDis > 1.15) {
                    printf("\033[0;31m");
                }
                
                if (dist[indexToShow[i]][indexToShow[j] - indexToShow[i]] != distance[i][j]) {
//                    printf("difference\t%f\t%f\n", dist[indexToShow[i]][indexToShow[j] - indexToShow[i]], distance[i][j]);
                    
                    if (relDis > maxDist) {
                        maxDist = relDis;
                        muchDist[0] = taxa[i];
                        muchDist[1] = taxa[j];
                        
                    }
//                    printf("\t\t%f\n", relDis);
                    
                }
//                printf("%f\t", relDis);
                printf("%f\t", dist[indexToShow[i]][indexToShow[j] - indexToShow[i]]);
                printf("\033[0m");
            } else {
                printf("%f\t", dist[indexToShow[j]][indexToShow[i] - indexToShow[j]]);
            }
        }
        printf("\n");
    }
    
    printf("\nTaxa: %s\t%s\n", muchDist[0], muchDist[1]);
    
}
