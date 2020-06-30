#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "node.h"
#include "tree.h"
#include "njst.h"

void stepByStep(const char *geneFile, const char *speciesFile);
static int compare( const void* a, const void* b);
static int compare2( const void* a, const void* b);

static int compare( const void* a, const void* b) {
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

static int compare2( const void* a, const void* b) {
     int int_a = * ( (double*) a );
     int int_b = * ( (double*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

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
    int *appearances = (int *) calloc(sizeof(int), numberOfTaxa);
    
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
        
        int *appearancesInTree = (int *) calloc(sizeof(int), numberOfTaxa);
        
        // Find taxon in taxa
        int *indexInTaxa = (int *) calloc(sizeof(int), size);
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < numberOfTaxa; k++) {
                if (strcmp(taxaInTree[j], taxa[k]) == 0) {
                    indexInTaxa[j] = k;
                    appearancesInTree[k] = 1;
                    break;
                }
            }
        }
        
        for (int j = 0; j < numberOfTaxa; j++) {
            appearances[j] += appearancesInTree[j];
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
                    if (taxaDistances[ii][jj][2 * i + 1] == 0) {
                        taxaDistances[ii][jj][2 * i + 1] = currentDistance;
                        taxaDistances[ii][jj][2 * i] = 1;
                    } else if (currentDistance < taxaDistances[ii][jj][2 * i + 1]) {
                        taxaDistances[ii][jj][2 * i + 1] = currentDistance;
                    }
                    
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
                double *allDiss = (double *) calloc(sizeof(double), numberOfTrees);
                for (int k = 0; k < numberOfTrees; k++) {
                    if (taxaDistances[i][j][2 * k] != 0) {
                        allDiss[count] = taxaDistances[i][j][2 * k + 1];
                        double currentDistance = taxaDistances[i][j][2 * k + 1];
                        distance[i][j] += currentDistance;
                        count += taxaDistances[i][j][2 * k];
                    }
                }
                if (count != 0) {
//                    qsort(allDiss, count, sizeof(double), compare2);
                    distance[i][j] = distance[i][j] / count;
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
    
    for (int i = 0; i < numberOfSpeciesTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            if (strcmp(speciesTaxa[i], taxa[j]) == 0) {
                indexToShow[i] = j;
                break;
            }
        }
    }
    
    printf("\nGene Trees Distance Matrix\n");
    
//    printf("\t");
//    for (int i = 0; i < numberOfTaxa; i++) {
//        printf("%s\t\t", taxa[i]);
//    }
    
    printf("\n\n");
    for (int i = 0; i < numberOfTaxa; i++) {
        printf("%.15s\t\t", taxa[indexToShow[i]]);
        for (int j = 0; j < numberOfTaxa; j++) {
            printf("%.2f\t", distance[indexToShow[i]][indexToShow[j]]);
        }
        printf("\n");
    }
    
    printf("\nSpeices Tree Distance Matrix\n");
    
//    printf("\t");
//    for (int i = 0; i < numberOfSpeciesTaxa; i++) {
//        printf("%s\t\t", speciesTaxa[indexToShow[i]]);
//    }
    
    double maxDist = 0;
    int *muchDist = (int *) calloc(sizeof(int), 2);
    int over = 0;
    int under = 0;
    int exact = 0;
    printf("\n\n");
    for (int i = 0; i < numberOfSpeciesTaxa; i++) {
        printf("%.15s\t\t", speciesTaxa[i]);
        for (int j = 0; j < numberOfSpeciesTaxa; j++) {
            if (i <= j) {
                int ii = i;
                int jj = j - i;
                double relDis = 0;
                if (i != j) {
                    relDis = distance[indexToShow[i]][indexToShow[j]] / dist[ii][jj];
                }
//                if (relDis < 1 && relDis != 0) {
//                    relDis = 1 / relDis;
//                }
                if (relDis > 1) {
                    printf("\033[0;31m");
                    over++;
                } else if (relDis < 1) {
                    printf("\033[0;32m");
                    under++;
                } else {
                    exact++;
                }

                if (dist[ii][jj] != distance[indexToShow[i]][indexToShow[j]]) {
                    if (relDis > maxDist) {
                        maxDist = relDis;
                        muchDist[0] = indexToShow[i];
                        muchDist[1] = indexToShow[j];
                    }
                }
                printf("%.2f\t", dist[ii][jj]);
                
            } else {
                int ii = j;
                int jj = i - j;
                double relDis = 0;
                if (indexToShow[i] != indexToShow[j]) {
                    relDis = distance[indexToShow[i]][indexToShow[j]] / dist[ii][jj];
                }
//                if (relDis < 1 && relDis != 0) {
//                    relDis = 1 / relDis;
//                }
                if (relDis > 1) {
                    printf("\033[0;31m");
                    over++;
                } else if (relDis < 1) {
                    printf("\033[0;32m");
                    under++;
                }

                if (dist[ii][jj] != distance[indexToShow[i]][indexToShow[j]]) {
                    if (relDis > maxDist) {
                        maxDist = relDis;
                        muchDist[0] = indexToShow[i];
                        muchDist[1] = indexToShow[j];
                    }
                }
                printf("%.2f\t", dist[ii][jj]);
            }
            printf("\033[0m");
        }
        printf("\n");
    }
    
    
    
//    if (muchDist[0] > muchDist[1]) {
//        int tmp = muchDist[0];
//        muchDist[0] = muchDist[1];
//        muchDist[1] = tmp;
//    }
    
    double relDis = dist[indexToShow[muchDist[0]]][indexToShow[muchDist[1]] - indexToShow[muchDist[0]]] / distance[muchDist[0]][muchDist[1]];
    if (relDis < 1 && relDis != 0) {
        relDis = 1 / relDis;
    }
    
//    printf("\nTaxa: %s\t%s\t%d\t%f\n", taxa[muchDist[0]], taxa[muchDist[1]], (int) dist[indexToShow[muchDist[0]]][indexToShow[muchDist[1]] - indexToShow[muchDist[0]]], relDis);
    
    int counter = 0;
    for (int i = 0; i < numberOfTrees; i++) {
        int ii = indexToShow[muchDist[0]];
        int jj = indexToShow[muchDist[1]] - indexToShow[muchDist[0]];
        if (jj < ii) {
            int tmp = ii;
            ii = jj;
            jj = tmp;
        }
        if (taxaDistances[muchDist[0]][muchDist[1]][2 * i] != 0) {
            counter++;
//            double relDiss = taxaDistances[muchDist[0]][muchDist[1]][2 * i + 1] / dist[ii][jj];
        }
    }
    int *diss = (int *) calloc(sizeof(int), counter);
    counter = 0;
    for (int i = 0; i < numberOfTrees; i++) {
        int ii = indexToShow[muchDist[0]];
        int jj = indexToShow[muchDist[1]] - indexToShow[muchDist[0]];
        if (jj < ii) {
            int tmp = ii;
            ii = jj;
            jj = tmp;
        }
        if (taxaDistances[muchDist[0]][muchDist[1]][2 * i] != 0) {
            diss[counter] = taxaDistances[muchDist[0]][muchDist[1]][2 * i + 1];
            counter++;
        }
    }
    qsort(diss, counter, sizeof(int), compare);
    printf("over:\t %d, under:\t %d, exact:\t%d\n", over, under, exact);
    
    printf("number of Taxa: %d \n", numberOfTaxa);
    for (int i = 0; i < numberOfTaxa; i++) {
        printf("%s\t%d\n", taxa[indexToShow[i]], appearances[indexToShow[i]]);
    }
//    for (int i = 0; i < counter; i++) {
//        printf("%d\n", diss[i]);
//    }
    
}
