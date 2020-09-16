#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "parse.h"
#include "node.h"
#include "tree.h"
#include "cluster.h"
#include "root.h"

const int arrayExtends = 64;

void makeTree(struct node **finalTree, const char *input, char mini, char ustar, char norm, char weight, char average, char median, char mostCommon, char cluster, char root, char astralTag, char notCountTag, char branchLength);
int compare_dbl (const void *a, const void *b);

int compare_dbl (const void *a, const void *b) {
    if (*(double *)a > *(double *)b) {
        return 1;
    } else if (*(double *)a < *(double *)b) {
        return -1;
    } else {
        return 0;
    }
}



void makeTree(struct node **finalTree, const char *input, char mini, char ustar, char norm, char weight, char average, char median, char mostCommon, char cluster, char toRoot, char astralTag, char notCountTag, char branchLength) {
    
    int numberOfTrees = 0;
    struct node **trees = 0;
    
    int numberOfSpecies = 0;
    char **speciesNames = 0;
    

    readFileToTrees(&trees, input, &numberOfTrees, &speciesNames, &numberOfSpecies);
    
    // distance matrix and
    // indices where per tree distances per taxa start
    double ***distances = (double ***) calloc(sizeof(double **), numberOfSpecies - 1);
    int ***distancesIndices = (int ***) calloc(sizeof(int **), numberOfSpecies - 1);
    for (int i = 0; i < numberOfSpecies - 1; i++) {
        distances[i] = (double **) calloc(sizeof(double *), numberOfSpecies - 1 - i);
        distancesIndices[i] = (int **) calloc(sizeof(int *), numberOfSpecies - 1 - i);
        for (int j = 0; j < numberOfSpecies - 1 - i; j++) {
            distances[i][j] = (double *) calloc(sizeof(double), arrayExtends);
            distancesIndices[i][j] = (int *) calloc(sizeof(int), numberOfTrees + 1);
        }
    }
    
    
    for (int treeNumber = 0; treeNumber < numberOfTrees; treeNumber++) {
        
        if (toRoot == 1) {
            astralRoot(&(trees[treeNumber]));
        } else if (toRoot == 2) {
            mad(&(trees[treeNumber]));
        }

        int treeSize = trees[treeNumber]->numberOfLeaves;
        
        char **treeSpeciesNames = (char **) calloc(sizeof(char *), treeSize);
        for (int i = 0; i < treeSize; i++) {
            treeSpeciesNames[i] = (char *) calloc(sizeof(char), maxNameLength);
        }
        double **treeDistances = (double **) calloc(sizeof(double *), treeSize - 1);
        for (int i = 0; i < treeSize - 1; i++) {
            treeDistances[i] = (double *) calloc(sizeof(double), treeSize - 1 - i);
        }
        
        leafToLeafDistance(trees[treeNumber], treeDistances, treeSpeciesNames, branchLength, astralTag, notCountTag);
        
        int *speciesIndices = (int *) calloc(sizeof(int), treeSize);
        for (int i = 0; i < treeSize; i++) {
            for (int j = 0; j < numberOfSpecies; j++) {
                if (strcmp(treeSpeciesNames[i], speciesNames[j]) == 0) {
                    speciesIndices[i] = j;
                    break;
                }
            }
        }
        
        double normFactor = 1.0;
        if (norm == 1) {
            normFactor = (double) treeSize;
        } else if (norm == 2) {
            normFactor = log((double) treeSize) / log(2.0);
        }
        
        int weightFactor = 0;
        if (weight == 0) {
            weightFactor = 1;
        } else if (weight == 1) {
            weightFactor = treeSize;
        } else if (weight == 2) {
            for (int i = 0; i < treeSize; i++) {
                int add = 1;
                for (int j = i + 1; j < treeSize; j++) {
                    if (speciesIndices[i] == speciesIndices[j]) {
                        add = 0;
                        break;
                    }
                }
                weightFactor += add;
            }
        }
        
        for (int i = 0; i < numberOfSpecies - 1; i++) {
            for (int j = 0; j < numberOfSpecies - 1 - i; j++) {
                distancesIndices[i][j][treeNumber + 1] = distancesIndices[i][j][treeNumber];
            }
        }
        
        // copy treeDistances to distances
        for (int i = 0; i < treeSize - 1; i++) {
            for (int j = 0; j < treeSize - 1 - i; j++) {
                if (treeDistances[i][j] <= 0) {
                    continue;
                }
                int ii = speciesIndices[i];
                int jj = speciesIndices[j + i + 1];
                
                if (ii == jj) {
                    continue;
                }
                if (jj < ii) {
                    int tmp = ii;
                    ii = jj;
                    jj = tmp;
                }
                
                jj = jj - 1 - ii;
                
                if (mini == 1) {
                    int index = distancesIndices[ii][jj][treeNumber];
                    if ((distancesIndices[ii][jj][treeNumber + 1] > index) && (weightFactor * treeDistances[i][j] / normFactor < distancesIndices[ii][jj][index])) {
//                        for (int repeat = 0; repeat < weightFactor; repeat++) {
                        if (index != 0 && (index % arrayExtends == 0)) {
                            double *newArray = (double *) calloc(sizeof(double), (index/arrayExtends + 1) * arrayExtends);
                            for (int copyI = 0; copyI < index; copyI++) {
                                newArray[copyI] = distances[ii][jj][copyI];
                            }
                            free(distances[ii][jj]);
                            distances[ii][jj] = newArray;
                        }
                        distances[ii][jj][index] = weightFactor * treeDistances[i][j] / normFactor;
                        index += 1;
                        if (index != 0 && (index % arrayExtends == 0)) {
                            double *newArray = (double *) calloc(sizeof(double), (index/arrayExtends + 1) * arrayExtends);
                            for (int copyI = 0; copyI < index; copyI++) {
                                newArray[copyI] = distances[ii][jj][copyI];
                            }
                            free(distances[ii][jj]);
                            distances[ii][jj] = newArray;
                        }
                        distances[ii][jj][index] = weightFactor;
//                        }
                    } else if (distancesIndices[ii][jj][treeNumber + 1] == index) {
//                        for (int repeat = 0; repeat < weightFactor; repeat++) {
                        int index = distancesIndices[ii][jj][treeNumber + 1];
                        distancesIndices[ii][jj][treeNumber + 1] += 2;
                        if (index != 0 && (index % arrayExtends == 0)) {
                            double *newArray = (double *) calloc(sizeof(double), (index/arrayExtends + 1) * arrayExtends);
                            for (int copyI = 0; copyI < index; copyI++) {
                                newArray[copyI] = distances[ii][jj][copyI];
                            }
                            free(distances[ii][jj]);
                            distances[ii][jj] = newArray;
                        }
                        distances[ii][jj][index] = weightFactor * treeDistances[i][j] / normFactor;
                        index += 1;
                        if (index != 0 && (index % arrayExtends == 0)) {
                            double *newArray = (double *) calloc(sizeof(double), (index/arrayExtends + 1) * arrayExtends);
                            for (int copyI = 0; copyI < index; copyI++) {
                                newArray[copyI] = distances[ii][jj][copyI];
                            }
                            free(distances[ii][jj]);
                            distances[ii][jj] = newArray;
                        }
                        distances[ii][jj][index] = weightFactor;
//                        }
                    }
                } else {
//                    for (int repeat = 0; repeat < weightFactor; repeat++) {
                    int index = distancesIndices[ii][jj][treeNumber + 1];
                    distancesIndices[ii][jj][treeNumber + 1] += 2;
                    if (index != 0 && (index % arrayExtends == 0)) {
                        double *newArray = (double *) calloc(sizeof(double), (index/arrayExtends + 1) * arrayExtends);
                        for (int copyI = 0; copyI < index; copyI++) {
                            newArray[copyI] = distances[ii][jj][copyI];
                        }
                        free(distances[ii][jj]);
                        distances[ii][jj] = newArray;
                    }
                    distances[ii][jj][index] = weightFactor * treeDistances[i][j] / normFactor;
                    index += 1;
                    if (index != 0 && (index % arrayExtends == 0)) {
                        double *newArray = (double *) calloc(sizeof(double), (index/arrayExtends + 1) * arrayExtends);
                        for (int copyI = 0; copyI < index; copyI++) {
                            newArray[copyI] = distances[ii][jj][copyI];
                        }
                        free(distances[ii][jj]);
                        distances[ii][jj] = newArray;
                    }
                    distances[ii][jj][index] = weightFactor;
//                    }
                }
            }
        }
        
        for (int i = 0; i < treeSize; i++) {
            free(treeSpeciesNames[i]);
        }
        if (treeSpeciesNames != 0) {
            free(treeSpeciesNames);
        }
        for (int i = 0; i < treeSize - 1; i++) {
            free(treeDistances[i]);
        }
        if (treeDistances != 0) {
            free(treeDistances);
        }
        if (speciesIndices != 0) {
            free(speciesIndices);
        }
        freeTree(trees[treeNumber]);
    }
    
    double **finalDistance = (double **) calloc(sizeof(double *), numberOfSpecies);
    for (int i = 0; i < numberOfSpecies; i++) {
        finalDistance[i] = (double *) calloc(sizeof(double), numberOfSpecies);
    }
    
    if (median != 0) {
        if (weight != 0) {
            printf("median and weighting not compatible\n");
        }
        // changed from adding multiple time to add only once with weight because much faster
        for (int i = 0; i < numberOfSpecies - 1; i++) {
            for (int j = 0; j < numberOfSpecies - 1 - i; j++) {
                qsort(distances[i][j], (size_t) distancesIndices[i][j][numberOfTrees], sizeof(double), compare_dbl);
                int k = 0;
                while (distances[i][j][k] <= 0) {
                    k++;
                }
                finalDistance[i][j + 1 + i] = distances[i][j][(int) (k + ((distancesIndices[i][j][numberOfTrees] - k) / 2))];
                finalDistance[j + 1 + i][i] = distances[i][j][(int) (k + ((distancesIndices[i][j][numberOfTrees] - k) / 2))];
            }
        }
    } else if (mostCommon != 0) {
        if (weight != 0) {
            printf("weighting and mostCommon not compatible\n");
        }
        // changed from adding multiple time to add only once with weight because much faster
        for (int i = 0; i < numberOfSpecies - 1; i++) {
            for (int j = 0; j < numberOfSpecies - 1 - i; j++) {
                int max = 0;
                // find max
                for (int k = 0; k < distancesIndices[i][j][numberOfTrees]; k++) {
                    if ((int) distances[i][j][k] > max) {
                        max = (int) distances[i][j][k];
                    }
                }
                int *numberOfAppearances = (int *) calloc((size_t) (max + 1), sizeof(int));
                for (int k = 0; k < distancesIndices[i][j][numberOfTrees]; k++) {
                    int index = (int) distances[i][j][k];
                    if (index > 0) {
                        numberOfAppearances[index] += 1;
                    }
                }
                int maxAppearance = 0;
                int indexMaxAppearance = 0;
                for (int k = 0; k < max; k++) {
                    if (numberOfAppearances[k] > maxAppearance) {
                        maxAppearance = numberOfAppearances[k];
                        indexMaxAppearance = k;
                    }
                }
                if (numberOfAppearances != 0) {
                    free(numberOfAppearances);
                }
        
                finalDistance[i][j + 1 + i] = indexMaxAppearance;
                finalDistance[j + 1 + i][i] = indexMaxAppearance;
            }
        }
    } else if (average != 0) {
        for (int i = 0; i < numberOfSpecies - 1; i++) {
            for (int j = 0; j < numberOfSpecies - 1 - i; j++) {
                double thisDist = 0.0;
                double numberOfDist = 0.0;
                if (ustar == 0) {
                    for (int k = 0; k < distancesIndices[i][j][numberOfTrees]; k += 2) {
                        thisDist += distances[i][j][k];
                        numberOfDist += distances[i][j][k];
                    }
                } else if (ustar != 0) {
                    for (int k = 0; k < numberOfTrees; k++) {
                        double thisTreeDist = 0.0;
                        double numberOfTreeDist = 0.0;
                        for (int l = distancesIndices[i][j][k]; l < distancesIndices[i][j][k + 1]; l += 2) {
                            thisTreeDist += distances[i][j][l];
                            numberOfTreeDist += distances[i][j][l + 1];
                        }
                        if (numberOfTreeDist > 0) {
                            thisDist += thisTreeDist / numberOfTreeDist;
                            numberOfDist += 1.0;
                        }
                    }
                }
                finalDistance[i][j + 1 + i] = thisDist / numberOfDist;
                finalDistance[j + 1 + i][i] = thisDist / numberOfDist;
            }
        }
    }
    
    struct node *root = (struct node *) calloc(1, sizeof(struct node));

    if (cluster == 0) {
        //NJ
        NJ(&root, speciesNames, finalDistance, numberOfSpecies);
    } else if (cluster == 1) {
        //UPGMA
        UWPGMA(&root, speciesNames, finalDistance, numberOfSpecies, 0);
    } else if (cluster == 2) {
        //WPGMA
        UWPGMA(&root, speciesNames, finalDistance, numberOfSpecies, 1);
    }

//    printTree(root, 4);
    (*finalTree) = root;
    
//    printf("\t");
//    for (int i = 0; i < numberOfSpecies; i++) {
//        printf("%s\t", speciesNames[i]);
//    }
//    printf("\n");
//    for (int i = 0; i < numberOfSpecies; i++) {
//        printf("%s\t", speciesNames[i]);
//        for (int j = 0; j < numberOfSpecies; j++) {
//            printf("%.3f\t", finalDistance[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
    for (int i = 0; i < numberOfSpecies; i++) {
        free(finalDistance[i]);
    }
    if (finalDistance != 0) {
        free(finalDistance);
    }
    for (int i = 0; i < numberOfSpecies - 1; i++) {
        for (int j = 0; j < numberOfSpecies - 1 - i; j++) {
            free(distancesIndices[i][j]);
            free(distances[i][j]);
        }
        free(distancesIndices[i]);
        free(distances[i]);
    }
    if (distances != 0) {
        free(distances);
    }
    if (distancesIndices != 0) {
        free(distancesIndices);
    }
    for (int i = 0; i < numberOfSpecies; i++) {
        free(speciesNames[i]);
    }
    if (speciesNames != 0) {
        free(speciesNames);
    }
    if (trees != 0) {
        free(trees);
    }
}
