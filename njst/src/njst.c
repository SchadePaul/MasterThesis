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

void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames);
void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root, double miniPairs, double quartil, int closeFriends);
static int compare( const void* a, const void* b);
void test(struct node **speciesTree, const char *filename);

static int compare( const void* a, const void* b) {
     int int_a = * ( (double*) a );
     int int_b = * ( (double*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

void test(struct node **speciesTree, const char *filename) {
    int numberOfTrees = 0;
    char **taxa;
    int numberOfTaxa = 0;
    struct node **trees;
        
    readFileToTrees(&trees, filename, &numberOfTrees, &taxa, &numberOfTaxa);
        
        
    // Array of taxa distances per tree
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
        leafToLeafDistance(trees[i], dist, taxaInTree, 0);
        
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
                    if (taxaDistances[ii][jj][2 * i + 1] == 0 || currentDistance < taxaDistances[ii][jj][2 * i + 1]) {
                        taxaDistances[ii][jj][2 * i + 1] = currentDistance;
                        taxaDistances[ii][jj][2 * i] = 1;
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
    free(trees);
    
    double **distance = (double **) calloc(sizeof(double*), numberOfTaxa);
    for (int i = 0; i < numberOfTaxa; i++) {
        distance[i] = (double *) calloc(sizeof(double), numberOfTaxa);
    }

//    printf("farDistances\n");
    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            if (i < j) {
                int count = 0;
                for (int k = 0; k < numberOfTrees; k++) {
                    if (taxaDistances[i][j][2 * k] != 0) {
                        double currentDistance = taxaDistances[i][j][2 * k + 1];
                        distance[i][j] += currentDistance;
                        count += taxaDistances[i][j][2 * k];
//                        if ((i == 5 && j == 13) || (i == 8 && j == 34)) {
//                            if (currentDistance > 7) {
//                                if (i == 5) {
//                                    printf("K L \t%d\t %.2f\n", k, currentDistance);
//                                } else {
//                                    printf("I H \t%d\t %.2f\n", k, currentDistance);
//                                }
//
//                            }
//                        }
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
    
//    printf("\n\n");
//
//    printf("\t");
//    for (int i = 0; i < numberOfTaxa; i++) {
//        printf("%s\t", taxa[i]);
//    }
//    printf("\n\n");
//
//    for (int i = 0; i < numberOfTaxa; i++) {
//        printf("%s\t", taxa[i]);
//        double min = DBL_MAX;
//        for (int j = 0; j < numberOfTaxa; j++) {
//            if (distance[i][j] < min && distance[i][j] != 0) {
//                min = distance[i][j];
//            }
//        }
//        for (int j = 0; j < numberOfTaxa; j++) {
//            if (distance[i][j] == min) {
//                printf("\033[1m\033[31m");
//            }
//            printf("%.2f\t", distance[i][j]);
//            printf("\033[0m");
//        }
//        printf("\n");
//    }
//    printf("\n");
    
    makeTree2(distance, numberOfTaxa, speciesTree, taxa);
//    makeTreeFromDistanceArray(distance, numberOfTaxa, speciesTree, taxa);
        
//    printTree(*speciesTree, 4);
    
    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            free(taxaDistances[i][j]);
        }
        free(taxaDistances[i]);
        free(taxa[i]);
    }
    free(taxaDistances);
    free(taxa);

}

void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root, double miniPairs, double quartil, int closeFriends) {
    // Read file and get important data (#Trees, #Taxa, Taxa)
    int numberOfTrees = 0;
    char **taxa;
    int numberOfTaxa = 0;
    struct node **trees;
    
    readFileToTrees(&trees, filename, &numberOfTrees, &taxa, &numberOfTaxa);
    
    
    // Array of taxa distances per tree
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
        
        // Prework on the trees
        if (root == 1) {
            tagAndRoot(trees[i]);
        } else if (root == 2) {
            mad(&trees[i]);
        }
        
        if (tag) {
            scoreAndTag2(trees[i]);
        }
        
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
        int *indexx = calloc(sizeof(int), 1);
        deletedTaggedDistance(trees[i], dist, indexx);
        
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
        
        int **allIndices = (int **) calloc(sizeof(int *), size);
        int *indexInMinDiances = (int *) calloc(sizeof(int), size);
        int numberOfDifferentTaxa = 0;
        for (int j = 0; j < size; j++) {
            allIndices[j] = (int *) calloc(sizeof(int), 2);
            for (int k = 0; k <= j; k++) {
                if (allIndices[k][0] == indexInTaxa[j]) {
                    if (allIndices[k][1] == 0) {
                        numberOfDifferentTaxa += 1;
                    }
                    allIndices[k][1] += 1;
                    indexInMinDiances[j] = k;
                    break;
                } else if (allIndices[k][0] == 0 && allIndices[k][1] == 0) {
                    numberOfDifferentTaxa += 1;
                    allIndices[k][0] = indexInTaxa[j];
                    allIndices[k][1] = 1;
                    indexInMinDiances[j] = k;
                    break;
                }
            }
        }
        
        // Set weight factor
        int weightFactor = 1;
        if (weight == 1) {
            weightFactor = size;
        } else if (weight == 2) {
            weightFactor = numberOfDifferentTaxa;
        }
        
        double normFactor = 1.0;
        if (norm == 1) {
            normFactor = size;
        } else if (norm == 2) {
            normFactor = log(size);
        } else if (norm == 3) {
            normFactor = numberOfDifferentTaxa;
        }
        
        
        
        // miniPairs
        if (miniPairs != 0) {
                    
            int ***minDistances = (int ***) calloc(sizeof(int **), numberOfDifferentTaxa);
            for (int j = 0; j < numberOfDifferentTaxa; j++) {
                minDistances[j] = (int **) calloc(sizeof(int *), numberOfDifferentTaxa);
                for (int k = 0; k < numberOfDifferentTaxa; k++) {
                    minDistances[j][k] = (int *) calloc(sizeof(int *), 3);
                }
            }
            int adddd = 0;
            while (adddd == 0) {
                adddd = 1;
                for (int j = 0; j < size; j++) {
                    for (int k = 0; k < size - j; k++) {
                        if (dist[j][k] > 0) {
                            int ii = indexInMinDiances[j];
                            int jj = indexInMinDiances[k + j];
                            if (ii == jj) {
                                continue;
                            }
                            if (jj < ii) {
                                int tmp = jj;
                                jj = ii;
                                ii = tmp;
                            }
                            if (minDistances[ii][jj][2] == 0) {
                                minDistances[ii][jj][0] = j;
                                minDistances[ii][jj][1] = k + j;
                                minDistances[ii][jj][2] = dist[j][k];
                            } else if (dist[j][k] < minDistances[ii][jj][2]) {
                                minDistances[ii][jj][0] = j;
                                minDistances[ii][jj][1] = k + j;
                                minDistances[ii][jj][2] = dist[j][k];
                            }
                        }
                    }
                }
                for (int j = 0; j < numberOfDifferentTaxa; j++) {
                    for (int k = 0; k < numberOfDifferentTaxa; k++) {
                        int ii = allIndices[j][0];
                        int jj = allIndices[k][0];
                        if (jj < ii) {
                            int tmp = jj;
                            jj = ii;
                            ii = tmp;
                        }
                        if (minDistances[j][k][2] != 0 && j != k) {
                            adddd = 0;
			    int numberOfPairs = allIndices[j][1];
			    if (allIndices[k][1] < numberOfPairs) {
				numberOfPairs = allIndices[k][1];
			    }
                            if (taxaDistances[ii][jj][2 * i] == 0 || minDistances[j][k][2] / normFactor < miniPairs * taxaDistances[ii][jj][2 * i + 1] / taxaDistances[ii][jj][2 * i]) {
				taxaDistances[ii][jj][2 * i] += 1;
                                taxaDistances[ii][jj][2 * i + 1] += minDistances[j][k][2] / normFactor;
                            }
                            for (int l = 0; l < size; l++) {
                                for (int m = 0; m < size - l; m++) {
                                    if ((l == minDistances[j][k][0] && indexInTaxa[m + l] == indexInTaxa[minDistances[j][k][1]]) || (m + l == minDistances[j][k][1] && indexInTaxa[l] == indexInTaxa[minDistances[j][k][0]]) || (l == minDistances[j][k][1] && indexInTaxa[m + l] == indexInTaxa[minDistances[j][k][0]]) || (m + l == minDistances[j][k][0] && indexInTaxa[l] == indexInTaxa[minDistances[j][k][1]])) {
                                        dist[l][m] = 0;
                                    }
                                }
                            }
                            minDistances[j][k][0] = 0;
                            minDistances[j][k][1] = 0;
                            minDistances[j][k][2] = 0;
                        }
                    }
                }
            }
        } else {
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
                        if (square == 1) {
                            currentDistance = pow(dist[j][k], 2.0);
                        }
                        if (closeFriends != 0) {
                            if (currentDistance > closeFriends) {
                                currentDistance = closeFriends;
                            }
                        }
                        if (mini) {
                            if (taxaDistances[ii][jj][2 * i + 1] > currentDistance * weightFactor / normFactor) {
                                taxaDistances[ii][jj][2 * i + 1] = currentDistance * weightFactor / normFactor;
                            } else if (taxaDistances[ii][jj][2 * i + 1] == 0) {
                                taxaDistances[ii][jj][2 * i + 1] = currentDistance * weightFactor / normFactor;
                                taxaDistances[ii][jj][2 * i] = 1 * weightFactor;
                            }
                        } else if (closeFriends != 0) {
                            if (taxaDistances[ii][jj][2 * i + 1] > currentDistance * weightFactor / normFactor) {
                                taxaDistances[ii][jj][2 * i + 1] = currentDistance * weightFactor / normFactor;
                            } else if (taxaDistances[ii][jj][2 * i + 1] == 0) {
                                taxaDistances[ii][jj][2 * i + 1] = currentDistance * weightFactor / normFactor;
                                taxaDistances[ii][jj][2 * i] = 1 * weightFactor;
                            }
                        } else {
                            taxaDistances[ii][jj][2 * i + 1] += currentDistance * weightFactor / normFactor;
                            taxaDistances[ii][jj][2 * i] += 1 * weightFactor;
                        }
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
    free(trees);
    
    double **distance = (double **) calloc(sizeof(double*), numberOfTaxa);
    for (int i = 0; i < numberOfTaxa; i++) {
        distance[i] = (double *) calloc(sizeof(double), numberOfTaxa);
    }

    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            if (i < j) {
                int count = 0;
                double *sorted = (double *) calloc(sizeof(double), numberOfTrees);
                for (int k = 0; k < numberOfTrees; k++) {
                    if (taxaDistances[i][j][2 * k] != 0) {
                        double currentDistance = taxaDistances[i][j][2 * k + 1];
                        if (quartil != 0) {
                            sorted[count] = currentDistance;
                        } else {
                            if (square == 1) {
                                currentDistance = sqrt(currentDistance);
                            }
                            if (ustar) {
                                distance[i][j] += currentDistance / taxaDistances[i][j][2 * k];
                                count++;
                            } else {
                                distance[i][j] += currentDistance;
                                count += taxaDistances[i][j][2 * k];
                            }
                        }
                    }
                }
                if (quartil != 0) {
                    qsort(sorted, count, sizeof(double), compare);
                    distance[i][j] = sorted[(int) quartil * count];
                } else {
                    if (count != 0) {
                        distance[i][j] = distance[i][j] / count;
                    }
                }
                free(sorted);
            } else {
                distance[i][j] = distance[j][i];
            }
        }
    }
//    makeTree2(distance, numberOfTaxa, speciesTree, taxa);
    makeTreeFromDistanceArray(distance, numberOfTaxa, speciesTree, taxa);
    
    for (int i = 0; i < numberOfTaxa; i++) {
        for (int j = 0; j < numberOfTaxa; j++) {
            free(taxaDistances[i][j]);
        }
        free(taxaDistances[i]);
        free(taxa[i]);
    }
    free(taxaDistances);
    free(taxa);
}

void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames) {
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
