#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

static void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames);
void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root, double miniPairs);

void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root, double miniPairs) {
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
        }
        
//        printf("%s\t%d\t%d\t\t", current->name, current->numberOfChildren, current->numberOfLeaves);
//        int childs = current->numberOfChildren;
//        printf("tag: %d\t\t tag2: ", current->tag);
//        for (int i = 0 ; i < (childs * childs - childs) / 2; i++) {
//            printf("%d\t", current->tag2[i]);
//        }
//        printf("\n");
//        while (1) {
//            while (current->firstChild != 0) {
//                current = current->firstChild;
//                printf("%s\t%d\t%d\t\t", current->name, current->numberOfChildren, current->numberOfLeaves);
//                int childs = current->numberOfChildren;
//                printf("tag: %d\t\t tag2: ", current->tag);
//                for (int i = 0 ; i < (childs * childs - childs) / 2; i++) {
//                    printf("%d\t", current->tag2[i]);
//                }
//                printf("\n");
//            }
//            while (current->nextSibling == 0) {
//                current = current->parent;
//                if (current == trees[i]) {
//                    break;
//                }
//            }
//            if (current == trees[i]) {
//                break;
//            }
//            current = current->nextSibling;
//            printf("%s\t%d\t%d\t\t", current->name, current->numberOfChildren, current->numberOfLeaves);
//            int childs = current->numberOfChildren;
//            printf("tag: %d\t\t tag2: ", current->tag);
//            for (int i = 0 ; i < (childs * childs - childs) / 2; i++) {
//                printf("%d\t", current->tag2[i]);
//            }
//            printf("\n");
//        }
//        printf("\n");
        if (tag) {
//            char **names = (char **) calloc(sizeof(char *), size);
//            for (int j = 0; j < size; j++) {
//                names[j] = (char *) calloc(sizeof(char), maxNameLength);
//            }
            scoreAndTag2(trees[i]);
//            for (int j = 0; j < size; j++) {
//                free(names[j]);
//            }
//            free(names);
        }
        
        
        // Set weight factor
        int weightFactor = 1;
        if (weight == 1) {
            weightFactor = size;
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
        leafToLeafDistance(trees[i], dist, taxaInTree, norm, branchLength);
        
        int *indexx = calloc(sizeof(int), 1);
        deletedTaggedDistance(trees[i], dist, indexx);
        
//        for (int i = 0; i < size; i++) {
//            for (int j = 0; j < size; j++) {
//                if (i < j) {
//                    printf("%d\t", (int) dist[i][j - i]);
//                } else {
//                    printf("0\t");
//                }
//            }
//            printf("\n");
//        }
//        printf("\n");
        
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
        
        // miniPairs
        if (miniPairs != 0) {
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
                            if (taxaDistances[ii][jj][2 * i] == 0 || minDistances[j][k][2] < miniPairs * taxaDistances[ii][jj][2 * i + 1] / taxaDistances[ii][jj][2 * i]) {
				taxaDistances[ii][jj][2 * i] += 1;
                                taxaDistances[ii][jj][2 * i + 1] += minDistances[j][k][2];
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
                        if (mini) {
                            if (taxaDistances[ii][jj][2 * i + 1] > currentDistance * weightFactor) {
                                taxaDistances[ii][jj][2 * i + 1] = currentDistance * weightFactor;
                            } else if (taxaDistances[ii][jj][2 * i + 1] == 0) {
                                taxaDistances[ii][jj][2 * i + 1] = currentDistance * weightFactor;
                                taxaDistances[ii][jj][2 * i] = 1 * weightFactor;
                            }
                        } else {
                            taxaDistances[ii][jj][2 * i + 1] += currentDistance * weightFactor;
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
                for (int k = 0; k < numberOfTrees; k++) {
                    if (taxaDistances[i][j][2 * k] != 0) {
                        double currentDistance = taxaDistances[i][j][2 * k + 1];
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
                if (count != 0) {
                    distance[i][j] = distance[i][j] / count;
                }
            } else {
                distance[i][j] = distance[j][i];
            }
        }
    }
//    printf("\t");
//    for (int i = 0; i < numberOfTaxa; i++) {
//        printf("%s\t", taxa[i]);
//    }
//    
//    printf("\n\n");
//    for (int i = 0; i < numberOfTaxa; i++) {
//        printf("%s\t", taxa[i]);
//        for (int j = 0; j < numberOfTaxa; j++) {
//            printf("%d\t", (int) distance[i][j]);
//        }
//        printf("\n");
//    }
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
