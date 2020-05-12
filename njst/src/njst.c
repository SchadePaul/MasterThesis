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
void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root);

void inferSpeciesTreeFromGeneTrees(struct node **speciesTree, const char *filename, int mini, int ustar, int norm, int branchLength, int weight, int square, int tag, int root) {
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
        if (tag) {
            char **names = (char **) calloc(sizeof(char *), size);
            for (int j = 0; j < size; j++) {
                names[j] = (char *) calloc(sizeof(char), maxNameLength);
            }
            scoreAndTag(trees[i], names);
            for (int j = 0; j < size; j++) {
                free(names[j]);
            }
            free(names);
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


//void taggedNJFromFile(struct node **root, const char *filename, int norm, int rooted, int weighted, int square) {
//    int numberOfTrees;
//    char **allLeafNames;
//    int numberOfLeafNames = 0;
//    struct node **tree;
//    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 1);
//    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
//    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), 2);
//        }
//    }
//
//    for (int i = 0; i < numberOfTrees; i++) {
//        if (rooted) {
//            tagAndRoot(&(tree[i]));
//        } else {
//            scoreAndTag(tree[i]);
//        }
//
//
//        // Number of Leaves in this tree
//        int size = tree[i]->numberOfLeaves;
//
//        // Half-Matrix for distances between leaves
//        double **dist = (double **) calloc(sizeof(double *), size);
//        for (int j = 0; j < size; j++) {
//            dist[j] = (double *) calloc(sizeof(double), size - j);
//        }
//
//        // Array of leave names
//        char **name = (char **) calloc(sizeof(char *), size);
//        for (int j = 0; j < size; j++) {
//            name[j] = (char *) calloc(sizeof(char), maxNameLength);
//        }
//
//        // Compute leaf distances
//        leafToLeafDistance(tree[i], dist, size, name, norm, 0);
//
//        // Add to overall distance matrix
//
//        // Find name in all Names
//        int *index = (int *) calloc(sizeof(int), size);
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < numberOfLeafNames; k++) {
//                if (strcmp(name[j], allLeafNames[k]) == 0) {
//                    index[j] = k;
//                    break;
//                }
//            }
//        }
//
//        int multiply = 1;
//        if (weighted == 1) {
//            multiply = size;
////            printf("weighted:\t%d\n",size);
//        }
//
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < size - j; k++) {
//                if (dist[j][k] > 0) {
//                    int ii = index[j];
//                    int jj = index[k + j];
//                    if (jj < ii) {
//                        int tmp = ii;
//                        ii = jj;
//                        jj = tmp;
//                    }
//                    if (square == 1) {
//                        allLeafDistances[ii][jj][0] += 1 * multiply;
//                        allLeafDistances[ii][jj][1] += dist[j][k] * dist[j][k] * multiply;
//                    } else {
//                        allLeafDistances[ii][jj][0] += 1 * multiply;
//                        allLeafDistances[ii][jj][1] += dist[j][k] * multiply;
//                    }
//                }
//            }
//        }
//        free(index);
//
//        for (int j = 0; j < size; j++) {
//            free(dist[j]);
//            free(name[j]);
//        }
//        free(name);
//        free(dist);
//
//    }
//
//    double **distance = (double **) calloc(sizeof(double*), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        distance[i] = (double *) calloc(sizeof(double), numberOfLeafNames);
//    }
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            if (i < j) {
//                if (allLeafDistances[i][j][0] != 0) {
//                    if (square == 1) {
//                        distance[i][j] = sqrt(allLeafDistances[i][j][1]);
//                    } else {
//                        distance[i][j] = allLeafDistances[i][j][1];
//                    }
//                    distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
//                }
//            } else {
//                distance[i][j] = distance[j][i];
//            }
//        }
//    }
//
//    makeTreeFromDistanceArray(distance, numberOfLeafNames, root, allLeafNames);
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            free(allLeafDistances[i][j]);
//        }
//        free(allLeafDistances[i]);
//        free(allLeafNames[i]);
//    }
//    free(allLeafDistances);
//    free(allLeafNames);
//
//}
//
//void miniNJFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength) {
//
//    int numberOfTrees;
//    char **allLeafNames;
//    int numberOfLeafNames = 0;
//    struct node **tree;
//    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 0);
//
//    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
//    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), numberOfTrees + 1);
//        }
//    }
//
//    for (int i = 0; i < numberOfTrees; i++) {
//        // Number of Leaves in this tree
//        int size = tree[i]->numberOfLeaves;
//
//        // Half-Matrix for distances between leaves
//        double **dist = (double **) calloc(sizeof(double *), size);
//        for (int j = 0; j < size; j++) {
//            dist[j] = (double *) calloc(sizeof(double), size - j);
//        }
//
//        // Array of leave names
//        char **name = (char **) calloc(sizeof(char *), size);
//        for (int j = 0; j < size; j++) {
//            name[j] = (char *) calloc(sizeof(char), maxNameLength);
//        }
//
//        // Compute leaf distances
//        leafToLeafDistance(tree[i], dist, size, name, norm, branchLength);
//
//        // Add to overall distance matrix
//
//        // Find name in all Names
//        int *index = (int *) calloc(sizeof(int), size);
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < numberOfLeafNames; k++) {
//                if (strcmp(name[j], allLeafNames[k]) == 0) {
//                    index[j] = k;
//                    break;
//                }
//            }
//        }
//
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < size - j; k++) {
//                int ii = index[j];
//                int jj = index[k + j];
//                if (jj < ii) {
//                    int tmp = ii;
//                    ii = jj;
//                    jj = tmp;
//                }
//                int multiply = 1;
//                if (weighted == 1) {
//                    multiply = size;
//                }
//
//                if (allLeafDistances[ii][jj][i + 1] > dist[j][k] * multiply) {
//                    allLeafDistances[ii][jj][i + 1] = dist[j][k] * multiply;
//                } else if (allLeafDistances[ii][jj][i + 1] == 0) {
//                    allLeafDistances[ii][jj][i + 1] = dist[j][k] * multiply;
//                    allLeafDistances[ii][jj][0] += 1 * multiply;
//                }
//            }
//        }
//        free(index);
//
//        for (int j = 0; j < size; j++) {
//            free(dist[j]);
//            free(name[j]);
//        }
//        free(name);
//        free(dist);
//
//    }
//
//    double **distance = (double **) calloc(sizeof(double*), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        distance[i] = (double *) calloc(sizeof(double), numberOfLeafNames);
//    }
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            if (i < j) {
//                int kMax = numberOfTrees;
//                for (int k = 0; k < kMax; k++) {
//                    distance[i][j] += allLeafDistances[i][j][k + 1];
//                }
//                if (allLeafDistances[i][j][0] != 0) {
//                    distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
//                }
//            } else {
//                distance[i][j] = distance[j][i];
//            }
//        }
//    }
//
//    makeTreeFromDistanceArray(distance, numberOfLeafNames, root, allLeafNames);
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            free(allLeafDistances[i][j]);
//        }
//        free(allLeafDistances[i]);
//        free(allLeafNames[i]);
//    }
//    free(allLeafDistances);
//    free(allLeafNames);
//
//}
//
//void ustar(struct node **root, const char *filename, int tag, int rooted) {
//    int numberOfTrees;
//    char **allLeafNames;
//    int numberOfLeafNames = 0;
//    struct node **tree;
//    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 0);
//
//    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
//    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), 2 * numberOfTrees);
//        }
//    }
//
//    for (int i = 0; i < numberOfTrees; i++) {
//
//        if (rooted) {
//            tagAndRoot(&(tree[i]));
//        } else if (tag) {
//            scoreAndTag(tree[i]);
//        }
//
//        // Number of Leaves in this tree
//        int size = tree[i]->numberOfLeaves;
//
//        // Half-Matrix for distances between leaves
//        double **dist = (double **) calloc(sizeof(double *), size);
//        for (int j = 0; j < size; j++) {
//            dist[j] = (double *) calloc(sizeof(double), size - j);
//        }
//
//        // Array of leave names
//        char **name = (char **) calloc(sizeof(char *), size);
//        for (int j = 0; j < size; j++) {
//            name[j] = (char *) calloc(sizeof(char), maxNameLength);
//        }
//
//        // Compute leaf distances
//        leafToLeafDistance(tree[i], dist, size, name, 0, 0);
//
//        // Add to overall distance matrix
//
//        // Find name in all Names
//        int *index = (int *) calloc(sizeof(int), size);
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < numberOfLeafNames; k++) {
//                if (strcmp(name[j], allLeafNames[k]) == 0) {
//                    index[j] = k;
//                    break;
//                }
//            }
//        }
//
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < size - j; k++) {
//                int ii = index[j];
//                int jj = index[k + j];
//                if (jj < ii) {
//                    int tmp = ii;
//                    ii = jj;
//                    jj = tmp;
//                }
//                if (dist[j][k] > 0) {
//                    allLeafDistances[ii][jj][2 * i + 1] += dist[j][k];
//                    allLeafDistances[ii][jj][2 * i] += 1;
//                }
//
//            }
//        }
//        free(index);
//
//        for (int j = 0; j < size; j++) {
//            free(dist[j]);
//            free(name[j]);
//        }
//        free(name);
//        free(dist);
//
//    }
//
//    double **distance = (double **) calloc(sizeof(double*), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        distance[i] = (double *) calloc(sizeof(double), numberOfLeafNames);
//    }
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            if (i < j) {
//                int kMax = numberOfTrees;
//                int count = 0;
//                for (int k = 0; k < kMax; k++) {
//                    if (allLeafDistances[i][j][2 * k] != 0) {
//                        distance[i][j] += allLeafDistances[i][j][2 * k + 1] / allLeafDistances[i][j][2 * k];
//                        count++;
//                    }
//                }
//                if (count != 0) {
//                    distance[i][j] = distance[i][j] / count;
//                }
//            } else {
//                distance[i][j] = distance[j][i];
//            }
//        }
//    }
//
//    makeTreeFromDistanceArray(distance, numberOfLeafNames, root, allLeafNames);
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            free(allLeafDistances[i][j]);
//        }
//        free(allLeafDistances[i]);
//        free(allLeafNames[i]);
//    }
//    free(allLeafDistances);
//    free(allLeafNames);
//
//}
//
//void njstFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength, int square) {
//
//    int numberOfTrees;
//    char **allLeafNames;
//    int numberOfLeafNames = 0;
//    struct node **tree;
//    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 0);
//
//    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
//    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), 2);
//        }
//    }
//
//    for (int i = 0; i < numberOfTrees; i++) {
//        // Number of Leaves in this tree
//        int size = tree[i]->numberOfLeaves;
//
//        // Half-Matrix for distances between leaves
//        double **dist = (double **) calloc(sizeof(double *), size);
//        for (int j = 0; j < size; j++) {
//            dist[j] = (double *) calloc(sizeof(double), size - j);
//        }
//
//        // Array of leave names
//        char **name = (char **) calloc(sizeof(char *), size);
//        for (int j = 0; j < size; j++) {
//            name[j] = (char *) calloc(sizeof(char), maxNameLength);
//        }
//
//        // Compute leaf distances
//        leafToLeafDistance(tree[i], dist, size, name, norm, branchLength);
//
//        // Add to overall distance matrix
//
//        // Find name in all Names
//        int *index = (int *) calloc(sizeof(int), size);
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < numberOfLeafNames; k++) {
//                if (strcmp(name[j], allLeafNames[k]) == 0) {
//                    index[j] = k;
//                    break;
//                }
//            }
//        }
//
//        int multiply = 1;
//        if (weighted == 1) {
//            multiply = size;
//        }
//
//        for (int j = 0; j < size; j++) {
//            for (int k = 0; k < size - j; k++) {
//
//                int ii = index[j];
//                int jj = index[k + j];
//                if (jj < ii) {
//                    int tmp = ii;
//                    ii = jj;
//                    jj = tmp;
//                }
//
//
//                if (square == 1) {
//                    allLeafDistances[ii][jj][0] += 1 * multiply;
//                    allLeafDistances[ii][jj][1] += dist[j][k] * dist[j][k] * multiply;
//                } else {
//                    allLeafDistances[ii][jj][0] += 1 * multiply;
//                    allLeafDistances[ii][jj][1] += dist[j][k] * multiply;
//                }
//
//            }
//        }
//        free(index);
//
//        for (int j = 0; j < size; j++) {
//            free(dist[j]);
//            free(name[j]);
//        }
//        free(name);
//        free(dist);
//
//    }
//
//    double **distance = (double **) calloc(sizeof(double*), numberOfLeafNames);
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        distance[i] = (double *) calloc(sizeof(double), numberOfLeafNames);
//    }
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            if (i < j) {
//                if (allLeafDistances[i][j][0] != 0) {
//                    if (square == 1) {
//                        distance[i][j] = sqrt(allLeafDistances[i][j][1]);
//                    } else {
//                        distance[i][j] = allLeafDistances[i][j][1];
//                    }
//                    distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
//                }
//            } else {
//                distance[i][j] = distance[j][i];
//            }
//        }
//    }
//
//    makeTreeFromDistanceArray(distance, numberOfLeafNames, root, allLeafNames);
//
//    for (int i = 0; i < numberOfLeafNames; i++) {
//        for (int j = 0; j < numberOfLeafNames; j++) {
//            free(allLeafDistances[i][j]);
//        }
//        free(allLeafDistances[i]);
//        free(allLeafNames[i]);
//    }
//    free(allLeafDistances);
//    free(allLeafNames);
//}
