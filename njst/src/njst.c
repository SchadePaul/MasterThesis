#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>

void sortSmall(double* distances, int from, int to) {
    int swap = 1;
    while (swap == 1) {
        swap = 0;
        for (int i = from; i < to; i++) {
            if (distances[i] > distances[i + 1]) {
                double tmp = distances[i];
                distances[i] = distances[i + 1];
                distances[i + 1] = tmp;
                swap = 1;
            }
        }
    }
}

void sort(double* distances, int from, int to, int seed) {
    if ((to - from) < 15) {
        sortSmall(distances, from, to);
        return;
    }
    int pivotIndex = (seed % (to - from)) + from;
    double pivot = distances[pivotIndex];
    int left = from;
    int right = to;
    int same = 0;
    while (1) {
        if (distances[left] <= pivot && left < to) {
            if (distances[left] == pivot) {
                same++;
            }
            left++;
        }
        if (distances[right] > pivot && right > from) {
            right--;
        }
        if (left >= right) {
            break;
        }
        double tmp = distances[left];
        distances[left] = distances[right];
        distances[right] = tmp;
    }
    if (same == to - from) {
        return;
    }
    sort(distances, from, right, seed + 1);
    sort(distances, left, to, seed + 1);
    return;
}

void njstFromFile(struct node **root, const char *filename, int branchLength, int minNJst, int normDistance, double quartil) {
    
    // read file to array of chars
    char **newickTree;
    int numberOfTrees;
    readFileToArray(filename, &newickTree, &numberOfTrees);
    
    if (errno != 0) {
        return;
    }
    
    
    // array of chars to trees
    char **allLeafNames;
    int numberOfLeafNames = 0;
    struct node **tree = (struct node **) calloc(sizeof(struct node *), numberOfTrees);
    for (int i = 0; i < numberOfTrees; i++) {
        // parse newick format to tree, add leaf names to allLeafNamesArray
        newickTreeToTree(newickTree[i], &tree[i], &allLeafNames, &numberOfLeafNames);
        free(newickTree[i]);
    }
    free(newickTree);
    
    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
        for (int j = 0; j < numberOfLeafNames; j++) {
            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), numberOfTrees + 1);
        }
    }
    
    for (int i = 0; i < numberOfTrees; i++) {
        // Number of Leaves in this tree
        int size = tree[i]->numberOfLeaves;
        
        // Half-Matrix for distances between leaves
        double **dist = (double **) calloc(sizeof(double *), size);
        for (int j = 0; j < size; j++) {
            dist[j] = (double *) calloc(sizeof(double), size - j);
        }
        
        // Array of leave names
        char **name = (char **) calloc(sizeof(char *), size);
        for (int j = 0; j < size; j++) {
            name[j] = (char *) calloc(sizeof(char), maxNameLength);
        }
        
        // Compute leaf distances
        leafToLeafDistance(tree[i], dist, size, name, normDistance, branchLength);
        
        // Add to overall distance matrix
        
        // Find name in all Names
        int *index = (int *) calloc(sizeof(int), size);
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < numberOfLeafNames; k++) {
                if (strcmp(name[j], allLeafNames[k]) == 0) {
                    index[j] = k;
                    break;
                }
            }
        }
        
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size - j; k++) {
                int ii = index[j];
                int jj = index[k + j];
                if (jj < ii) {
                    int tmp = ii;
                    ii = jj;
                    jj = tmp;
                }
                if (minNJst) {
                    // minNJst
                    if (allLeafDistances[ii][jj][i + 1] > dist[j][k]) {
                        allLeafDistances[ii][jj][i + 1] = dist[j][k];
                    } else if (allLeafDistances[ii][jj][i + 1] == 0) {
                        allLeafDistances[ii][jj][i + 1] = dist[j][k];
                        allLeafDistances[ii][jj][0] += 1;
                    }
                } else {
                    // NJst
                    allLeafDistances[ii][jj][0] += 1;
                    if (((int) allLeafDistances[ii][jj][0] % numberOfTrees) == 0) {
                        double *tmp = allLeafDistances[ii][jj];
                        allLeafDistances[ii][jj] = (double *) calloc(sizeof(double), ((tmp[0] / numberOfTrees) + 1) * numberOfTrees);
                        for (int index = 0; index < tmp[0]; index++) {
                            allLeafDistances[ii][jj][index] = tmp[index];
                        }
                        free(tmp);
                    }
                    allLeafDistances[ii][jj][(int) allLeafDistances[ii][jj][0]] = dist[j][k];
                    
                }
            }
        }
        free(index);
        
        for (int j = 0; j < size; j++) {
            free(dist[j]);
            free(name[j]);
        }
        free(name);
        free(dist);
    }
    
    
    
    
    
    
    double **distance = (double **) calloc(sizeof(double*), numberOfLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        distance[i] = (double *) calloc(sizeof(double), numberOfLeafNames);
    }
    
    for (int i = 0; i < numberOfLeafNames; i++) {
        for (int j = 0; j < numberOfLeafNames; j++) {
            if (i < j) {
                if (quartil != 0) {
                    if (minNJst) {
                        sort(allLeafDistances[i][j], 1, numberOfTrees + 1, 1);
                    } else {
                        sort(allLeafDistances[i][j], 1, (int) allLeafDistances[i][j][0], 1);
                    }
                }
                
                int kMax = allLeafDistances[i][j][0];
                if (minNJst) {
                    kMax = numberOfTrees;
                }
                
                // take median
                if (quartil != 0) {
                    if (minNJst) {
                        distance[i][j] = allLeafDistances[i][j][(numberOfTrees - (int) allLeafDistances[i][j][0]) + (int) (allLeafDistances[i][j][0] * quartil)];
                    } else {
                        distance[i][j] = allLeafDistances[i][j][(int) (allLeafDistances[i][j][0] * quartil)];
                    }
                } else {
                    for (int k = 0; k < kMax; k++) {
                        distance[i][j] += allLeafDistances[i][j][k + 1];
                    }
                    if (allLeafDistances[i][j][0] != 0) {
                        distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
                    }
                }
            } else {
                distance[i][j] = distance[j][i];
            }
        }
    }
    
    makeTreeFromDistanceArray(distance, numberOfLeafNames, root, allLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        for (int j = 0; j < numberOfLeafNames; j++) {
            free(allLeafDistances[i][j]);
        }
        free(allLeafDistances[i]);
        free(allLeafNames[i]);
    }
    free(allLeafDistances);
    free(allLeafNames);
}
