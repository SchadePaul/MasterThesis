#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>


void readFileToTrees(struct node ***trees, const char *filename, int *numberOfTrees, char ***allLeafNames, int *numberOfLeafNames, int tagged) {
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


void taggedNJFromFile(struct node **root, const char *filename, int norm, int rooted, int weighted, int square) {
    int numberOfTrees;
    char **allLeafNames;
    int numberOfLeafNames = 0;
    struct node **tree;
    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 1);
    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
        for (int j = 0; j < numberOfLeafNames; j++) {
            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), 2);
        }
    }
    
    for (int i = 0; i < numberOfTrees; i++) {
        if (rooted) {
            tagAndRoot(&(tree[i]));
        } else {
            scoreAndTag(tree[i]);
        }
        
        
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
        leafToLeafDistance(tree[i], dist, size, name, norm, 0);

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

        int multiply = 1;
        if (weighted == 1) {
            multiply = size;
//            printf("weighted:\t%d\n",size);
        }
        
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size - j; k++) {
                if (dist[j][k] > 0) {
                    int ii = index[j];
                    int jj = index[k + j];
                    if (jj < ii) {
                        int tmp = ii;
                        ii = jj;
                        jj = tmp;
                    }
                    if (square == 1) {
                        allLeafDistances[ii][jj][0] += 1 * multiply;
                        allLeafDistances[ii][jj][1] += dist[j][k] * dist[j][k] * multiply;
                    } else {
                        allLeafDistances[ii][jj][0] += 1 * multiply;
                        allLeafDistances[ii][jj][1] += dist[j][k] * multiply;
                    }
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
                if (allLeafDistances[i][j][0] != 0) {
                    if (square == 1) {
                        distance[i][j] = sqrt(allLeafDistances[i][j][1]);
                    } else {
                        distance[i][j] = allLeafDistances[i][j][1];
                    }
                    distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
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

void miniNJFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength) {
    
    int numberOfTrees;
    char **allLeafNames;
    int numberOfLeafNames = 0;
    struct node **tree;
    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 0);
    
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
        leafToLeafDistance(tree[i], dist, size, name, norm, branchLength);

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
                int multiply = 1;
                if (weighted == 1) {
                    multiply = size;
                }

                if (allLeafDistances[ii][jj][i + 1] > dist[j][k] * multiply) {
                    allLeafDistances[ii][jj][i + 1] = dist[j][k] * multiply;
                } else if (allLeafDistances[ii][jj][i + 1] == 0) {
                    allLeafDistances[ii][jj][i + 1] = dist[j][k] * multiply;
                    allLeafDistances[ii][jj][0] += 1 * multiply;
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
                int kMax = numberOfTrees;
                for (int k = 0; k < kMax; k++) {
                    distance[i][j] += allLeafDistances[i][j][k + 1];
                }
                if (allLeafDistances[i][j][0] != 0) {
                    distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
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

void njstFromFile(struct node **root, const char *filename, int norm, int weighted, int branchLength, int square) {
    
    int numberOfTrees;
    char **allLeafNames;
    int numberOfLeafNames = 0;
    struct node **tree;
    readFileToTrees(&tree, filename, &numberOfTrees, &allLeafNames, &numberOfLeafNames, 0);
    
    // All distance Array [i][j][0] sum of distances, [i][j][1] number of summands
    double ***allLeafDistances = (double ***) calloc(sizeof(double **), numberOfLeafNames);
    for (int i = 0; i < numberOfLeafNames; i++) {
        allLeafDistances[i] = (double **) calloc(sizeof(double *), numberOfLeafNames);
        for (int j = 0; j < numberOfLeafNames; j++) {
            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), 2);
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
        leafToLeafDistance(tree[i], dist, size, name, norm, branchLength);

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

        int multiply = 1;
        if (weighted == 1) {
            multiply = size;
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

                
                if (square == 1) {
                    allLeafDistances[ii][jj][0] += 1 * multiply;
                    allLeafDistances[ii][jj][1] += dist[j][k] * dist[j][k] * multiply;
                } else {
                    allLeafDistances[ii][jj][0] += 1 * multiply;
                    allLeafDistances[ii][jj][1] += dist[j][k] * multiply;
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
                if (allLeafDistances[i][j][0] != 0) {
                    if (square == 1) {
                        distance[i][j] = sqrt(allLeafDistances[i][j][1]);
                    } else {
                        distance[i][j] = allLeafDistances[i][j][1];
                    }
                    distance[i][j] = distance[i][j] / allLeafDistances[i][j][0];
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
