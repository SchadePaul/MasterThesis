#include "node.h"
#include "parse.h"
#include "tree.h"
#include "neighborJoining.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>

int static const expander = 200;

void njstFromFile(struct node **root, const char *filename) {
    
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
            allLeafDistances[i][j] = (double *) calloc(sizeof(double *), numberOfTrees);
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
        leafToLeafDistance(tree[i], dist, size, name);
        
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
//                allLeafDistances[ii][jj][0] += 1;
//                if ((int) allLeafDistances[ii][jj][0] % expander == 0) {
//                    double *tmp = allLeafDistances[ii][jj];
//                    allLeafDistances[ii][jj] = calloc(sizeof(double), ((tmp[0] / expander) + 1) * expander);
//                    for (int l = 0; l < tmp[0]; l++) {
//                        allLeafDistances[ii][jj][l] = tmp[l];
//                    }
//                    free(tmp);
//                }
                if (allLeafDistances[ii][jj][i] > dist[j][k] || allLeafDistances[ii][jj][i] == 0) {
                    allLeafDistances[ii][jj][i] = dist[j][k];
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
                int counter = 0;
                for (int k = 0; k < numberOfTrees; k++) {
                    distance[i][j] += allLeafDistances[i][j][k];
                    counter += (allLeafDistances[i][j][k] != 0);
                }
                distance[i][j] = distance[i][j] / (double) counter;
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
