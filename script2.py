import dendropy
import subprocess
import math
import time
import matplotlib.pyplot as plt

fileSpeciesTree = open("test_100/01/outAst", "r")
fileOutputTree = open("test_100/01/s_tree.trees", "r")
strFileOutputTree = fileOutputTree.read()
strFileSpeciesTree = fileSpeciesTree.read()
fileOutputTree.close()
fileSpeciesTree.close()
            
dataSet = dendropy.DataSet()
dataSet.read(data=strFileSpeciesTree, schema = "newick")
dataSet.read(data=strFileOutputTree, schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)
distance = dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[1][0], is_bipartitions_updated=False)
print(distance)
        
