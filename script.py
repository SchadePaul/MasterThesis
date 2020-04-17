import dendropy
import subprocess
import math
import time
import matplotlib.pyplot as plt

nameOfDirectory = "test"
inputName = "allTrees"
outputName = "output"
numberOfGeneTrees = 1000
numberOfSets = 10
numberOfLeaves = [15,25,50,75,100,150,200]
distances = []
times = []
fileResults = open("results", "w")
for leaves in numberOfLeaves:
    print("leaves: " + str(leaves))
#    cmd = subprocess.Popen(["../externalWork/SimPhy_1.0.2/bin/simphy_mac64", "-sl", "f:" + str(leaves), "-rs", str(numberOfSets), "-rl", "f:" + str(numberOfGeneTrees), "-rg", "1", "-sb", "f:0.000000005", "-sd", "f:0", "-st", "ln:21.25,0.2", "-so", "f:1", "-si", "f:1", "-sp", "f:470000000", "-su", "ln:-21.9,0.1", "-hh", "f:1", "-hs", "ln:1.5,1", "-hl", "ln:1.551533,0.6931472", "-hg", "ln:1.5,1", "-cs", "9644", "-v", "3", "-o", nameOfDirectory + "_" + str(leaves), "-ot", "0", "-op", "1", "-lb", "f:0.00000000049", "-ld", "f:0.00000000049", "-lt", "f:0"])
#    cmd.communicate()
    exponentGeneTree = int(math.log10(numberOfGeneTrees))
    exponentSet = int(math.log10(numberOfSets))
    for set in range(1, numberOfSets + 1):
        stringFillerSet = ""
        for i in range(1, exponentSet + 1):
            if (set < 10 ** i):
                stringFillerSet = stringFillerSet + "0"
        file = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "w")
        for genteTree in range(1, numberOfGeneTrees + 1):
            stringFillerGeneTree = ""
            for i in range(1, exponentGeneTree + 1):
                if (genteTree < 10 ** i):
                    stringFillerGeneTree = stringFillerGeneTree + "0"
            file2 = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/g_trees" + stringFillerGeneTree + str(genteTree) + ".trees" , "r")
            ignore = False
            while True:
                c = file2.read(1)
                if not c:
                    file2.close()
                    break
                if (c == '_'):
                    ignore = True
                elif (c == ':'):
                    ignore = False
                if not ignore:
                    file.write(c)
        file.close()
        currentTime = int(round(time.time() * 1000))
        cmd2 = subprocess.Popen(["./njst/main", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName])
        cmd2.communicate()
        executionTime = int(round(time.time() * 1000)) - currentTime
        print(executionTime)
        fileSpeciesTree = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/s_tree.trees", "r")
        fileOutputTree = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName, "r")
        strFileOutputTree = fileOutputTree.read()
        strFileSpeciesTree = fileSpeciesTree.read()
        fileOutputTree.close()
        fileSpeciesTree.close()
            
        dataSet = dendropy.DataSet()
        dataSet.read(data=strFileSpeciesTree, schema = "newick")
        dataSet.read(data=strFileOutputTree, schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)
        distance = dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[1][0], is_bipartitions_updated=False)
        distances.append(distance)
        times.append(executionTime)
        fileResults.write(str(leaves) + "\t" + str(numberOfGeneTrees) + "\t" + str(executionTime) + "\t" + str(distance) + "\n")
        
fileResults.close()

averageDistance = []
averageTime = []
for i in range(len(numberOfLeaves)):
    averageDistance.append(0)
    averageTime.append(0)
    for j in range(numberOfSets):
        averageDistance[i] += distances[i * numberOfSets + j]
        averageTime[i] += times[i * numberOfSets + j]
    averageDistance[i] = averageDistance[i] / numberOfSets
    averageTime[i] = averageTime[i] / numberOfSets
    
plt.figure(0)
plt.plot(numberOfLeaves, averageDistance)
plt.savefig("distance.pdf")
plt.figure(1)
plt.plot(numberOfLeaves, averageTime)
plt.savefig("time.pdf")

#for j in range(1, number2 + 1):
#    strr = ""
#    if (j < 10):
#        strr = "0"
#
#    file = open(inputName + strr + str(j), "w")
#    for i in range(1, number + 1):
#        ignore = False
#        strrr = ""
#        if (i < 10):
#            strrr = "000"
#        elif (i < 100):
#            strrr = "00"
#        elif (i < 1000):
#            strrr = "0"
#
#
#        file2 = open(name + "/" + strr + str(j) + "/g_trees" + strrr + str(i) + ".trees", "r")
#        while True:
#            c = file2.read(1)
#            if not c:
#                file2.close()
#                break
#            if (c == '_'):
#                ignore = True
#            elif (c == ':'):
#                ignore = False
#            if not ignore:
#                file.write(c)
#    file.close()
#    cmd2 = subprocess.Popen(["time", "./njst/main", "-i", inputName + strr + str(j), "-o", outputName + strr + str(j)])
#    cmd2.communicate()
#
#    fileSpeciesTree = open(name + "/" + strr + str(j) + "/s_tree.trees", "r")
#    strFileSpeciesTree = fileSpeciesTree.read()
#    fileSpeciesTree.close()
#    fileOutputTree = open(outputName + strr + str(j), "r")
#    strFileOutputTree = fileOutputTree.read()
#    fileOutputTree.close()
#
#    dataSet = dendropy.DataSet()
#    dataSet.read(data=strFileSpeciesTree, schema = "newick")
#    dataSet.read(data=strFileOutputTree, schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)
#
#    d += dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[1][0], is_bipartitions_updated=False)
#    print(d)
#print(str(d/number2))
