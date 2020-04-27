import dendropy
import subprocess
import math
import numpy
import time
import matplotlib.pyplot as plt


inputName = "allTrees"
outputName = "output"
numberOfGeneTrees = 500
nameOfDirectory = "test_" + str(numberOfGeneTrees)
numberOfSets = 100
numberOfLeaves = [15,25,50,75,100,150,200]
legendNames = ["NJst", "NJst weighted", "NJst normed by tree size", "MiniNJ", "MiniNJ weighted", "MiniNJ normed by tree size"]
fileResults = open("results", "w")
exponentGeneTree = int(math.log10(numberOfGeneTrees))
exponentSet = int(math.log10(numberOfSets))

for leaves in numberOfLeaves:
    cmd = subprocess.Popen(["../externalWork/SimPhy_1.0.2/bin/simphy_mac64", "-sl", "f:" + str(leaves - 1), "-rs", str(numberOfSets), "-rl", "f:" + str(numberOfGeneTrees), "-rg", "1", "-sb", "f:0.000000005", "-sd", "f:0", "-st", "ln:21.25,0.2", "-so", "f:1", "-si", "f:1", "-sp", "f:470000000", "-su", "ln:-21.9,0.1", "-hh", "f:1", "-hs", "ln:1.5,1", "-hl", "ln:1.551533,0.6931472", "-hg", "ln:1.5,1", "-cs", "9644", "-v", "3", "-o", nameOfDirectory + "_" + str(leaves), "-ot", "0", "-op", "1", "-lb", "f:0.00000000049", "-ld", "f:0.00000000049", "-lt", "f:0"])
    cmd.communicate()

allDistances = [[[0 for k in range (len(numberOfLeaves))] for j in range(len(legendNames))] for i in range(2)]
for leaves in numberOfLeaves:
    distancesForSets = [[0 for j in range(numberOfSets)] for i in range(len(legendNames))]
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
#            currentTime = int(round(time.time() * 1000))
        inputStr = nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName
        outputStr = nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName
        
        cmd2 = subprocess.Popen(["./njst/main", "-i", inputStr, "-o", outputStr + "0"])
        cmd2.communicate()

        cmd2 = subprocess.Popen(["./njst/main", "-w", "-i", inputStr, "-o", outputStr + "1"])
        cmd2.communicate()
        
        cmd2 = subprocess.Popen(["./njst/main", "-n", "1", "-i", inputStr, "-o", outputStr + "2"])
        cmd2.communicate()

        cmd2 = subprocess.Popen(["./njst/main", "-m", "-i", inputStr, "-o", outputStr + "3"])
        cmd2.communicate()

        cmd2 = subprocess.Popen(["./njst/main", "-m", "-w", "-i", inputStr, "-o", outputStr + "4"])
        cmd2.communicate()
        
        cmd2 = subprocess.Popen(["./njst/main", "-m", "-n", "1", "-i", inputStr, "-o", outputStr + "5"])
        cmd2.communicate()


        fileSpeciesTree = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/s_tree.trees", "r")
        strFileSpeciesTree = fileSpeciesTree.read()
        fileSpeciesTree.close()

        dataSet = dendropy.DataSet()
        dataSet.read(data=strFileSpeciesTree, schema = "newick")

        for i in range(len(legendNames)):

            fileOutputTree = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + str(i), "r")
            strFileOutputTree = fileOutputTree.read()
            fileOutputTree.close()
            dataSet.read(data=strFileOutputTree, schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)
            distancesForSets[i][set - 1] = dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[i + 1][0], is_bipartitions_updated=False) / (2 * len(dataSet.tree_lists[0].taxon_namespace) - 6)
    for i in range(len(legendNames)):
        allDistances[0][i][numberOfLeaves.index(leaves)] = numpy.mean(distancesForSets[i])
        allDistances[1][i][numberOfLeaves.index(leaves)] = numpy.std(distancesForSets[i])


fileResults.write("Number of Sets: " + str(numberOfSets) + ", Gene Families: " + str(numberOfGeneTrees) + ", number of Species: " + str(numberOfLeaves) + "\n")

for i in range(2):
    if (i == 0):
        fileResults.write("average:\n")
    elif (i == 1):
        fileResults.write("standard deviation:\n")
    for j in range(len(legendNames)):
        fileResults.write(str(legendNames[j]) + "\n")
        for k in range(len(numberOfLeaves)):
            fileResults.write(str(allDistances[i][j][k]) + ",")
        fileResults.write("\n")
    fileResults.write("\n")
fileResults.write("\n")
fileResults.close()

for i in range(len(legendNames)):
    plt.errorbar(numberOfLeaves, allDistances[0][i], yerr=allDistances[1][i], fmt='x')
plt.legend(legendNames)
plt.xlabel("Number of species")
plt.ylabel("RF-Distance")
plt.title("RF-Distance for " + str(numberOfGeneTrees) + " gene families")
plt.savefig("RF-Distance for " + str(numberOfGeneTrees) + " gene families.pdf")