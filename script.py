import dendropy
import subprocess
import math
import time
import matplotlib.pyplot as plt

nameOfDirectory = "test_500"
inputName = "allTrees"
outputName = "output"
numberOfGeneTrees = 500
numberOfSets = 10
numberOfLeaves = [15,25,50,75,100,150,200]
fileResults = open("results", "w")
exponentGeneTree = int(math.log10(numberOfGeneTrees))
exponentSet = int(math.log10(numberOfSets))

#for leaves in numberOfLeaves:
#    cmd = subprocess.Popen(["../externalWork/SimPhy_1.0.2/bin/simphy_mac64", "-sl", "f:" + str(leaves), "-rs", str(numberOfSets), "-rl", "f:" + str(numberOfGeneTrees), "-rg", "1", "-sb", "f:0.000000005", "-sd", "f:0", "-st", "ln:21.25,0.2", "-so", "f:1", "-si", "f:1", "-sp", "f:470000000", "-su", "ln:-21.9,0.1", "-hh", "f:1", "-hs", "ln:1.5,1", "-hl", "ln:1.551533,0.6931472", "-hg", "ln:1.5,1", "-cs", "9644", "-v", "3", "-o", nameOfDirectory + "_" + str(leaves), "-ot", "0", "-op", "1", "-lb", "f:0.00000000049", "-ld", "f:0.00000000049", "-lt", "f:0"])
#    cmd.communicate()

for leaves in numberOfLeaves:
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
        cmd2 = subprocess.Popen(["./njst/main", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "0"])
        cmd2.communicate()
        
        cmd2 = subprocess.Popen(["./njst/main", "-n", "1", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "1"])
        cmd2.communicate()
        
        cmd2 = subprocess.Popen(["./njst/main", "-n", "2", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "2"])
        cmd2.communicate()
        
        cmd2 = subprocess.Popen(["./njst/main", "-n", "3", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "3"])
        cmd2.communicate()
        
        cmd2 = subprocess.Popen(["./njst/main", "-n", "4", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "4"])
        cmd2.communicate()
#        cmd2 = subprocess.Popen(["./njst/main", "-m", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "2"])
#        cmd2.communicate()
#
#        cmd2 = subprocess.Popen(["./njst/main", "-b", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "4"])
#        cmd2.communicate()
#
#        cmd2 = subprocess.Popen(["./njst/main", "-n", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "1"])
#        cmd2.communicate()
#
#        cmd2 = subprocess.Popen(["./njst/main", "-b", "-m", "-n", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "7"])
#        cmd2.communicate()
#
#        cmd2 = subprocess.Popen(["./njst/main", "-b", "-m", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "6"])
#        cmd2.communicate()
#
#        cmd2 = subprocess.Popen(["./njst/main", "-m", "-n", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "3"])
#        cmd2.communicate()
#
#        cmd2 = subprocess.Popen(["./njst/main", "-b", "-n", "-i", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + inputName, "-o", nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + "5"])
#        cmd2.communicate()
                
#            executionTime = int(round(time.time() * 1000)) - currentTime


        fileSpeciesTree = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/s_tree.trees", "r")
        strFileSpeciesTree = fileSpeciesTree.read()
        fileSpeciesTree.close()
            
        dataSet = dendropy.DataSet()
        dataSet.read(data=strFileSpeciesTree, schema = "newick")
            
        distance = []
            
        for i in range (5):
            
            fileOutputTree = open(nameOfDirectory + "_" + str(leaves) + "/" + stringFillerSet + str(set) + "/" + outputName + str(i), "r")
            strFileOutputTree = fileOutputTree.read()
            fileOutputTree.close()
            dataSet.read(data=strFileOutputTree, schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)
            distance.append(dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[i + 1][0], is_bipartitions_updated=False))
        print(str(leaves) + "\t" + str(numberOfGeneTrees) + "\t" + str(distance[0]) + "\t" + str(distance[1]) + "\t" + str(distance[2]) + "\t" + str(distance[3]) + "\t" + str(distance[4]) + "\n")
        fileResults.write(str(leaves) + "\t" + str(numberOfGeneTrees) + "\t" + str(distance[0]) + "\t" + str(distance[1]) + "\t" + str(distance[2]) + "\t" + str(distance[3]) + "\t" + str(distance[4]) + "\n")
        
fileResults.close()
