import dendropy
import subprocess
name = "test_new"
inputName = "allTrees"
outputName = "output"
number = 1000
number2 = 10
number3 = 100
cmd = subprocess.Popen(["../externalWork/SimPhy_1.0.2/bin/simphy_mac64", "-sl", "f:" + str(number3), "-rs", str(number2), "-rl", "f:" + str(number), "-rg", "1", "-sb", "f:0.000000005", "-sd", "f:0", "-st", "ln:21.25,0.2", "-so", "f:1", "-si", "f:1", "-sp", "f:470000000", "-su", "ln:-21.9,0.1", "-hh", "f:1", "-hs", "ln:1.5,1", "-hl", "ln:1.551533,0.6931472", "-hg", "ln:1.5,1", "-cs", "9644", "-v", "3", "-o", name, "-ot", "0", "-op", "1", "-lb", "f:0.00000000049", "-ld", "f:0.00000000049", "-lt", "f:0"])
cmd.communicate()
d = 0
for j in range(1, number2 + 1):
    strr = ""
    if (j < 10):
        strr = "0"
        
    file = open(inputName + strr + str(j), "w")
    for i in range(1, number + 1):
        ignore = False
        strrr = ""
        if (i < 10):
            strrr = "000"
        elif (i < 100):
            strrr = "00"
        elif (i < 1000):
            strrr = "0"
        
            
        file2 = open(name + "/" + strr + str(j) + "/g_trees" + strrr + str(i) + ".trees", "r")
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
    cmd2 = subprocess.Popen(["time", "./njst/main", "-i", inputName + strr + str(j), "-o", outputName + strr + str(j)])
    cmd2.communicate()

    fileSpeciesTree = open(name + "/" + strr + str(j) + "/s_tree.trees", "r")
    strFileSpeciesTree = fileSpeciesTree.read()
    fileSpeciesTree.close()
    fileOutputTree = open(outputName + strr + str(j), "r")
    strFileOutputTree = fileOutputTree.read()
    fileOutputTree.close()

    dataSet = dendropy.DataSet()
    dataSet.read(data=strFileSpeciesTree, schema = "newick")
    dataSet.read(data=strFileOutputTree, schema = "newick", taxon_namespace=dataSet.tree_lists[0].taxon_namespace)

    d += dendropy.calculate.treecompare.symmetric_difference(dataSet.tree_lists[0][0], dataSet.tree_lists[1][0], is_bipartitions_updated=False)
    print(d)
print(str(d/number2))
