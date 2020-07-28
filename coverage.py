import dendropy
import numpy as np

names = ["ensembl_98_ncrna_vertebrates", "cyano_empirical"]
speciesName = "species_mapped.newick"
geneNames = ["geneTrees_mapped.newick"]

for name in names:
  print("\n" + name + "\n")
  for geneName in geneNames:
    file = open(name + "_" + geneName + "_cov", "w")
    print(geneName + "\n")
    string = ""
    string += (name + "\n")
    string += (geneName + "\n")
    
    speciesTree = open(name + "/" + speciesName, "r").read()
    geneTrees = open(name + "/" + geneName, "r").read()

    dataSet = dendropy.DataSet()
    dataSet.read(data=speciesTree, schema = "newick")
    namespace = [""] * len(dataSet.tree_lists[0].taxon_namespace)

    for i in range(len(namespace)):
      namespace[i] = str(dataSet.tree_lists[0].taxon_namespace[i])[1:-1]

    #Number of trees it appears in
    coverage = np.zeros(len(namespace))
    
    #Avg. number of appearances in trees it appears
    av_cover = np.zeros(len(namespace))
    
    #Appearances together with others
    crossCoverage = np.zeros(shape = (len(namespace), len(namespace)))
    
    rel_crossCoverage = np.zeros(shape = (len(namespace), len(namespace)))

    size = 0
    for geneTree in geneTrees.split(";"):
      size += 1
      for i in range(len(namespace)):
        if (geneTree.find(namespace[i] + ":") != -1):
          coverage[i] += 1
          av_cover[i] += geneTree.count(namespace[i] + ":")
          for j in range(len(namespace)):
            if (geneTree.find(namespace[j] + ":") != -1):
              crossCoverage[i][j] += 1

    size -= 1
    
    print(str(int(size)) + "\n")
    string += ("size: " + str(int(size)) + "\n")
    
    for i in range(len(namespace)):
        for j in range(len(namespace)):
            rel_crossCoverage[i][j] = crossCoverage[i][j] / coverage[i]
        av_cover[i] = av_cover[i] / coverage[i]
        name = namespace[i]
        cov = str(int(coverage[i]))
        if (coverage[i] < 1000):
            cov = "0" + cov
        av_cov = "{:.2f}".format(av_cover[i])
        if (av_cover[i] < 10):
            av_cov = "0" + av_cov
        print(name + "\t" + cov + "\t" + av_cov)

    string += ("\ncoverage:\n")
    for i in range(len(namespace)):
        if (coverage[i] < 1000):
            string += (namespace[i] + "\t0" + str(int(coverage[i])) + "\n")
        else:
            string += (namespace[i] + "\t" + str(int(coverage[i])) + "\n")
    
    string += ("\nav_coverage:\n")
    
    for i in range(len(namespace)):
        if (av_cover[i] < 10):
            string += (namespace[i] + "\t0" + "{:.2f}".format(av_cover[i]) + "\n")
        else:
            string += (namespace[i] + "\t{:.2f}".format(av_cover[i]) + "\n")

    string += ("\ncross_cov\n\t")
    
    for i in range(len(namespace)):
        string += (namespace[i] + "\t")
        
    for i in range(len(namespace)):
        string += ("\n" + namespace[i] + "\t")
        for j in range(len(namespace)):
            string += (str(int(crossCoverage[i][j])) + "\t")
            
    
    string += ("\nrel_cross_cov\n\t")
    
    for i in range(len(namespace)):
        string += (namespace[i] + "\t")
        
    for i in range(len(namespace)):
        string += ("\n" + namespace[i] + "\t")
        for j in range(len(namespace)):
            string += ("{:.2f}".format(rel_crossCoverage[i][j]) + "\t")
    
    file.write(string)
    file.close()
