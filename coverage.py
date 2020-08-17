import dendropy
from dendropy import treecalc
import numpy as np
import matplotlib.pyplot as plt


names = ["ensembl_98_ncrna_vertebrates", "cyano_empirical"]
shortNames = ["vertebrates", "cyano"]
speciesName = "species_mapped.newick"
geneNames = ["geneTrees_mapped.newick"]
shorGeneNames = ["gene"]

for ii, name in enumerate(names):
  
  for jj, geneName in enumerate(geneNames):

    # Read species and gene trees
    speciesTree = open(name + "/" + speciesName, "r").read()
    geneTrees = open(name + "/" + geneName, "r").read()

    dataSet = dendropy.DataSet()
    dataSet.read(data=speciesTree, schema = "newick")
    namespace = [""] * len(dataSet.tree_lists[0].taxon_namespace)

    for i in range(len(namespace)):
      namespace[i] = str(dataSet.tree_lists[0].taxon_namespace[i])[1:-1]

    # Number of trees it appears in
    coverage = np.zeros(len(namespace))
    
    # Avg. number of appearances in trees it appears
    av_cover = np.zeros(len(namespace))
    
    # TreeSize
    tree_sizes = []
    
#    #Appearances together with others
#    crossCoverage = np.zeros(shape = (len(namespace), len(namespace)))
#    rel_crossCoverage = np.zeros(shape = (len(namespace), len(namespace)))

    size = 0
    for j,geneTree in enumerate(geneTrees.split(";")):
      size += 1
      tree_sizes.append(0)
      for i in range(len(namespace)):
        if (geneTree.find(namespace[i] + ":") != -1):
          coverage[i] += 1
          count = geneTree.count(namespace[i] + ":")
          av_cover[i] += count
          tree_sizes[j] += count
#          for j in range(len(namespace)):
#            if (geneTree.find(namespace[j] + ":") != -1):
#              crossCoverage[i][j] += 1

    size -= 1
    
    # First Line Size
    string = str(size) + "\n\n"
    
    # Number species per tree
    for i in range(size):
        string += str(tree_sizes[i]) + "\n"
        
    string += "\n"
    
    for i in range(len(namespace)):
        av_cover[i] = av_cover[i] / coverage[i]
        string += (namespace[i] + "\t" + str(coverage[i]) + "\t" + str(av_cover[i]) + "\n")
    

    # File to write
    file = open(name + "_" + geneName + "_cov", "w")
    file.write(string)
    file.close()


    # plot
    
    # Histogramm tree size
    
    plt.figure(num = ii * len(geneNames) * 3 + jj * 3 + 0)
    plt.hist(tree_sizes, bins = 100, color = '#c47e21')
    plt.title("Tree sizes" + "\n" + shortNames[ii] + "  " + shorGeneNames[jj])
    plt.ylabel("Absolut number of trees")
    plt.xlabel("Number of leaves per tree")
    plt.savefig("tree_sizes" + "_" + name + "_" + geneName + ".pdf", format="pdf")
    
    plt.figure(num = ii * len(geneNames) * 3 + jj * 3 + 1, figsize = (6.4 * len(namespace) / 32, 4.8))
    plt.bar(namespace, height = coverage, color = '#c43a21')
    plt.title("Species coverage" + "\n" + shortNames[ii] + "  " + shorGeneNames[jj])
    plt.ylabel("Number of trees")
    plt.xlabel("Species")
    plt.savefig("coverage" + "_" + name + "_" + geneName + ".pdf", format="pdf")
    
    plt.figure(num = ii * len(geneNames) * 3 + jj * 3 + 2, figsize = (6.4 * len(namespace) / 32, 4.8))
    plt.bar(namespace, height = av_cover, color = '#14c9b4')
    plt.title("Average coverage per tree" + "\n" + shortNames[ii] + "  " + shorGeneNames[jj])
    plt.ylabel("Average number of appearances per tree")
    plt.xlabel("Species")
    plt.savefig("av_coverage" + "_" + name + "_" + geneName + ".pdf", format="pdf")
    
