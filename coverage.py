import dendropy
import numpy as np

names = ["cyano_empirical"]
speciesName = "species_mapped"
geneNames = ["geneTrees_mapped.newick"]

for name in names:
  for geneName in geneNames:
    speciesTree = open(name + "/" + speciesName, "r").read()
    geneTrees = open(name + "/" + geneName, "r").read()

    dataSet = dendropy.DataSet()
    dataSet.read(data=speciesTree, schema = "newick")
    namespace = [""] * len(dataSet.tree_lists[0].taxon_namespace)

    for i in range(len(namespace)):
      namespace[i] = str(dataSet.tree_lists[0].taxon_namespace[i])[1:-1]

    coverage = np.zeros(len(namespace))
    crossCoverage = np.zeros(shape = (len(namespace), len(namespace)))

    size = 0
    for geneTree in geneTrees.split(";"):
      size += 1
      for i in range(len(namespace)):
        if (geneTree.find(namespace[i] + ":") != -1):
          coverage[i] += 1
          for j in range(len(namespace)):
            if (geneTree.find(namespace[j] + ":") != -1):
              crossCoverage[i][j] += 1

    size -= 1
    print(size)
    out_str = ""
    for i in range(len(namespace)):
      out_str += namespace[i] + "\t"
    out_str + "\n"
    out_str + "\n"
    for i in range(len(namespace)):
      out_str += "{:.3f}".format(coverage[i] / size) + "\t"
    out_str + "\n"
    out_str + "\n"
    for i in range(len(coverage)):
      for j in range(len(coverage)):
        out_str += "{:.3f}".format((crossCoverage[i][j] / crossCoverage[i][i])) + "\t"
        
      out_str += "\n"
    print(out_str)
