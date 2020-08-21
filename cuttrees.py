import random

specs = ["10", "20"]
trees = ["25", "100", "250"]
rates = ["0", "1.9e-10", "2.7e-10", "4.9e-10", "5.2e-11"]

names = ["baum.newick"]
#for tree in trees:
#  for spec in specs:
#    for rate in rates:
#      names.append(tree + "_" + spec + "_" + rate + "_" + rate + "_470000000")

for name in names:
  for i in range(1,2):
    num = str(i)
    if i < 10:
      num = "0" + num
    newTree = ""
    geneTrees = open(name, "r").read()
#    geneTrees = open(name + "/" + num + "/geneTrees.newick", "r").read()
    i = 0
    for geneTree in geneTrees.split(";"):
      i += 1
      appears = geneTree.count(")")
      if appears == 0:
        continue
      toRemove = random.randint(1, appears)
      isRoot = (toRemove == appears)
      num = 0
      index = 0
      for j, char in enumerate(geneTree):
        if char == ')':
          num += 1
        if num == toRemove:
          index = j
          break
      braketCount = 0
      commaCount = 0
      while True:
        if geneTree[index] == ')':
          braketCount += 1
        if geneTree[index] == '(':
          braketCount -= 1
          if braketCount == 0:
            break
        if (geneTree[index] == ',' and braketCount == 1):
          commaCount += 1
        index -= 1
      toRemove = random.randint(0, commaCount)
      
      if (commaCount > 1):
        commaCount = 0
        braketCount = 0
        while True:
          if geneTree[index] == ',' and braketCount == 1:
            commaCount += 1
          if geneTree[index] == '(':
            braketCount += 1
          if geneTree[index] == ')':
            braketCount -= 1
            if braketCount == 0:
              break
          if (commaCount == toRemove or (commaCount == 1 and geneTree[index] == ',' and toRemove == 0 and braketCount == 1)):
            geneTree = geneTree[:index] + geneTree[index + 1:]
          else:
            index += 1
      else:
        commaCount = 0
        braketCount = 1
        geneTree = geneTree[:index] + geneTree[index + 1:]
        while True:
          if (geneTree[index] == ',' and braketCount == 1):
            commaCount += 1
          if geneTree[index] == '(':
            braketCount += 1
          if geneTree[index] == ')':
            braketCount -= 1
            if braketCount == 0:
              geneTree = geneTree[:index] + geneTree[index + 1:]
              if not isRoot:
                while True:
                  if geneTree[index] == ',' or geneTree[index] == ')':
                    break
                  geneTree = geneTree[:index] + geneTree[index + 1:]
              break
          if (commaCount == toRemove or (braketCount == 1 and geneTree[index] == ',')):
            geneTree = geneTree[:index] + geneTree[index + 1:]
          else:
            index += 1
      
      newTree += geneTree + ";"
    
    geneTrees_new = open(name + "_cut.newick", "w")
    geneTrees_new.write(newTree)
    geneTrees_new.close()
