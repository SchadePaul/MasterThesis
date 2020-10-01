import os
import matplotlib.pyplot as plt

colors = ['#7b14c9', '#1d14c9', '#14c9b4', '#23c914', '#faf61e', '#c47e21', '#c43a21']
dir = "."
files = [o for o in os.listdir(dir) if os.path.isfile(os.path.join(dir,o))]
files = sorted(files)

file = 0
for filename in files:
    print(filename)
    if not "25_1000_4.9e-10_4.9e-10_470000000_cut" in filename:
        continue
    treesizes = []
    coverage = []
    av_coverage = []
    fileRes = open(filename, "r").read()
    for i,line in enumerate(fileRes.split("\n")):
        if line == "":
            continue
        if i == 0:
            for size in line.split("\t"):
                if not size == "":
                    treesizes.append(int(size))
        else:
            coverage.append(0.0)
            av_coverage.append(0.0)
            for j,cov in enumerate(line.split("\t")):
                if j == 0:
                    continue
                if not cov == "":
                    if j%2 == 1:
                        coverage[i - 1] += float(cov)
                    else:
                        av_coverage[i - 1] += float(cov)
            av_coverage[i - 1] = av_coverage[i - 1] / 50.0
            
    plt.figure(num = file * 3 + 0)
    plt.hist(treesizes, bins = 50, range=[0, 250], color = [colors[file % len(colors)]])
    plt.ylim(0,3500)
    plt.title("Tree sizes")
    plt.ylabel("Absolut number of trees")
    plt.xlabel("Number of leaves per tree")
    plt.savefig("tree_sizes_cut_" + str(file) + ".pdf", format="pdf")

    plt.figure(num = file * 3 + 1)
    plt.bar(range(len(coverage)), height = coverage, color = [colors[file % len(colors)]])
    plt.ylim(0,17000)
    plt.xticks(range(len(coverage)),["" for o in range(len(coverage))])
    plt.title("Species coverage")
    plt.ylabel("Number of trees")
    plt.xlabel("Species")
    plt.savefig("coverage_cut_" + str(file) + ".pdf", format="pdf")
    
    plt.figure(num = file * 3 + 2)
    plt.bar(range(len(av_coverage)), height = av_coverage, color = [colors[file % len(colors)]])
    plt.ylim(0,3.7)
    plt.xticks(range(len(coverage)),["" for o in range(len(coverage))])
    plt.title("Average coverage per tree")
    plt.ylabel("Average number of appearances per tree")
    plt.xlabel("Species")
    plt.savefig("av_coverage_cut_" + str(file) + ".pdf", format="pdf")
    
    plt.close('all')
    
    file += 1
