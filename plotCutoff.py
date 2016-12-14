import numpy as np
import matplotlib.pyplot as plt
import os

def plotNumpyArray(table, Title):
    print table
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(table[:, 0], table[:, 1:12], linewidth=2.0)
    plt.legend([str(n)+" cutoff" for n in [10,25,50,75,100,150,200,250,300,400,500]], loc=2)
    plt.savefig("/home/tmhbxx3/archive/"+Title, dpi=600)

    return


path = os.getcwd()

cutoffFiles = [file for file in os.listdir(path) if file.find("cutoff") != -1]

candidates = {10:1, 25:2, 50:3, 75:4, 100:5, 150:6, 200:7, 250:8, 300:9, 400:10, 500:11}

file = open(cutoffFiles[0], "r")
count = 0
for line in file.readlines():
    count+=1

file.close()

coverage = np.zeros((count, 12))
coverage[:, 0] = np.arange(1, count+1)

regionLength = np.zeros((count, 12))
regionLength[:, 0] = np.arange(1, count+1)

region = np.zeros((count, 12))
region[:, 0] = np.arange(1, count+1)

for f in cutoffFiles:
    cutoff = int(f[:f.find("cutoff")])
    colNumber = candidates[cutoff]
    table = np.genfromtxt(f, delimiter=',')
    coverage[:colNumber] = table[:, 1]
    regionLength[: colNumber] = table[:, 2]
    region[: colNumber] = table[:, 3]

plotNumpyArray(coverage, "coverage")
plotNumpyArray(regionLength,"region length")
plotNumpyArray(region, "region")

# table = np.genfromtxt("foo.csv", delimiter=',')
# plotNumpyArray(table)