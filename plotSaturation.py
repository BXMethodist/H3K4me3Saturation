import matplotlib.pyplot as plt
import numpy as np





def plotCoverageToSample(table, cutoff=0):
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(table[:, 0], table[:, [2,3]], 'r-', linewidth=2.0)
    plt.ylabel("Coverage of H3K4me3", fontname="Times New Roman")
    plt.xlabel("Number of H3K4me3 Chip-Seq Sample", fontname="Times New Roman")
    plt.title("Coverage VS Sample Number", fontname="Times New Roman")
    plt.savefig("./CoverageToSample", dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None,
            frameon=None)


def plotRegionToSample(table, cutoff=0):
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(table[:, 0], table[:, 2], 'g-', linewidth=2.0)
    plt.ylabel("Number of H3K4me3 Region", fontname="Times New Roman")
    plt.xlabel("Number of H3K4me3 Chip-Seq Sample", fontname="Times New Roman")
    plt.title("Region Number VS Sample Number", fontname="Times New Roman")
    plt.savefig("./RegionToSample", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None,
                frameon=None)


def plotLengthToSample(table, cutoff=0):
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(table[:, 0], table[:, 3], 'b-', linewidth=2.0)
    plt.ylabel("Length of H3K4me3 Region", fontname="Times New Roman")
    plt.xlabel("Number of H3K4me3 Chip-Seq Sample", fontname="Times New Roman")
    plt.title("Region Length VS Sample Number", fontname="Times New Roman")
    plt.savefig("./LengthToSample", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None,
                frameon=None)




my_data = np.genfromtxt('H3K4me3Together.csv', delimiter=',')

plotCoverageToSample(my_data)
plotRegionToSample(my_data)
plotLengthToSample(my_data)

