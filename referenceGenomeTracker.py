import numpy as np
import os
import csv
import matplotlib.pyplot as plt


def genomeSize(path="./hg19_chr_sizes.txt"):
    genome = {}

    file = open(path, "r")
    for line in file.readlines():
        chrName, chrSize = line.split("\t")
        size = int(chrSize.rstrip())
        if size%10 == 0:
            size = size/10
        else:
            size = size/10+1
        genome[chrName] = np.zeros(size)
    file.close()
    return genome


def plotConverge(table, xaxis, yaxis_start, yaxis_end,  xTitle, yTitle, Title, cutoff=0, color="r-"):
    table = np.asarray(table)
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.plot(table[:, xaxis], table[:, yaxis_start:yaxis_end], color, linewidth=2.0)
    plt.ylabel(yTitle, fontname="Times New Roman")
    plt.xlabel(xTitle, fontname="Times New Roman")
    plt.title(Title +" with cutoff" + str(cutoff), fontname="Times New Roman")
    Title = Title.replace(" ", "_")
    plt.savefig(Title + "with_cutoff" + str(cutoff), dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None,
            frameon=None)



class H3K4me3Saturation:
    def __init__(self, iterations):
        self.genome = genomeSize()
        self.iterations = iterations
        self.coverage = None
        self.region = None
        self.regionLength = None
        self.numberSample = None


    def initialization(self, sampleNumber):
        self.coverage = [[i+1] for i in range(sampleNumber)]
        self.region = [[i+1] for i in range(sampleNumber)]
        self.regionLength = [[i + 1] for i in range(sampleNumber)]


    def saturated(self, path, sampleSequence, cutoff=0):

        #really need pandas to read heteregous data table.

        # df = pd.read_excel(path)
        # rowNum = df.shape[0]
        #
        # for i in range(rowNum):
        #     start = df.ix[i, :]['start']/10
        #     end = df.ix[i, :]['end']/10
        #     chrName = df.ix[i, :]['chr']
        #     self.genome[chrName][start-1:end] = 1

        file = open(path, "rb")

        for line in file.readlines():
            info = line.split("\t")
            if info[1] == "start":
                continue
            start = int(info[1])/10
            end = int(info[2])/10
            height = float(info[6])
            chrName = info[0]
            if height >= cutoff:
                if chrName in self.genome:
                    self.genome[chrName][start-1:end] = 1
        totalCoverage = 0
        totalIsland = 0

        for value in self.genome.values():
            newCoverage = np.sum(value)
            totalCoverage += newCoverage
            if newCoverage == 0:
                continue
            else:
                sign = ((np.roll(value, 1) - value) != 0).astype(int)
                sign[0] = 0
                islandNumber = np.sum(sign)
                if value[0] == 1:
                    islandNumber+=1
                if value[-1] == 1:
                    islandNumber+=1
                totalIsland += islandNumber/2
        if totalIsland == 0:
            avgLength = 0
        else:
            avgLength = totalCoverage*1.0/totalIsland*10

        self.coverage[sampleSequence].append(totalCoverage*10)
        self.region[sampleSequence].append(totalIsland)
        self.regionLength[sampleSequence].append(avgLength)


    def converge(self, prev, current, convergeCap = 10):
        prevCoverage, prevNumber = prev
        curCoverage, curNumber = current

        converge = 0
        if prevCoverage*1.1 > curCoverage and (prevNumber*0.9<curNumber):
            converge += 1
        else:
            converge = 0

        return converge > convergeCap


    def reset(self):
        self.genome = genomeSize()


    def draw(self, cutoff):

        plotConverge(self.coverage, 0, 1, self.iterations+1 , "Number of Sample", "Coverage",
                     "H3K4me3 Coverage VS Number of Sample", cutoff)

        plotConverge(self.region, 0, 1, self.iterations+1, "Number of Sample", "Region Number",
                     "H3K4me3 Region Number VS Number of Sample", cutoff)

        plotConverge(self.regionLength, 0, 1, self.iterations+1, "Number of Sample", "Region Length",
                     "H3K4me3 Region Length VS Number of Sample", cutoff)




    def trainMap(self, directories, cutoff=0):
        listFiles = os.listdir(directories)

        self.numberSample = len(listFiles)

        self.initialization(self.numberSample)

        n = 0

        while n <=self.iterations:
            n+=1
            np.random.shuffle(listFiles)
            seq = 0
            for file in listFiles:
                self.saturated(directories+'/'+file, seq, cutoff=cutoff)
                seq += 1

            output = open(str(n)+"H3K4me3.csv", "wb")
            writer = csv.writer(output)
            for i in range(self.numberSample):
                writer.writerow([self.coverage[i][0],
                                 self.coverage[i][n],
                                 self.region[i][n],
                                 self.regionLength[i][n]])

            output.close()

            self.reset()

        self.draw(cutoff)












