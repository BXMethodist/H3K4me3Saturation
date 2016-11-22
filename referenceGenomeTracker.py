# import pandas as pd
import numpy as np
import os
import csv


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


class H3K4me3Saturation:
    def __init__(self):
        self.genome = genomeSize()
        self.results = []


    def saturated(self, path):

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
            chrName = info[0]
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
                islandNumber = np.sum(((np.roll(value, 1) - value) != 0).astype(int))
                if value[0] == 1:
                    islandNumber+=1
                if value[-1] == 1:
                    islandNumber+=1
                totalIsland += islandNumber/2
        if totalIsland == 0:
            avgLength = 0
        else:
            avgLength = totalCoverage*1.0/totalIsland*10

        self.results.append([totalCoverage*10, totalIsland, avgLength])


    def converge(self, prev, current, convergeCap = 10):
        prevCoverage, prevNumber = prev
        curCoverage, curNumber = current

        converge = 0
        if prevCoverage*1.1 > curCoverage and (prevNumber*0.9<curNumber):
            converge += 1
        else:
            converge = 0

        return converge > convergeCap


    def trainMap(self, directories):
        listFiles = os.listdir(directories)

        n = 0

        while n <=0:
            n+=1
            np.random.shuffle(listFiles)

            for file in listFiles:
                self.saturated(directories+'/'+file)

            output = open(n+"H3K4me3.csv", "wb")
            writer = csv.writer(output)
            for row in self.results:
                writer.writerow(row)

            output.close()

            self.genome = genomeSize()








