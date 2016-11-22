import pandas as pd, numpy as np
import os
import csv


def genomeSize(path="./hg19_chr_sizes.txt"):
    genome = {}

    file = open(path, "r")
    for line in file.readlines():
        chrName, chrSize = line.split("\t")
        genome[chrName] = np.zeros(int(chrSize.rstrip()))
    file.close()
    return genome


class H3K4me3Saturation:
    def __init__(self):
        self.genome = genomeSize()
        self.results = []


    def saturated(self, path):
        df = pd.read_excel(path)
        rowNum = df.shape[0]

        for i in range(rowNum):
            start = df.ix[i, :]['start']
            end = df.ix[i, :]['end']
            chrName = df.ix[i, :]['chr']
            self.genome[chrName][start-1:end] = 1

        totalCoverage = 0
        totalIsland = 0
        totalIslandLength = 0

        for value in self.genome.values():
            totalCoverage += np.count_nonzero(value == 1)
            breakpointsPosition = np.where(((np.roll(value, 1) - value) != 0).astype(int)==1)

            if value[0] == 1:
                breakpointsPosition = np.insert(breakpointsPosition, 0, 0)
            if value[-1] == 1:
                breakpointsPosition = np.insert(breakpointsPosition, -1, value.shape[0])
            islandLength = (breakpointsPosition - np.roll(breakpointsPosition, 1))[0][1::2]
            totalIsland += islandLength.shape[0]
            totalIslandLength += np.sum(islandLength)

        avgLength = totalIslandLength*1.0/totalIsland

        self.results.append([totalCoverage, totalIsland, avgLength])


    def converge(self, prev, current, convergeCap = 10):
        prevCoverage, prevNumber = prev
        curCoverage, curNumber = current

        converge = 0
        if prevCoverage*1.1 > curCoverage and (prevNumber*0.9<curNumber):
            converge += 1
        else:
            converge = 0

        return converge > convergeCap


    def trainReferenceMap(self, directories):
        listFiles = os.listdir(directories)
        listFiles = np.random.shuffle(listFiles)

        for file in listFiles:
            self.saturated(file)

        output = open("H3K4me3.csv", "wb")
        writer = csv.writer(output)
        for row in self.results:
            writer.writerow(row)

        output.close()








