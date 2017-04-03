import numpy as np
import os
import csv
import matplotlib.pyplot as plt
from refMapUtils import *

class refMap:
    ### take a directory of wig files, generate reference map
    ### based on number of iterations, generate saturation map for coverage, average peak size and number of average peak number
    ### use plotSaturation.py to make figure
    def __init__(self, iterations, genome_size_path="/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt"):
        self.genome = genome_size(genome_size_path)
        self.non_noise_genome = genome_size(genome_size_path)
        self.noise_genome = genome_size(genome_size_path)
        self.iterations = iterations
        self.coverage = None
        self.region = None
        self.regionLength = None
        self.numberSample = None
        self.genome_size_path=genome_size_path

    def initialization(self, sampleNumber):
        self.coverage = [[0] * self.iterations for i in range(sampleNumber)]
        self.region = [[0] * self.iterations for i in range(sampleNumber)]
        self.regionLength = [[0] * self.iterations for i in range(sampleNumber)]

        self.noise_coverage = [[0] * self.iterations for i in range(sampleNumber)]
        self.noise_region = [[0] * self.iterations for i in range(sampleNumber)]
        self.noise_regionLength = [[0] * self.iterations for i in range(sampleNumber)]

        self.non_noise_coverage = [[0] * self.iterations for i in range(sampleNumber)]
        self.non_noise_region = [[0] * self.iterations for i in range(sampleNumber)]
        self.non_noise_regionLength = [[0] * self.iterations for i in range(sampleNumber)]

    def saturated(self, path, sampleSequence, iteration, cutoff=0):
        file = open(path, "rb")

        for line in file.readlines():
            info = line.split("\t")
            if info[1] == "start":
                continue
            start = int(info[1])/10
            end = int(info[2])/10
            height = float(info[6])
            # print height
            chrName = info[0]
            if height >= cutoff:
                if chrName in self.genome:
                    self.genome[chrName][start - 1:end] = 1
                    self.non_noise_genome[chrName][start-1:end] = 1
            elif height < cutoff:
                if chrName in self.genome:
                    self.noise_genome[chrName][start-1:end] = 1
        file.close()
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

        self.coverage[sampleSequence][iteration] = totalCoverage * 10
        self.region[sampleSequence][iteration] = totalIsland
        self.regionLength[sampleSequence][iteration] = avgLength

        totalCoverage = 0
        totalIsland = 0
        for value in self.non_noise_genome.values():
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

        self.non_noise_coverage[sampleSequence][iteration] = totalCoverage * 10
        self.non_noise_region[sampleSequence][iteration] = totalIsland
        self.non_noise_regionLength[sampleSequence][iteration] = avgLength

        totalCoverage = 0
        totalIsland = 0
        for value in self.noise_genome.values():
            newCoverage = np.sum(value)
            totalCoverage += newCoverage
            if newCoverage == 0:
                continue
            else:
                sign = ((np.roll(value, 1) - value) != 0).astype(int)
                sign[0] = 0
                islandNumber = np.sum(sign)
                if value[0] == 1:
                    islandNumber += 1
                if value[-1] == 1:
                    islandNumber += 1
                totalIsland += islandNumber / 2

        if totalIsland == 0:
            avgLength = 0
        else:
            avgLength = totalCoverage * 1.0 / totalIsland * 10

        self.noise_coverage[sampleSequence][iteration] = totalCoverage * 10
        self.noise_region[sampleSequence][iteration] = totalIsland
        self.noise_regionLength[sampleSequence][iteration] = avgLength

    def converge(self, prev, current, convergeCap = 10):
        prevCoverage, prevNumber = prev
        curCoverage, curNumber = current

        converge = 0
        if prevCoverage * 1.1 > curCoverage and (prevNumber * 0.9 < curNumber):
            converge += 1
        else:
            converge = 0

        return converge > convergeCap


    def reset(self):
        self.genome = genome_size(self.genome_size_path)
        self.non_noise_genome = genome_size(self.genome_size_path)
        self.noise_genome = genome_size(self.genome_size_path)


    def saveRefMap(self, cutoff):
        refmap = {}
        for chr, vector in self.genome.items():
            if vector[0] == 1:
                vector = np.insert(vector, 0, 0)
            if vector[-1] == 1:
                vector = np.append(vector, 0)

            # sign change mark the start of the peak and the end to the peak, the end mark is exclusive
            # and this is the index of vector with step = 10, not real genome position
            signchange = ((np.roll(vector, 1) - vector) != 0).astype(int)
            peaksindex = np.where(signchange == 1)[0]

            rowNumber = peaksindex.shape[0] / 2
            colNumber = 2
            peaksindex = peaksindex.reshape((rowNumber, colNumber))

            refmap[chr] = peaksindex

        output = open(str(cutoff) + "_refmap.csv", "w")
        writer = csv.writer(output)
        for chr, index in refmap.items():
            output.write(">" + chr + "\n")
            for i in range(index.shape[0]):
                writer.writerow(index[i, :])
        output.close()


    def trainMap(self, directories, outputname, preffix= None, surfix=None, cutoff=0, saveRefMap=True, individual=False):
        """
        :param directories: directory of danpas peak calling results xlsx(txt format) files
        :param surfix: file suffix
        :param cutoff: peak height cutoff
        :param saveRefMap: whether save the reference map
        :param individual: whether it is a individual value, or accumulative value
        :return: the result table file address, the table index is the sample number, columns are coverage, average length,
                number of peaks
        """
        if individual:
            outputname += "_individual"

        listFiles = os.listdir(directories)

        if not directories.endswith("/"):
            directories += "/"

        if surfix is not None:
            listFiles = [x for x in listFiles if x.endswith(surfix)]

        if preffix is not None:
            listFiles = [x for x in listFiles if x.startswith(preffix)]

        self.numberSample = len(listFiles)

        self.initialization(self.numberSample)

        n = 0

        while n < self.iterations:
            if self.iterations > 1:
                np.random.shuffle(listFiles)
            seq = 0
            for file in listFiles:
                if individual:
                    self.reset()
                self.saturated(directories + file, seq, n, cutoff=cutoff)
                seq += 1

            if saveRefMap and n == 0:
                self.saveRefMap(cutoff)

            self.reset()

            n += 1

        table = np.zeros((self.numberSample, 4))
        table[:, 0] = np.arange(1, self.numberSample + 1)
        table[:, 1] = np.mean(self.coverage, axis=1)
        table[:, 2] = np.mean(self.regionLength, axis=1)
        table[:, 3] = np.mean(self.region, axis=1)

        np.savetxt(outputname + str(cutoff) +".csv", table, delimiter=",")

        return outputname + str(cutoff) +".csv"

