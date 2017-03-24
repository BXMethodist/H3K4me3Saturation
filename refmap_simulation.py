import numpy as np, pandas as pd, os
import matplotlib.pyplot as plt
from Wig import Wig
from refMapUtils import save_obj


def callpeak(sampleID, genome, cutoff):
    colnames = ['chr', 'start', 'end', 'center', 'width_above_cutoff', 'total_signal', 'height', 'height_logP']
    peaks = []
    for key, value in genome.iteritems():
        cur_chromosome = np.copy(value)
        peak = []
        cur_chromosome[cur_chromosome-cutoff<=0] = 0
        cur_chromosome[cur_chromosome-cutoff>0] = 1
        start = None
        end = None
        for i in range(cur_chromosome.shape[0]):
            if cur_chromosome[i] == 1 and start is not None:
                continue
            elif cur_chromosome[i] == 0 and start is not None:
                end = i
                peaks.append([key, start, end, 0, 0, 0, 0, 0])
                start  = None
                end = None
            elif cur_chromosome[i] == 0 and start is None:
                continue
            elif cur_chromosome[i] == 1 and start is None:
                start = i
        if start is not None:
            end = i
            if start != end:
                peaks.append([key, start, end, 0, 0, 0, 0, 0])
            start = None
            end = None
    df = pd.DataFrame(peaks, index=None, columns=colnames)
    df.to_csv(sampleID+str(cutoff)+'.csv', sep='\t', header=True, index=None)

def callcoverage(wig, cutoffs):
    """
    :param wig: wigfile object
    :param cutoff: a list contains different cutoff values
    :return: a dataframe, index is cutoffs, columns are 'noise_coverage', 'noise_signals,
    'non_noise_coverage', 'non_noise_signals', (coverage * step is equals to the original ones'
    """
    array = np.zeros((len(cutoffs), 4))

    df = pd.DataFrame(array, index=cutoffs,
                      columns=['noise_coverage', 'noise_signals', 'non_noise_coverage', 'non_noise_signals'])

    for i in range(len(cutoffs)):
        cutoff = cutoffs[i]
        for value in wig.genome.values():
            cur_non_noise = value.signals[value.signals>=cutoff]
            cur_non_zero = value.signals[value.signals>0]
            df.ix[cutoff, 2] += len(cur_non_noise)
            df.ix[cutoff, 0] += len(cur_non_zero) - len(cur_non_noise)
            cur_non_noise_signals = np.sum(cur_non_noise)
            df.ix[cutoff, 1] += np.sum(cur_non_zero) - cur_non_noise_signals
            df.ix[cutoff, 3] += cur_non_noise_signals
    df.to_csv(wig.file_name+'.csv')
    return

def simulation():
    f = open('hg19_genome_size.txt', 'r')

    info = f.readlines()

    genome = {}

    f.close()
    table = []
    for line in info:
        line = line.strip()
        line = line.split()
        line[1] = int(line[1])/10
        genome[line[0]] = np.zeros(line[1])
        table.append(line)


    genome_size = pd.DataFrame(table)

    total_length = genome_size.iloc[:, 1].sum()

    total_reads_number = 25000000

    cutoffs = [x for x in range(2, 20, 2)]

    for iteration in range(300):
        for key in genome.keys():
            # genome[key] =np.zeros(100)
            genome[key] = np.zeros(genome[key].shape[0])
        for key, v in genome.iteritems():
            cur_chrome_len = v.shape[0]
            # cur_reads_number = 30
            cur_reads_number = int(total_reads_number*(cur_chrome_len*1.0/total_length))
            cur_indexes = np.random.randint(0, cur_chrome_len-20, size=cur_reads_number)
            for i in cur_indexes:
                v[i:i+20] += 1
        for cutoff in cutoffs:
            sampleID = './simulation_results/simulated_sample'+str(iteration)
            callpeak(sampleID, genome, cutoff)


def noise_simulation(cutoffs):
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    n = 0
    for i in range(len(wigFiles))[300:]:
        print wigFiles[i]
        cur_wig = Wig(wigPath+wigFiles[i])
        callcoverage(cur_wig, cutoffs)
        n+=1
    return n

noise_simulation([x for x in range(10,310, 10)])