import numpy as np, pandas as pd
import matplotlib.pyplot as plt

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
#     break
#
# import matplotlib.pyplot as plt
# print key
#
# plt.plot(np.arange(100), v)
# plt.show()
