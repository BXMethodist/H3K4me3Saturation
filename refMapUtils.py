import numpy as np, pickle
from Wig import Wig
import os
import wigChrom
import pandas as pd


def genome_size(path="/archive/tmhkxc48/ref_data/hg19/hg19.chrom.sizes.xls", step=10):
    genome = {}
    genome_size_file = open(path, "r")
    for line in genome_size_file.readlines():
        chr_name, chr_size = line.split("\t")
        vector_size = int(chr_size.rstrip())
        if vector_size % step == 0:
            vector_size /= step
        else:
            vector_size = vector_size / step + 1
        genome[chr_name] = np.zeros(vector_size)
    genome_size_file.close()
    return genome


def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    f = open(name + '.pkl', 'rb')
    result = pickle.load(f)
    f.close()
    return result

def coverage_saturation(wig_file_path, start_cutoff, end_cutoff, step_cutoff):
    wig_file = Wig(wig_file_path)
    results = []
    for cutoff in range(start_cutoff, end_cutoff, step_cutoff):
        coverage = wig_file.get_coverage(cutoff)
        results.append((cutoff, coverage))
    return results


def sample_vs_cutoff(cutoffs, wigPath="/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/", partition=None):
    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    if partition is not None:
        block_size = len(wigFiles)/20
        start = partition * block_size
        end = (partition + 1) * block_size
        wigFiles = wigFiles[start:end]

    avg_size = np.zeros((len(cutoffs), len(wigFiles)))
    coverage = np.zeros((len(cutoffs), len(wigFiles)))
    peak_number = np.zeros((len(cutoffs), len(wigFiles)))

    for i in range(len(wigFiles)):
        wigFile = wigFiles[i]
        wig_obj = Wig(wigPath+wigFile)

        for j in range(len(cutoffs)):
            cutoff = cutoffs[j]
            coverage[j, i] = wig_obj.get_coverage(cutoff)
            peak_number[j, i] = wig_obj.get_peak_number(cutoff)
            avg_size[j, i] = coverage[j, i]/peak_number[j, i]

    if partition is None:
        np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/avg_peak_size_vs_cutoff.txt", avg_size, delimiter="\t")
        np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/coverage_vs_cutoff.txt", coverage, delimiter="\t")
        np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/peak_number_vs_cutoff.txt", peak_number, delimiter="\t")
    else:
        np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/avg_peak_size_vs_cutoff_" + str(partition) + ".txt",
                   avg_size, delimiter="\t")
        np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/coverage_vs_cutoff_" + str(partition) + ".txt",
                   coverage, delimiter="\t")
        np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/peak_number_vs_cutoff_" + str(partition) + ".txt",
                   peak_number, delimiter="\t")

    return coverage, peak_number, avg_size

def combine_patition_result(list_files, name, path="./refmap_cutoff/"):
    result_df = pd.read_csv(path+list_files[0], sep="\t", index_col=None, header=None)
    for df_path in list_files[1:]:
        print df_path
        df = pd.read_csv(path+df_path, sep="\t", index_col=None, header=None)
        result_df = pd.concat([result_df, df], axis=1)

    result_df.to_csv(name+".txt", sep="\t")
    return result_df


def combine_wig(wig1, wig2):
    for chr_name, chromosome_wig1 in wig1.genome.items():
        if chr_name not in wig2.genome:
            print chr_name, " is not in this wig!?"
        else:
            chromosome_wig2 = wig2.genome[chr_name]
            chromosome_wig1.signals = chromosome_wig1.signals + chromosome_wig2.signals

    for chr_name, chromosome_wig2 in wig2.genome.items():
        if chr_name not in wig1.genome:
            chromosome_wig2 = wig2.genome[chr_name]
            print chr_name, " is not in the wig1!!!"
            wig1.genome[chr_name] = wigChrom.WigChrom(chr_name, chromosome_wig2.start, chromosome_wig2.size,
                                                      chromosome_wig2.step, chromosome_wig2.span, chromosome_wig2.fixed)
            wig1.genome[chr_name].signals = chromosome_wig2.signals
    return wig1

def super_wig():
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    # superwig = Wig(wigPath+wigFiles[0])

    for i in range(0, len(wigFiles)):
        print wigFiles[i]
        cur_wig = Wig(wigPath+wigFiles[i])
    #     superwig = combine_wig(superwig, cur_wig)
        save_obj(cur_wig, "./wigs/"+cur_wig.file_name)
    return


if __name__ == "__main__":
    # super_wig()
    total_wig_signals = {}
    wigs = [x[:-4] for x in os.listdir('./wigs') if x.endswith('.pkl') and x[:-4] != '']
    # print wigs
    for wig in wigs:
        # wig_pkl = load_obj('./wigs/'+wig)
        # total_signals=0
        # for key, value in wig_pkl.genome.iteritems():
        #     total_signals += np.sum(value.signals)
        # total_wig_signals[wig] = total_signals
        # print total_signals
        total_wig_signals[wig] = 5000000000

    FRIP10 = []
    FRIP100 = []
    path10 = '/home/tmhbxx3/archive/KFH3K4me3/10cutoff/pooled/'
    path100 = '/home/tmhbxx3/archive/KFH3K4me3/100cutoff/pooled/'
    for wig in wigs:
        peaks_xls10 = path10+'archive_tmhkxc48_BroadH3K4me3_broadpeak201401_H3K4me3_dregion_pooled_'+wig+'.peaks.xls'
        peaks_xls100 = path100 + 'archive_tmhkxc48_BroadH3K4me3_broadpeak201401_H3K4me3_dregion_pooled_' + wig + '.peaks.xls'
        df10 = pd.read_csv(peaks_xls10, sep='\t')
        df100 = pd.read_csv(peaks_xls100, sep='\t')
        sum10 = df10['total_signal'].sum()
        sum100 = df100['total_signal'].sum()
        # print sum10, sum100
        # print peaks_xls10
        FRIP10.append((wig, sum10/total_wig_signals[wig]))
        FRIP100.append((wig, sum100 / total_wig_signals[wig]))

        print 'complete', wig

    result_df10= pd.DataFrame(FRIP10)
    result_df100 = pd.DataFrame(FRIP100)

    result_df10.to_csv('FRIP_336_10_cutoff.csv', index=None)

    result_df100.to_csv('FRIP_336_100_cutoff.csv', index=None)
    pass