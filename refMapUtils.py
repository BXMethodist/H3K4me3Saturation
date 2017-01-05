import numpy as np, pickle
from Wig import Wig
import os


def genome_size(path="/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt", step=10):
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
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def coverage_saturation(wig_file_path, start_cutoff, end_cutoff, step_cutoff):
    wig_file = Wig(wig_file_path)
    results = []
    for cutoff in range(start_cutoff, end_cutoff, step_cutoff):
        coverage = wig_file.get_coverage(cutoff)
        results.append((cutoff, coverage))
    return results


def sample_vs_cutoff(cutoffs, wigPath="/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"):
    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    avg_size = np.zeros((len(cutoffs), len(wigFiles)))
    coverage = np.zeros((len(cutoffs), len(wigFiles)))
    peak_number = np.zeros((len(cutoffs), len(wigFiles)))

    for i in range(len(wigFiles)):
        wigFile = wigFiles[i]
        wig_obj = Wig(wigPath+wigFile)
        for cutoff in cutoffs:
            coverage[cutoff, i] = wig_obj.get_coverage(cutoff)
            peak_number[cutoff, i] = wig_obj.get_peak_number(cutoff)
            avg_size[cutoff, i] = coverage[cutoff, i]/peak_number[cutoff, i]

    np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/avg_peak_size_vs_cutoff.txt", avg_size, delimiter="\t")
    np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/coverage_vs_cutoff.txt", coverage, delimiter="\t")
    np.savetxt("/home/tmhbxx3/archive/refmap_saturation/code/peak_number_vs_cutoff.txt", peak_number, delimiter="\t")

    return coverage, peak_number, avg_size


def combine_wig(wig1, wig2):
    pass


def super_wig():
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    superwig = Wig(wigFiles[0])

    for i in range(1, len(wigFiles)):
        cur_wig = Wig(wigFiles[i])
        superwig = combine_wig(superwig, cur_wig)
    save_obj(superwig, "superwig")





if __name__ == "__main__":
    cutoffs = [x for x in range(10,300,10)]
    sample_vs_cutoff(cutoffs)