import numpy as np
from Wig import Wig


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


def coverage_saturation(wig_file_path, start_cutoff, end_cutoff, step_cutoff):
    wig_file = Wig(wig_file_path)
    results = []
    for cutoff in range(start_cutoff, end_cutoff, step_cutoff):
        coverage = wig_file.get_coverage(cutoff)
        results.append((cutoff, coverage))
    return results




if __name__ == "__main__":
    print coverage_saturation("", 5, 501, 5)