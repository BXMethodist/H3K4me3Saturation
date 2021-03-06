from clusterUtils import *
import os
import wigChrom


class Wig:
    ## This is the class for input a list Wig files and generate the partitions for the wig files
    def __init__(self, wig_file, genome_size_path="/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt"):
        # address of genomesize file on server "/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt"
        self.genome_size = genome_size_chrom(genome_size_path)
        self.genome = {}
        self.initiate(wig_file)
        self.splitted_chroms = {}

        wig_file_name = wig_file.split("/")[-1]

        self.file_name = wig_file_name[:-4]

    def initiate(self, wig_file):
        wig = open(wig_file, "r")
        chr_name = None
        start = None
        step = None
        span = None

        cur_position = 0

        for line in wig.readlines():
            cur_line = line.rstrip().split()
            if len(cur_line) > 1:
                if cur_line[1].startswith("chrom="):
                    chr_name = cur_line[1][cur_line[1].find("=")+1:]
                if cur_line[2].startswith("start="):
                    start = int(cur_line[2][cur_line[2].find("=")+1:])
                if cur_line[3].startswith("step="):
                    step = int(cur_line[3][cur_line[3].find("=") + 1:])
                if cur_line[4].startswith("span="):
                    span = int(cur_line[4][cur_line[4].find("=") + 1:])
                if chr_name in self.genome_size:
                    size = self.genome_size[chr_name]
                    self.genome[chr_name] = wigChrom.WigChrom(chr_name, start, size, step, span)
                else:
                    chr_name = "unknown"
                cur_position = 0 + start - 1
            else:
                if chr_name == "unknown":
                    continue
                if cur_position < self.genome[chr_name].signals.shape[0]:
                    self.genome[chr_name].signals[cur_position] = float(cur_line[0])
                cur_position += 1
        wig.close()

    def split_chr(self, chr_name, split_vector_size=10000):
        # split chromesome based on the input size, this vector size, not the actual chromesome position
        # actual chromesome position = step * vector_size
        chromosome = self.genome[chr_name]

        self.splitted_chroms[chr_name] = {}

        if chromosome.signals.shape[0] % split_vector_size == 0:
            number_split = chromosome.signals.shape[0]/split_vector_size
        else:
            number_split = chromosome.signals.shape[0]/split_vector_size + 1
        for i in range(number_split):
            start = i*split_vector_size
            end = (i+1)*split_vector_size #if (i+1)*split_vector_size < chromosome.signals.shape[0] else -1

            self.splitted_chroms[chr_name][start] = chromosome.signals[start:end]

    def save_split_chr(self, split_vector_size, path="/home/tmhbxx3/archive/H3K4me3_Encode_wig_splits/"):
        # save the splited array to corresponding folder, the first line of file will always be
        # fixedStep chrom=chrX start=1  step=10 span=10
        # create map size for splited file
        # this is where I store the file on sever path="/home/tmhbxx3/archive/WigChrSplits/"
        if not path.endswith("/"):
            path += "/"

        for chr_name, chromosome in self.genome.items():
            self.split_chr(chr_name, split_vector_size)
            output_path = path + chr_name + "_"
            # print self.splitted_chroms

            for key, value in self.splitted_chroms[chr_name].items():
                vector_size = value.shape[0]
                start = str(key*chromosome.step)
                end = str(key*chromosome.step + vector_size*chromosome.step)
                step = str(chromosome.step)
                final_output_path = output_path + start + "_" + end + "_" + step + "/"

                if not os.path.isdir(final_output_path):
                    os.system("mkdir " + final_output_path)

                output_file_name = final_output_path + self.file_name + "_" + chr_name + "_" + start \
                                        + "_" + end + "_" + step + ".txt"
                np.savetxt(output_file_name, value, delimiter=",", fmt='%1.1f')

    def get_coverage(self, cutoff):
        coverage = 0
        for chr_name, chromosome in self.genome.items():
            coverage += (chromosome.signals>=cutoff).sum() * chromosome.step
        return coverage

    def get_avg_peak_size(self, cutoff):
        coverage = self.get_coverage(cutoff)
        peak_number = self.get_peak_number(cutoff)
        return coverage/peak_number

    def get_peak_number(self, cutoff):
        peak_number = 0
        for chr_name, chromosome in self.genome.items():
            peak_map = np.zeros(chromosome.signals.shape[0])
            peak_map[chromosome.signals>=cutoff] = 1
            sign = ((np.roll(peak_map, 1) - peak_map) != 0).astype(int)
            sign[0] = 0
            islandNumber = np.sum(sign)
            if peak_map[0] == 1:
                islandNumber += 1
            if peak_map[-1] == 1:
                islandNumber += 1
            peak_number += islandNumber / 2
        return peak_number







