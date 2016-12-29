import numpy as np, os, pandas as pd
from collections import defaultdict


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


def genome_size_chrom(path="/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt"):
    genome = {}
    genome_size_file = open(path, "r")
    for line in genome_size_file.readlines():
        chr_name, chr_size = line.split()
        genome[chr_name] = int(chr_size)
    genome_size_file.close()
    return genome

def get_split_chr(chr_name, start, end, vector_size=10000, step=10, merge=10, cutoff=100,
                  delimiter="\t", path="/home/tmhbxx3/archive/WigChrSplits"):
    if not path.endswith("/"):
        path += "/"

    if abs(start - end) >= vector_size*step*2:
        return

    if not os.path.isdir("./csv"):
        os.system("mkdir csv")

    output_path = path + chr_name + "_"

    unit_size = vector_size * step

    start_name = start - start % unit_size
    end_name = (end/unit_size)*unit_size + unit_size if end % unit_size != 0 else end

    second_file = None

    if end_name - start_name == unit_size:
        first_file = output_path + str(start_name) + "_" + str(end_name) + "_" + str(step)+ "/"
    else:
        first_file = output_path + str(start_name) + "_" + str(start_name + unit_size) + "_" + str(step) + "/"
        second_file = output_path + str(start_name + unit_size) + "_" + str(end_name) + "_" + str(step) + "/"

    print first_file
    print second_file

    sample_number = len(os.listdir(first_file))

    sample_size = (end - start)/step if (end - start) % step == 0 else (end - start)/step + 1

    result_array = np.zeros((sample_number, sample_size))

    if second_file is None:
        start_index = (start - start_name)/step
        end_index = (end - start_name)/step if (end - start_name) % step == 0 else (end - start_name)/step + 1

        file_names = os.listdir(first_file)
        for i in range(len(file_names)):
            cur_array = np.genfromtxt(first_file + file_names[i])
            result_array[i, :] = cur_array[start_index:end_index]
    else:
        start_index = (start - start_name) / step
        mid_index = vector_size - start_index
        end_index = (end - start_name - unit_size)/step if (end - start_name - unit_size) % step == 0 \
            else (end - start_name - unit_size)/step + 1

        file_names = os.listdir(first_file)
        for i in range(len(file_names)):
            cur_array = np.genfromtxt(first_file + file_names[i])

            result_array[i, 0:mid_index] = cur_array[start_index:]

        file_names = os.listdir(second_file)
        for i in range(len(file_names)):
            cur_array = np.genfromtxt(second_file + file_names[i])
            result_array[i, mid_index:] = cur_array[0:end_index]

    # fileter out the data below the cutoff
    copy_result = result_array.copy()

    copy_result[copy_result<=cutoff] = 0
    row_to_delete = []
    for i in range(copy_result.shape[0]):
        if np.sum(copy_result[i, :]) == 0:
            row_to_delete.append(i)

    result_array = np.delete(result_array, row_to_delete, 0)

    sample_names = [file_names[i][:file_names[i].find("_")] for i in range(len(file_names)) if i not in row_to_delete]

    #TO DO, use merge parameter to support generate smaller array if want larger step.

    ####

    df = pd.DataFrame(result_array, index=sample_names)

    output_elements = [chr_name, str(start), str(end)]
    output = "_".join(output_elements)

    df.to_csv("./csv/" + output + ".csv", sep=delimiter)

    return result_array

def quantile_normalization(array, axis=0):
    array = np.asarray(array)
    if axis == 0:
        array = array.T
    elif axis == 1:
        pass
    else:
        return
    df = pd.DataFrame(array)
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    quantile_normalized = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    result = quantile_normalized.as_matrix()

    if axis == 0:
        return result.T
    elif axis == 1:
        return result

def mean_normalization(array, axis=1, target_signal=200):
    if array.ndim == 1:
        print "please convert to samples as array, shape (n_samples, ) which is 2d array"
        return

    n_sample = array.shape[0]
    factors = target_signal/(np.mean(array, axis=axis)).reshape(n_sample, 1)

    return np.multiply(array, factors)


def smooth_normalization(array, smooth=20, axis=1):
    df = pd.DataFrame(array)
    return df.rolling(min_periods=1, window=smooth, center=True, axis=axis).mean().values

def max_normalization(array, max_peak=500, axis=1):
    if array.ndim == 1:
        print "please convert to samples as array, shape (n_samples, ) which is 2d array"
        return

    n_sample = array.shape[0]
    factors = max_peak / (np.amax(array, axis=axis)).reshape(n_sample, 1)
    return np.multiply(array, factors)

def peak_combiner(peak_ref_map, combine_cut_off, step=10):
    combine_cut_off = combine_cut_off/step

    input_file = open(peak_ref_map, "r")

    region_map = defaultdict(list)

    start = None
    end = None
    chr_name = None
    for line in input_file.readlines():
        if line.startswith(">"):
            if start is not None:
                region_map[chr_name].append((start, end))
            chr_name = line.rstrip()[1:]
            start = None
            end = None
        else:
            line = line.rstrip().split(",")
            x, y = int(line[0]), int(line[1])
            if start is None:
                start = x
                end = y
            elif x - end <= combine_cut_off:
                end = y
            elif x - end > combine_cut_off:
                region_map[chr_name].append((start, end))
                start = x
                end = y
    if start is not None:
        region_map[chr_name].append((start, end))

    file_name = peak_ref_map[peak_ref_map.rfind("/")+1: peak_ref_map.rfind(".")]

    input_path = peak_ref_map[:peak_ref_map.rfind("/")+1]

    output_file = open(input_path+file_name+"_combined.csv", "r")

    for key, value in region_map.iteritems():
        output_file.write(">"+key+"\n")
        for region in value:
            output_file.write(region[0]+","+region[1]+"\n")
    output_file.close()
    return region_map

if __name__ == "__main__":
    # get_split_chr("chr3", 187450000, 187470000, cutoff=25)
    peak_combiner("/home/tmhbxx3/archive/25_refmap.csv", 4000)
