import numpy as np, os, pandas as pd
from collections import defaultdict


def genome_size_chrom(path="/home/tmhbxx3/archive/ref_data/hg19/hg19_chr_sizes.txt"):
    """
    :param path: the location where is the genome size file
    :return: a dictionary in which key is the chromosome name (str), value is the chromosome size (int)
    """
    genome = {}
    genome_size_file = open(path, "r")
    for line in genome_size_file.readlines():
        chr_name, chr_size = line.split()
        genome[chr_name] = int(chr_size)
    genome_size_file.close()
    return genome


def get_split_chr(chr_name, start, end, vector_size=10000, step=10, merge=10, cutoff=25,
                  delimiter="\t", path="/home/tmhbxx3/archive/WigChrSplits"):
    """
    :param chr_name: name of the chromosome
    :param start: start position, number before step/bin
    :param end: end position, number before step/bin
    :param vector_size: default is 10000, this is decided by how you store the splitted wig files
    :param step: bin size, default is 10.
    :param merge: the number of bin want to merge, currently is not implemented
    :param cutoff: currently is 25. The signal cutoff.
    :param delimiter: the delimiter when splitted wig files used
    :param path: where the splitted wig files are stored
    :return: the output filename (without directory)
    """
    print chr_name, start, end
    if not path.endswith("/"):
        path += "/"

    if abs(start - end) >= vector_size*step*2:
        return

    if not os.path.isdir(path+"code/csv/"):
        os.system("mkdir "+path+"code/csv/")

    output_path = path + chr_name + "_"

    unit_size = vector_size * step

    start_name = start - start % unit_size
    end_name = (end/unit_size)*unit_size + unit_size if end % unit_size != 0 else end

    folder_names = []

    while start_name < end_name:
        folder_name = output_path + str(start_name) + "_" + str(start_name + unit_size) + "_" + str(step) + "/"
        if not os.path.exists(folder_name):
            allfolders = os.listdir(path)
            for folder in allfolders:
                if folder.find(chr_name + "_" + str(start_name) + "_") != -1:
                    folder_name = path + folder + '/'
        start_name += unit_size
        folder_names.append(folder_name)
        print folder_name

    sample_number = len(os.listdir(folder_names[0]))

    sample_size = (end - start)/step if (end - start) % step == 0 else (end - start)/step + 1

    result_array = np.zeros((sample_number, sample_size))

    # print result_array.shape

    start_index = start % unit_size / step
    end_index = end % unit_size/ step
    output_index = 0

    for i in range(len(folder_names)):
        if i == 0:
            cur_start = start_index
        else:
            cur_start = 0
        if i == len(folder_names) - 1:
            cur_end = end_index
        else:
            cur_end = vector_size

        # print cur_start, cur_end, output_index

        file_names = os.listdir(folder_names[i])
        for j in range(len(file_names)):
            cur_array = np.genfromtxt(folder_names[i]+file_names[j])

            # print cur_array.shape

            result_array[j, output_index:output_index+cur_end-cur_start] = cur_array[cur_start:cur_end]
        output_index += cur_end-cur_start

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

    # if sample_number - len(row_to_delete) > 5:

    df = pd.DataFrame(result_array, index=sample_names)

    output_elements = [chr_name, str(start), str(end)]
    output = "_".join(output_elements)

    df.to_csv(path+"code/csv/" + output + ".csv", sep=delimiter)

    return output+".csv"
    # else:
    #     return None


def quantile_normalization(array, axis=0):
    """
    :param array: a numpy array, row is the sample, column is the index in specific reference peak
    :param axis: Normalization axis
    :return: numpy array
    """
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
    """
    :param array:  numpy array
    :param axis: normalize axis
    :param target_signal: target mean, normalize all samples' mean to the target
    :return: numpy array
    """
    if array.ndim == 1:
        print "please convert to samples as array, shape (n_samples, ) which is 2d array"
        return

    n_sample = array.shape[0]
    factors = target_signal/(np.mean(array, axis=axis)).reshape(n_sample, 1)

    return np.multiply(array, factors)


def smooth_normalization(array, smooth=20, axis=1):
    """
    :param array: numpy array
    :param smooth: smooth target, note this is value after step which means 20 is actually 200 bp.
    :param axis:
    :return: numpy array
    """
    df = pd.DataFrame(array)
    return df.rolling(min_periods=1, window=smooth, center=True, axis=axis).mean().values


def max_normalization(array, max_peak=500, axis=1):
    """
    :param array: numpy array
    :param max_peak: normalize the maximum summit in all the sample to the same max_peak target.
    :param axis:
    :return: numpy array
    """
    if array.ndim == 1:
        print "please convert to samples as array, shape (n_samples, ) which is 2d array"
        return

    n_sample = array.shape[0]
    factors = max_peak / (np.amax(array, axis=axis)).reshape(n_sample, 1)
    return np.multiply(array, factors)


def peak_combiner(peak_ref_map, combine_cut_off, step=10):
    """
    TO DO: the logic for recombine the peak need to be redefined, currently, the region are combined too much

    :param peak_ref_map: a reference map containing chromosome name, start, end
    :param combine_cut_off: the cut_off of combining peak
    :param step: bin/step size
    :return: a combined reference map
    """
    n= 0
    combine_cut_off = combine_cut_off/step
    input_file = open(peak_ref_map, "r")
    region_map = defaultdict(list)

    start = None
    end = None
    chr_name = None
    max_peak = None
    for line in input_file.readlines():
        if line.startswith(">"):
            if start is not None:
                region_map[chr_name].append((start, end))
                n +=1
            chr_name = line.rstrip()[1:]
            start = None
            end = None
            max_peak = None
        else:
            line = line.rstrip().split(",")
            x, y = int(line[0]), int(line[1])
            cur_peak_width = y - x

            if start is None:
                start = x
                end = y
                max_peak = y - x
            elif x - end <= combine_cut_off: # and x - end <= max_peak:
                end = y
                if cur_peak_width > max_peak:
                    max_peak = cur_peak_width
            elif x - end > combine_cut_off or x - end > max_peak:
                region_map[chr_name].append((start, end))
                n += 1
                start = x
                end = y
                max_peak = y - x
    print n
    input_file.close()

    if start is not None:
        region_map[chr_name].append((start, end))

    file_name = peak_ref_map[peak_ref_map.rfind("/")+1: peak_ref_map.rfind(".")]

    input_path = peak_ref_map[:peak_ref_map.rfind("/")+1]

    output_file = open(input_path+file_name+"_combined.csv", "w")

    for key, value in region_map.iteritems():
        output_file.write(">"+key+"\n")
        for region in value:
            output_file.write(str(region[0])+","+str(region[1])+"\n")
    output_file.close()
    return region_map


def get_map(refmap, step=10, sep=",", finished_job=()):
    """
    Currently restrict to chr3
    :param refmap: the txt file containing the reference map, with chromosome name, start, end (value after step)
    :param step: step/bin size
    :param sep:
    :param finished_job: finished job names
    :return: the file names of the reference peak cluster table file name
    """
    finished_job = set(finished_job)

    input_file = open(refmap, "r")

    files_read_for_cluster = []
    chr_name = None

    region_map = defaultdict(list)

    for line in input_file.readlines():
        if line.startswith(">"):
            chr_name = line.rstrip()[1:]
        else:
            line = line.rstrip().split(sep)
            # print line
            start = int(line[0]) * step
            end = int(line[1]) * step
            region_map[chr_name].append((start, end))

    for chr_name in region_map.keys():
        for line in region_map[chr_name]:
            start, end = line

            output_name = chr_name + "_" + str(start) + "_" + str(end) + ".csv"

            if output_name not in finished_job:
                file_name = get_split_chr(chr_name, start, end)
                if file_name is not None:
                    files_read_for_cluster.append(file_name)
    return files_read_for_cluster


def remove_duplicate(pairs):
    """Remove duplicate pairs and return distinct dissimilar pairs

        Parameters
        ----------
        pairs: ndarray, contains pair of sample index.

        Return
        ----------
        a list of tuple containing pair of sample index.

        Notes
        -----
        This function should only be called by .fit function

        References
        ----------
        """

    distinct_set = set()

    for i in range(pairs.shape[0]):
        x, y = pairs[i]
        if (x, y) not in distinct_set and (y, x) not in distinct_set and x != y:
            distinct_set.add((x, y))
    distinct_set = list(distinct_set)
    return distinct_set


if __name__ == "__main__":
    # get_split_chr("chr3", 122562770, 122562820, cutoff=25)

    # prefix = "/home/tmhbxx3/archive/refmap_saturation/code/"
    # surffix = "_refmap.csv"
    # cutoffs = [75]
    # for cutoff in cutoffs:
    #     path = str(cutoff) + surffix
    #     peak_combiner(path, 3000)

    # get_split_chr('chr3', 3167850, 3169860)
    # widths = []
    f = open("75_refmap_combined.csv", 'r')
    info = [line.strip().split(",") for line in f.readlines() if not line.startswith(">")]
    widths = [(int(y) - int(x))*10 for x, y in info]

    print np.percentile(widths, 50)
    # import matplotlib.pyplot as plt
    # plt.xlim((0, 10000))
    # plt.ylim((0, 200))
    # plt.hist(widths, bins='auto')
    # plt.show()

