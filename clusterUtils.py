import numpy as np, os, pandas as pd


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


def get_split_chr(chr_name, start, end, vector_size=10000, step=10, cutoff=100,
                  delimiter="\t", path="/home/tmhbxx3/archive/WigChrSplits"):
    if not path.endswith("/"):
        path += "/"

    if abs(start - end) >= vector_size*step*2:
        return

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
    result_array[result_array<=cutoff] = 0
    row_to_delete = []
    for i in range(result_array.shape[0]):
        if np.sum(result_array[i, :]) == 0:
            row_to_delete.append(i)
    result_array = np.delete(result_array, row_to_delete, 0)

    sample_names = [file_names[i][:file_names[i].find("_")] for i in range(len(file_names)) if i not in row_to_delete]
    df = pd.DataFrame(result_array, index=sample_names)

    output_elements = [chr_name, str(start), str(end)]
    output = "_".join(output_elements)

    df.to_csv(output + ".csv", sep=delimiter)

    return result_array

if __name__ == "__main__":
    get_split_chr("chr3", 187450000, 187470000)