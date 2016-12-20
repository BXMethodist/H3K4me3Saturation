import numpy as np, os, pandas as pd


def find_optional_pos_refmap(path="./200_refmap.csv"):
    refmap_file = open("./200_refmap.csv", "r")

    chrX = False
    # largest = float("-inf")
    results = []
    position = None
    for line in refmap_file.readlines():
        if line.startswith(">chrX"):
            chrX = True
        elif line.startswith(">"):
            chrX = False
        if chrX and (not line.startswith(">")):
            line = line.rstrip().split(",")
            start, end = tuple(line)
            result = abs(int(start) - int(end))
            line.append(result)
            results.append(line)

    results = sorted(results, key=lambda x: x[2], reverse=True)
    refmap_file.close()
    return results

def genome_region(chr_name, start, end, step=10, wigPath="/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"):
    # start and end are actual genome position
    output_elements = [chr_name, str(start), str(end)]
    output = "_".join(output_elements)

    if not wigPath.endswith("/"):
        wigPath += "/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]

    start /= step
    end /= step if end%step == 0 else end/step+1

    colsize = end - start
    rowsize = len(wigFiles)
    data = np.zeros((rowsize, colsize))

    sample_names = [wigFiles[i][:-4] for i in range(len(wigFiles))]

    for i in range(len(wigFiles)):
        data[i, :] = get_position_array(wigPath+wigFiles[i], chr_name, start, end)

    df = pd.DataFrame(data, index=sample_names)

    df.to_csv(output+".csv")
    return df

def get_position_array(path, chr_name, start, end):
    wigFile = open(path, "r")

    cur_position = 0

    data_array = np.zeros(end-start)

    target = False

    for line in wigFile.readlines():
        cur_line = line.rstrip().split()
        if len(cur_line) > 1:
            if cur_line[1].startswith("chrom="):
                cur_chr_name = cur_line[1][cur_line[1].find("=") + 1:]
            if cur_line[2].startswith("start="):
                cur_start = int(cur_line[2][cur_line[2].find("=") + 1:])
            if cur_chr_name == chr_name:
                target = True
                cur_position = cur_start
            else:
                if target == True:
                    break
        else:
            if target:
                if start <= cur_position < end:
                    data_array[cur_position-start] = float(cur_line[0])
                cur_position += 1
                if cur_position >= end:
                    break

    wigFile.close()

    return data_array


if __name__ == "__main__":
    genome_region("chrX", 17393000, 17397000)
