from DistinctAffinityPropagation import DistinctAffinityPropagation
from region import Region
import pandas as pd, os, numpy as np, pickle
from clusterUtils import get_map
from visualizationUtils import plotSaturation #, heatmap

def region_cluster(number_sample_used,
                   cutoff=75,
                   list_files=None,
                   directory="/home/tmhbxx3/archive/WigChrSplits/code/csv/",
                   affinity=np.corrcoef,
                   verbose=True):
    if not os.path.isdir("./pictures"):
        os.system("mkdir pictures")
    if not os.path.isdir("./tempcsv"):
        os.system("mkdir tempcsv")
    if not os.path.isdir("./cluster_csv"):
        os.system("mkdir cluster_csv")

    if not directory.endswith("/"):
        directory += "/"

    if list_files is None:
        list_files = [ x for x in os.listdir(directory) if x.endswith(".csv")]

    regions = []

    for file_name in list_files:
        print file_name
        cluster = DistinctAffinityPropagation(number_sample_used, affinity=affinity)
        df = pd.read_csv(directory+file_name, sep="\t")

        pos_surfix = file_name[file_name.rfind("/") + 1:-4]

        data = df.as_matrix()
        data_values = data[:, 1:].astype(np.float64)
        sample_names = data[:, 0]

        if len(sample_names) == 0:
            continue

        max_values = np.max(data_values, axis=1)
        first_ten_percentile = np.percentile(max_values, 2)
        if cutoff is not None:
            cutoff = min(cutoff, first_ten_percentile)
        else:
            cutoff = first_ten_percentile
        left_index = np.where(max_values > cutoff)[0]
        data_values = data_values[left_index, :]
        sample_names = sample_names[left_index]

        if data_values.shape[0] == 1:
            continue

        cluster.fit(data_values, 0.8, 0.3, cutoff=cutoff)

        data_values = cluster.data

        # print cluster.cluster_centers_indices_
        # print [len(x) for x in cluster.labels_]

        np.set_printoptions(threshold=np.nan)

        total = 0

        for label in cluster.labels_:
            total += len(label)
        result = np.zeros((total, data_values.shape[1]))
        result_index = np.zeros(total)

        cur_pos = 0

        for i in range(len(cluster.labels_)):
            result[cur_pos:cur_pos + len(cluster.labels_[i])] = data_values[cluster.labels_[i]]
            result_index[cur_pos:cur_pos + len(cluster.labels_[i])] = i
            cur_pos += len(cluster.labels_[i])
        df = pd.DataFrame(result, index=result_index)

        #save the cluster information after normalization
        df.to_csv("./cluster_csv/" + pos_surfix + "_clusters.csv", sep="\t")

        chromosome, start, end = pos_surfix.split("_")
        # print sample_names
        peak = Region(chromosome, int(start), int(end),
                      representatives=cluster.representation_ ,
                      seeds=data_values[cluster.cluster_centers_indices_],
                      variants_members=data_values,
                      labels=cluster.labels_,
                      sample_names=sample_names)

        regions.append(peak)

        if verbose:
            if len(peak.variants) > 1 and peak.plot:
                for i in range(len(peak.variants)):
                    variant = peak.variants[i]
                    plotSaturation(pos_surfix + "_cluster" + str(i), variant, peak.sample_names, peak.variants_members)
                pass
            else:
                for i in range(len(peak.variants)):
                    variant = peak.variants[i]
                    plotSaturation(pos_surfix + "_cluster" + str(i), variant)
                pass

        # if len(cluster.labels_) > 1:
        #     i = 0
        #     for label in cluster.labels_:
        #         plot_data = data_values[label]
        #         region = plotSaturation(pos_surfix + "_cluster" + str(i), plot_data,
        #                        data_values[cluster.cluster_centers_indices_[i]],
        #                        cluster.representation_[i], len(cluster.labels_[i]))
        #         i += 1
        #         regions += region
        #         # df = pd.DataFrame(plot_data)chr3_187450000_187470000.csv
        #         # df.to_csv("./tempcsv/" + "cluster" + str(i) + ".csv")
        #         # heatmap("./tempcsv/" + "cluster" + str(i) + ".csv", pos_surfix + "_cluster" + str(i))
        #
        # ### TO DO:
        # else:
        #     # print cluster.representation_
        #     # print "norm_X is ", cluster.data
        #     # region = plotSaturation(pos_surfix + "_cluster0", data_values, [],
        #     #                cluster.representation_[0], len(cluster.labels_[0]))
        #     # regions += region
        #
        #     # diable single cluster photo for now
        #     pass
    return regions


# if __name__ == "__main__":
    # open a reference map
# map_path ="./75_refmap_combined.csv"
# finished_job = os.listdir("/home/tmhbxx3/archive/WigChrSplits/code/csv/")
# files_read_for_clusters = get_map(map_path, finished_job=finished_job)
#
# import matplotlib.pyplot as plt
# print plt.gcf().canvas.get_supported_filetypes()


# regions = region_cluster(300, directory='./csv', verbose=False)



#
# data = regions[0].variants[0].members[3, :]
#
# from visualizationUtils import plot_predict
# from predict import optimize_allocs
#
# allocs = optimize_allocs(data, regions[0].representatives)
#
# plot_predict(data, regions[0].representatives, allocs)

# regions = region_cluster()
#
# print regions
#
# import pickle
#
# with open('75refmap_combined_3kb_regions' + '.pkl', 'wb') as f:
#     pickle.dump(regions, f, pickle.HIGHEST_PROTOCOL)
#
# f.close()
    # region_cluster(directory="./csv")
