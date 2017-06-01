# This module contains the function to create the clusters and draw the pictures

from DistinctAffinityPropagation import DistinctAffinityPropagation
from region import Region
import pandas as pd, os, numpy as np, pickle
from clusterUtils import get_map, smooth_normalization
from visualizationUtils import plotSaturation
from multiprocessing import Process, Queue

def region_cluster(number_sample_used,
                   cutoff=10,
                   list_files=None,
                   directory="/home/tmhbxx3/archive/WigChrSplits/code/csv/",
                   affinity=np.corrcoef,
                   verbose=True,
                   example=False,
                   process=20,
                   outputdir='./pictures/',
                   hide=False):
    if not os.path.isdir("./pictures"):
        os.system("mkdir pictures")
    if not os.path.isdir("./tempcsv"):
        os.system("mkdir tempcsv")
    if not os.path.isdir("./cluster_csv"):
        os.system("mkdir cluster_csv")

    if not directory.endswith("/"):
        directory += "/"

    if list_files is None:
        list_files = [x for x in os.listdir(directory) if x.endswith(".csv")]

    list_files = ['chr3_187455150_187464690.csv']

    chunks = []
    cur_index = 0
    reminder = len(list_files) % process
    chunk_size = len(list_files) / process
    for i in range(process):
        if reminder > 0:
            chunks.append(list_files[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
            cur_index += 1
            reminder -= 1
        else:
            chunks.append(list_files[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

    total_chunk_size = 0
    for chunk in chunks:
        total_chunk_size += len(chunk)
    if total_chunk_size != len(list_files):
        print 'multiple processes chunk size is not correct'
        return None

    queue = Queue()
    processes = []

    for i in range(process):
        cur_chunk = chunks[i]
        p = Process(target=Region_Cluster_Process, args=(queue, number_sample_used, cutoff, cur_chunk, directory, np.corrcoef,
                   verbose, example, outputdir, hide))
        processes.append(p)
        p.start()

    regions = []

    for i in range(process):
        cur_regions = queue.get()
        regions += cur_regions

    for p in processes:
        p.join()

    return regions



def Region_Cluster_Process(queue,
                           number_sample_used,
                           cutoff,
                           list_files,
                           directory,
                           affinity,
                           verbose,
                           example,
                           outputdir,
                           hide):
    regions = []
    for file_name in list_files:
        cluster = DistinctAffinityPropagation(number_sample_used, affinity=affinity)
        df = pd.read_csv(directory+file_name, sep="\t")

        pos_surfix = file_name[file_name.rfind("/") + 1:-4]

        data = df.as_matrix()
        data_values = data[:, 1:].astype(np.float64)
        sample_names = data[:, 0]

        if len(sample_names) == 0:
            continue

        max_values = np.max(data_values, axis=1)
        # first_ten_percentile = np.percentile(max_values, 2)
        # if cutoff is not None:
        #     cutoff = min(cutoff, first_ten_percentile)
        # else:
        #     cutoff = first_ten_percentile
        left_index = np.where(max_values > cutoff)[0]
        data_values = data_values[left_index, :]
        sample_names = sample_names[left_index]

        if data_values.shape[0] == 1:
            continue

        cluster.fit(data_values, 0.8, 0.3, cutoff=cutoff)

        data_values = cluster.data

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
                    types = set()
                    variant = peak.variants[i]
                    for key in peak.transitions.keys():
                        if i in key:
                            for type in peak.transitions[key]:
                                types.add(type)
                    types = list(types)
                    plotSaturation(pos_surfix + "_cluster" + str(i), variant, peak.sample_names, peak.variants_members,
                                   types=types,
                                   verbose=example,
                                   outputdir=outputdir)
                    if hide:
                        plotSaturation(pos_surfix + "_cluster" + str(i), variant, peak.sample_names,
                                       peak.variants_members,
                                       types=types,
                                       verbose=example,
                                       outputdir='./web_pictures/',
                                       hide=hide)
            #     pass
            # else:
            #     for i in range(len(peak.variants)):
            #         variant = peak.variants[i]
            #         plotSaturation(pos_surfix + "_cluster" + str(i), variant, peak.sample_names, peak.variants_members)
            #     pass

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
    queue.put(regions)
    return


# if __name__ == "__main__":
    # open a reference map
# map_path ="./75_refmap_combined.csv"
# finished_job = os.listdir("/home/tmhbxx3/archive/WigChrSplits/code/csv/")
# files_read_for_clusters = get_map(map_path, finished_job=finished_job)
#
# import matplotlib.pyplot as plt
# print plt.gcf().canvas.get_supported_filetypes()


regions = region_cluster(300, directory='./csv', verbose=False, example=False, hide=False)




#### This part is to draw the prediction demo
# from visualizationUtils import plot_predict
from predict import optimize_allocs
# this part is for draw the demo for the example of scipy minimizer
df = pd.read_csv('./csv/chr3_187455150_187464690.csv', sep='\t',index_col=0)
from refMapUtils import load_obj
def get_wig_signals(wig, regions):
    cur_wig = load_obj(wig[:-4])
    results = []
    for region in regions:
        if region[0] in cur_wig.genome:
            cur_signals = cur_wig.genome[region[0]].get_signals(region[1], region[2])
            # print region
            # print cur_signals
            # print type(region[0]), type(region[1])
            return cur_signals

# wigs = ['./wigs/'+x for x in os.listdir('./wigs') if x.find('BI.CD4+')!=-1]

# signals = []
# index = []

# for w in wigs:
#     signals.append([x for x in get_wig_signals(w, [['chr3', 187455150, 187464690]])])
#     index.append(w[7:])
#
# df = pd.DataFrame(signals, index=index)

errors = []
correlations = []
results = []
for i in range(df.shape[0]):
    data = df.ix[i, :]

    data = smooth_normalization(data)
    # print data
# print data

    allocs = optimize_allocs(data.T[0], regions[0].representatives)
    # print allocs
    predicted = 0
    for j in range(regions[0].representatives.shape[0]):
        signal = (np.sum(data)/np.sum(regions[0].representatives[j])) * regions[0].representatives[j]
        # print np.sum(signal), np.sum(data)
        predicted += signal * allocs[j]
    # predicted = predicted * np.sum(data)/np.sum(predicted)

    # print np.sum(predicted), np.sum(data)
    # print data.T[0]
    # break
    # print data[:, 0].shape, predicted.shape
    # print data.T[0]-predicted
    # print (data.T[0]-predicted).astype(int)
    # print np.absolute(data.T[0]-predicted).astype(int)
    # print np.sum(predicted)
    # print data.T[0].sum()
    error = np.absolute(data.T[0]-predicted).sum()/np.sum(predicted)
    corr = np.corrcoef(predicted, data.T[0])[0,1]
    errors.append(error)
    correlations.append(corr)
    results.append((df.index[i], error, corr, allocs[0], allocs[1], 1-np.sum(allocs), abs(1-np.sum(allocs))))
    # break
    # print error, 'error rate'
    # print corr, 'correlation'
    # if corr < 0.7:
    #     print i, corr, error, error/np.sum(data.T[0])
    #     plot_predict(data, regions[0].representatives, allocs)
    #     break
    # print allocs, i, np.sum(predicted), np.sum(data), error
    # print allocs[1] > 0.8, allocs[1]
    # if allocs[1] > 0.8:
    #     print allocs, i, error
    # print error
    # if error < 0.4:
    #     print i
    #     print allocs
    # if i == 37:
    #     plot_predict(data, regions[0].representatives, allocs)
    #     break


df = pd.DataFrame(results)
df.columns = ['sample_name', 'error', 'correlation', 'alloc1', 'alloc2', 'error towards 1', 'abs error towards 1']

print np.corrcoef(df['correlation'], df['abs error towards 1'])
df.to_csv('error_correlation_prediction.csv', index=None)
print np.median(errors)
print np.median(correlations)
#####

# regions = region_cluster()
#
# print regions
#

# import pickle
#
# with open('./pkl/100refmap_100_10_2200' + '.pkl', 'wb') as f:
#     pickle.dump(regions, f, pickle.HIGHEST_PROTOCOL)
#
# f.close()
    # region_cluster(directory="./csv")

