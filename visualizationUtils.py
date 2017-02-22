import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage
from region import region, callpeakbycutoff
# import rpy2.robjects as robjects


def plotSaturation(title, array, original_seed, rep, parameters):
    fig = plt.figure(figsize=(8,2))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    chromosome, start, end, cluster_number = title.split("_")
    ax.set_xlabel("bp")
    ax.set_ylabel("Signal")
    ax.set_title(chromosome + " " + start + "~"+end+" "+"variant "+str(int(cluster_number[-1])+1))
    xvalues = np.arange(array.shape[1])*10

    ax.set_color_cycle(['grey'])
    for i in range(parameters):
        ax.plot(xvalues, array[i, :], linewidth=0.2, alpha=0.5)

    ax.set_color_cycle(['red'])
    ax.plot(xvalues, rep, linewidth=0.5)

    peaks = region(chromosome, int(start), int(end), rep).peaks

    for p in peaks:
        plt.axvspan(xmin=p.start-int(start), xmax=p.end-int(start), ymin=0, ymax=0.05, facecolor='red', alpha=0.5, edgecolor='red')

    #######

    if original_seed != []:
        ax.set_color_cycle(['blue'])
        ax.plot(xvalues, original_seed, linewidth=0.5)


    # if verbose:
    #     legend = ax.legend(parameters, loc='center right', bbox_to_anchor=(1.3, 0.5))

    fig.savefig("./pictures/"+title, dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', bbox_inches='tight')

    plt.close('all')
    return peaks

# def heatmap(path, name):
#     # df = pd.DataFrame(array, index=index)
#     # row_dist = pd.DataFrame(squareform(pdist(df, metric='euclidean')), index=index)
#     # row_clusters = linkage(df.values, method='complete', metric='euclidean')
#     # pd.DataFrame(row_clusters,
#     #              columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
#     #              index=['cluster %d' % (i + 1) for i in range(row_clusters.shape[0])])
#     # fig = plt.figure(n)
#     # fig.set_size_inches(10, 10, forward=True)
#     #
#     # ax = fig.add_subplot(111)
#     #
#     # cax = ax.matshow(df, interpolation='nearest', cmap='hot_r')
#     # fig.colorbar(cax)
#     #
#     # ax.set_xticklabels([''] + list(df.columns))
#     # ax.set_yticklabels([''] + list(df.index))
#     #
#     # plt.show()
#     robjects.r('''
#         if (!require("gplots")) {
#             install.packages("gplots", dependencies = TRUE)
#             library(gplots)
#         }
#         if (!require("RColorBrewer")) {
#             install.packages("RColorBrewer", dependencies = TRUE)
#             library(RColorBrewer)
#         }
#
#         data <- read.csv("'''+path+'''", comment.char="#")
#         rnames <- data[,1]                            # assign labels in column 1 to "rnames"
#         mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
#         rownames(mat_data) <- rnames                  # assign row names
#         my_palette <- colorRampPalette(c("blue", "yellow"))(n = 199)
#         # col_breaks = c(seq(-1,0,length=100),  # for red
#         # seq(0.01,1,length=100))           # for yellow
#         # seq(0.51,1,length=100))             # for green
#         png("./pictures/heatmaps_'''+name+'''.png",    # create PNG for the heat map
#         width = 15*300,        # 5 x 300 pixels
#         height = 15*300,
#         res = 600,            # 300 pixels per inch
#         pointsize = 8)        # smaller font size
#         heatmap.2(mat_data,
#         density.info="none",  # turns off density plot inside color legend
#         trace="none",         # turns off trace lines inside the heat map
#         margins =c(12,9),     # widens margins around plot
#         col=my_palette,       # use on color palette defined earlier
#         # breaks=col_breaks,    # enable color transition at specified limits
#         dendrogram="none",     # only draw a row dendrogram
#         Colv="NA")            # turn off column clustering
#
#         dev.off()               # close the PNG device
#
#     ''')



