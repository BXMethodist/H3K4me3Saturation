import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage
from region import region
# import rpy2.robjects as robjects


def plotSaturation(title, array, original_seed, rep, parameters):
    fig = plt.figure(figsize=(8,2))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_color_cycle(['grey'])
    for i in range(parameters):
        ax.plot(np.arange(array.shape[1]), array[i, :], linewidth=0.2)

    ax.set_color_cycle(['red'])
    ax.plot(np.arange(array.shape[1]), rep, linewidth=0.5)
    peaks = rep.copy()
    draw = np.max(rep)/10
    peaks[peaks < 200] = 0
    peaks[peaks > 200] = draw

    # ax.set_color_cycle(['grey'])
    ax.plot(np.arange(array.shape[1]), peaks, linewidth=1)
    plt.fill_between(np.arange(array.shape[1]), peaks, color='red')

    chromosome, start, end = title.split("_")[:-1]

    refregion = region(chromosome, start, end, signals=rep, peaks=peaks)

    if original_seed != []:
        ax.set_color_cycle(['blue'])
        ax.plot(np.arange(array.shape[1]), original_seed, linewidth=0.5)


    # if verbose:
    #     legend = ax.legend(parameters, loc='center right', bbox_to_anchor=(1.3, 0.5))

    fig.savefig("./pictures/"+title, dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', bbox_inches='tight')

    return refregion

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



