
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage
# from rpy2 import robjects
from cycler import cycler

def plotSaturation(title, variant, sample_names, data_values, types, verbose=False):
    width = variant.end-variant.start
    array = variant.members
    rep = variant.signals
    original_seed = variant.seed
    parameters = len(variant.members)

    fig = plt.figure(figsize=(8,2))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    title_font = {'fontname': 'Arial', 'size': '15', 'color': 'black', 'weight': 'bold',
                  'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
    axis_font = {'fontname': 'Arial', 'size': '15', 'weight': 'bold'}

    chromosome, start, end, cluster_number = title.split("_")
    ax.set_xlabel("bp", **axis_font)
    ax.set_ylabel("Signal", **axis_font)
    ax.set_title(chromosome + " " + start + "~"+end+" "+"variant "+str(int(cluster_number[-1])+1), **title_font)

    xvalues = np.arange(array.shape[1])*10

    # ax.set_prop_cycle(color=['grey'])
    # for i in range(parameters):
    #     ax.plot(xvalues, array[i, :], linewidth=0.1, alpha=0.5)

    ax.set_prop_cycle(color=['red'])
    # ax.plot(xvalues, rep, linewidth=0.5, label='representative', alpha=0.5)
    ax.fill_between(xvalues, 0, rep, alpha=0.5, color='red')

    peaks = variant.units

    for p in peaks:
        plt.axvspan(xmin=p.start-int(start)+width/100.0, xmax=p.end-int(start)-width/100.0, ymin=0, ymax=0.05,
                    facecolor='blue', alpha=0.5, edgecolor='black')

    #######

    # if original_seed != []:
    #     ax.set_prop_cycle(color=['blue'])
    #     ax.plot(xvalues, original_seed, linewidth=1, label='seed')

    # legend = ax.legend(loc='upper right', prop={'size': 2})

    # # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    # frame = legend.get_frame()
    # frame.set_facecolor('0.90')

    # Set the fontsize
    # for label in legend.get_texts():
    #     label.set_fontsize(5)
    #
    # for label in legend.get_lines():
    #     label.set_linewidth(0.75)  # the legend line width

    fig.set_tight_layout(True)

    variant_id = title.replace("_", ".")[3:]

    for type in types:
        fig.savefig("./pictures/"+type+'/'+variant_id+'.png', dpi=600, facecolor='w', edgecolor='w',
                    orientation='portrait', bbox_inches='tight')
    plt.close('all')


    if verbose:
        best_sample_index = variant.best_representative_sample()
        for i in range(len(best_sample_index)):
            index = best_sample_index[i]
            sample_name = sample_names[index]
            output_name = variant_id + " " + sample_name
            cur_data = data_values[index]
            plotBestSample(sample_name, output_name, cur_data, xvalues)
    return peaks

def plotBestSample(title, output_name, signals, xvalues):
    fig = plt.figure(figsize=(8, 2))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    title_font = {'fontname': 'Arial', 'size': '15', 'color': 'black', 'weight': 'bold',
                  'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
    axis_font = {'fontname': 'Arial', 'size': '15', 'weight': 'bold'}

    ax.set_xlabel("bp", **axis_font)
    ax.set_ylabel("Signal", **axis_font)
    ax.set_title(title, **title_font)

    ax.set_prop_cycle(color=['blue'])
    ax.fill_between(xvalues, 0, signals, alpha=1, color='blue')

    fig.set_tight_layout(True)

    fig.savefig("./pictures/" + output_name + '.png', dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', bbox_inches='tight')
    plt.close('all')

# def heatmap(path, name):
    # df = pd.DataFrame(array, index=index)
    # row_dist = pd.DataFrame(squareform(pdist(df, metric='euclidean')), index=index)
    # row_clusters = linkage(df.values, method='complete', metric='euclidean')
    # pd.DataFrame(row_clusters,
    #              columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
    #              index=['cluster %d' % (i + 1) for i in range(row_clusters.shape[0])])
    # fig = plt.figure(n)
    # fig.set_size_inches(10, 10, forward=True)
    #
    # ax = fig.add_subplot(111)
    #
    # cax = ax.matshow(df, interpolation='nearest', cmap='hot_r')
    # fig.colorbar(cax)
    #
    # ax.set_xticklabels([''] + list(df.columns))
    # ax.set_yticklabels([''] + list(df.index))
    #
    # plt.show()
    # robjects.r('''
    #     if (!require("gplots")) {
    #         install.packages("gplots", dependencies = TRUE)
    #         library(gplots)
    #     }
    #     if (!require("RColorBrewer")) {
    #         install.packages("RColorBrewer", dependencies = TRUE)
    #         library(RColorBrewer)
    #     }
    #
    #     data <- read.csv("'''+path+'''", comment.char="#")
    #     rnames <- data[,1]                            # assign labels in column 1 to "rnames"
    #     mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
    #     rownames(mat_data) <- rnames                  # assign row names
    #     my_palette <- colorRampPalette(c("blue", "yellow"))(n = 199)
    #     # col_breaks = c(seq(-1,0,length=100),  # for red
    #     # seq(0.01,1,length=100))           # for yellow
    #     # seq(0.51,1,length=100))             # for green
    #     png("./pictures/heatmaps_'''+name+'''.png",    # create PNG for the heat map
    #     width = 15*300,        # 5 x 300 pixels
    #     height = 15*300,
    #     res = 600,            # 300 pixels per inch
    #     pointsize = 8)        # smaller font size
    #     heatmap.2(mat_data,
    #     density.info="none",  # turns off density plot inside color legend
    #     trace="none",         # turns off trace lines inside the heat map
    #     margins =c(12,9),     # widens margins around plot
    #     col=my_palette,       # use on color palette defined earlier
    #     # breaks=col_breaks,    # enable color transition at specified limits
    #     dendrogram="none",     # only draw a row dendrogram
    #     Colv="NA")            # turn off column clustering
    #
    #     dev.off()               # close the PNG device
    #
    # # ''')

def plot_predict(data, representatives, allocs):
    fig = plt.figure(figsize=(8, 2))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel("bp")
    ax.set_ylabel("Signal")
    xvalues = np.arange(len(data)) * 10
    ax.set_prop_cycle(cycler('color', ['c', 'm', 'y', 'k']) +
                   cycler('lw', [1, 2, 3, 4]))
    for i in range(representatives.shape[0]):
        ax.plot(xvalues, representatives[i, :], linewidth=1, label='variant'+str(i+1)+" "+str(round(allocs[i],2)))

    ax.set_prop_cycle(color=['blue'])
    ax.plot(xvalues, data, linewidth=1, label='sample')

    legend = ax.legend(loc='best', prop={'size': 4})
    legend.draggable()
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize(10)

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width

    fig.savefig("./pictures/test.png", dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', bbox_inches='tight')

    plt.show()