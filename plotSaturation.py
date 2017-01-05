# This is used for plot results from different peak cutoff, coverage, avgLength, number of peaks/regions

import matplotlib.pyplot as plt
import numpy as np

import os

from mpl_toolkits.axes_grid.inset_locator import inset_axes

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def plotSaturation(title, array, parameters, n, subplot=True, x_reverse=False, diagonal=False, plot_legend=False):

    fig = plt.figure(n)
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    cm = plt.get_cmap('gist_rainbow')
    ax.set_color_cycle([cm(1. * i / len(parameters)) for i in range(len(parameters))])

    ax.plot(array[:, 0], array[:, 1:], linewidth=2.0)
    # ax.plot(array, linewidth=2.0)

    # ax.set_title(title+" VS Sample Number", fontname="Times New Roman")
    # ax.set_xlabel("Number of H3K4me3 Chip-Seq Sample", fontname="Times New Roman")
    # ax.set_ylabel(title, fontname="Times New Roman")

    ax.set_title(title, fontname="Times New Roman")
    ax.set_xlabel("Cutoff", fontname="Times New Roman")
    ax.set_ylabel("Coverage", fontname="Times New Roman")

    if plot_legend:
        legend = ax.legend(parameters, loc='center right', bbox_to_anchor=(1.3, 0.5))

    if subplot:
        inside_ax = inset_axes(ax,
                                width="35%",  # width = 30% of parent_bbox
                                height="35%",  # height : 1 inch
                                loc=5)
        inside_ax.set_color_cycle([cm(1. * i / len(parameters)) for i in range(0, len(parameters))])

        if x_reverse:
            inside_ax.set_xlim(ax.get_xlim()[::-1])

        inside_ax.plot(array[20:, 0], array[20:, 1:], linewidth=2.0)

        inside_ax.tick_params(labelsize=10)

    if x_reverse:
        ax.set_xlim(ax.get_xlim()[::-1])

    if diagonal:
        ax.plot(array[:, 0], array[:, 0], linewidth=2.0, color="orange")

    if plot_legend:
        fig.savefig("./pictures/"+title, dpi=600, facecolor='w', edgecolor='w',
                    orientation='portrait', bbox_extra_artists=(legend,), bbox_inches='tight')
    else:
        fig.savefig("./pictures/" + title, dpi=600, facecolor='w', edgecolor='w',
                    orientation='portrait')


def csvFiles():
    files = [csv for csv in os.listdir(os.getcwd()) if csv.endswith(".csv")]
    return files


if __name__ == "__main__":
    parameters = [10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500]

    csvFiles = csvFiles()
    csvArrays = [(csv[:csv.find(".csv")], np.genfromtxt(csv, delimiter=',')) for csv in csvFiles]
    csvArrays = sorted(csvArrays, key=lambda x: int(x[0]))

    coverage = np.zeros((csvArrays[0][1].shape[0], len(parameters)+1))
    avgLength = np.zeros((csvArrays[0][1].shape[0], len(parameters)+1))
    regionNumber = np.zeros((csvArrays[0][1].shape[0], len(parameters)+1))

    coverage[:, 0] = csvArrays[0][1][:, 0]
    avgLength[:, 0] = csvArrays[0][1][:, 0]
    regionNumber[:, 0] = csvArrays[0][1][:, 0]

    Titles = ["Total Reference Map Size", "Average Peak Size", "Total Peak Number"]

    for i in range(len(parameters)):
        curArray = csvArrays[i][1]
        coverage[:, i+1] = curArray[:, 1]
        avgLength[:, i+1] = curArray[:, 2]
        regionNumber[:, i + 1] = curArray[:, 3]



    # plotSaturation(Titles[0], coverage, parameters, 1)
    plotSaturation(Titles[1], avgLength, parameters, 2, False)
    plotSaturation(Titles[2], regionNumber, parameters, 3)





