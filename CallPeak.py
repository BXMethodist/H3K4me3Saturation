"""
This module is used for call the peaks for a numpy array
"""

import pandas as pd, numpy as np

def callpeak(genome_signals, cutoff, step):
    """
    :param signals: a dictionary, key is chromosome name, value is a 1-d np array
    :param cutoff:  peak height cutoff
    :return: a list of peaks including ['chr', 'start', 'end', 'center', 'width_above_cutoff', 'total_signal',
            'height', 'height_logP']
    """

    results = []
    for key, value in genome_signals.iteritems():
        cur_peak_locations = np.where(value>=cutoff)[0]
        start = None
        end = None

        for location in cur_peak_locations:
            if start is None:
                start = location
                end = location
            elif location - end == 1:
                end = location
            elif location - end > 1:
                if start != end:
                    center = start + np.argmax(value[start:end])
                    peak_width = (end - start) * step
                    cur_total_signals = np.sum(value[start:end])
                    height = np.max(value[start:end])
                    height_logP = np.log10(height)
                    peak_start = start * step
                    peak_end = end * step
                    results.append([key, peak_start, peak_end, center,
                                    peak_width, cur_total_signals, height, height_logP])
                start = location
                end = location
        if start is not None:
            if start != end:
                center = start + np.argmax(value[start:end])
                peak_width = (end - start) * step
                cur_total_signals = np.sum(value[start:end])
                height = np.max(value[start:end])
                height_logP = np.log10(height)
                peak_start = start * step
                peak_end = end * step
                results.append([key, peak_start, peak_end, center,
                                peak_width, cur_total_signals, height, height_logP])
            start = None
            end = None
    return results