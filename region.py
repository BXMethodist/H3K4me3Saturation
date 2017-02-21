
import numpy as np

class region():
    """
    Store the inforamtion of a reference peak:
        Parameters
        ----------
        chromosome: chromosome name
        start: start position
        end: end position
        step: step or bin size use to store the signal

        Attributes
        ----------
        chromosome name
        start
        end
        signals: np.array that keep the peaks signal
        step
        """

    def __init__(self, chromosome, start, end, signals, step=10, cutoff=200):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.cutoff = cutoff

        peaks_by_cutoff= callpeakbycutoff(self.chromosome, self.start, self.end, self.signals, cutoff=self.cutoff)

        self.peaks = []

        for peak in peaks_by_cutoff:
            cur_chr, cur_start, cur_end, cur_signals = peak
            cur_peaks = callpeak(cur_chr, cur_start, cur_end, cur_signals)
            self.peaks += cur_peaks


class peak():
    def __init__(self, chromosome, start, end, signals, step=10):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.width = end - start
        self.height = np.max(self.signals)


def callpeak(chromosome, start, end, signals):
    peaks = []

    frontiers = []

    frontiers.append((chromosome, start, end, signals))

    while len(frontiers) > 0:
        cur_chromosome, cur_start, cur_end, cur_signals = frontiers.pop(0)
        split_index = splitable(cur_chromosome, cur_start, cur_end, cur_signals)
        if split_index is None:
            peaks.append((cur_chromosome, cur_start, cur_end, cur_signals))
        else:
            left = (cur_chromosome, cur_start, cur_start+split_index, cur_signals[:split_index])
            right = (cur_chromosome, cur_start+split_index, cur_end, cur_signals[split_index:])
            frontiers.append(left)
            frontiers.append(right)

    # print peaks
    peaks_obj = []
    for p in peaks:
        chromosome, start, end, signals = p
        p_obj = peak(chromosome, start, end, signals)
        peaks_obj.append(p_obj)
    return peaks_obj

def splitable(chromosome, start, end, signals):
    local_min = np.r_[True, signals[1:] < signals[:-1]] & np.r_[signals[:-1] < signals[1:], True]
    local_max = np.r_[True, signals[1:] > signals[:-1]] & np.r_[signals[:-1] > signals[1:], True]

    min_index = np.where(local_min)[0]
    max_index = np.where(local_max)[0]

    left = max_index[0]
    right = max_index[-1]

    split_region = signals[left:right]
    if len(split_region) == 0:
        return None

    min = np.where(split_region==np.min(split_region))[0][0]

    left_part = signals[:left+min]
    right_part = signals[right-min:]

    left_peak_height = np.max(left_part)
    right_peak_height = np.max(right_part)
    min_height = split_region[min]

    if left_peak_height > 2 * min_height and right_peak_height > 2 * min_height:
        return left+min
    else:
        return None

def callpeakbycutoff(chromosome, start, end, signals, cutoff=200, step=10):
    peaks = []

    peaks_index = np.where(signals >= cutoff)[0]
    if len(peaks_index) == 0:
        return peaks
    elif len(peaks_index) == 1:
        index = peaks_index[0]
        peak = (chromosome, start+index, start+index+1, np.asarray([signals[index]]))
        peaks.append(peak)
        return peaks

    cur_start = peaks_index[0]
    prev = peaks_index[0]
    for i in range(1, len(peaks_index)):
        cur_index = peaks_index[i]
        if cur_index - prev == 1 and i != len(peaks_index):
            prev = cur_index
        elif cur_index - prev == 1 and i == len(peaks_index):
            peak = (chromosome, start+cur_start, start+cur_index, signals[cur_start:cur_index])
            peaks.append(peak)
        else:
            peak = (chromosome, start+cur_start, start+prev, signals[cur_start:prev])
            peaks.append(peak)
            prev = cur_index
            cur_start = cur_index
    return peaks


# a,b,c,d =('chr3', 66443698, 66443712, np.asarray([ 202.01670923,  204.06618201,  205.70825071,  206.70408383,
#         207.14969609,  207.46192666,  207.50757453,  207.31251443,
#         206.85803726,  206.26042879,  205.7258182 ,  205.06906288,
#         203.8992773 ,  202.21577459]))
#
# print callpeak(a,b,c,d)[0].signals