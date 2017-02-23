
import numpy as np
import math
from scipy.stats import pearsonr


class region():
    """
    Parameters
        ----------
        chromosome: chromosome name
        start: start position
        end: end position
        step: step or bin size use to store the signal
        variants: the representative of each clusters
        variants_members: all samples signals, ndarray, each row is a sample
        seeds: seed signals used for each cluster
        labels: the index (row number) of samples in each clusters.
    """
    def __init__(self, chromosome, start, end, representatives, variants_members, seeds, labels, step=10):
        self.chromosome = chromosome,
        self.start = start
        self.end = end
        self.representatives = representatives
        self.variants_members = variants_members
        self.seeds = seeds
        self.labels = labels
        self.step = step
        self.variants = []
        self.units = set()

        self.create_variants()
        self.merge_units_across_variants()
        self.merge_variants()

    def create_variants(self):
        for i in range(len(self.representatives)):
            cur_variant = variant(self.chromosome, self.start, self.end, self.representatives[i],
                                  self.variants_members[self.labels[i]], self.seeds[i])
            self.variants.append(cur_variant)

            for unit in cur_variant.units:
                self.units.add((unit.start, unit.end, i))

        return


    def merge_units_across_variants(self):
        # initiate variants object
        # call units in variants object
        # merge the units between variants
        # merge variants by two creterias 1. have the same locations of units, 2. correlation is bigger than xxx

        group_units_by_overlap = []
        while len(self.units) > 0:
            unit = self.units.pop()
            cur_overlap = self.get_overlap(unit)
            group_units_by_overlap.append(cur_overlap)

        group_need_to_merge = []

        for group in group_units_by_overlap:
            print group

            related_variants = set([x[2] for x in group])
            if len(related_variants) == 1:
                continue
            left_border = (min([x[0] for x in group])-self.start)/self.step
            right_border = (max([x[1] for x in group])-self.start)/self.step

            group_matrix = []
            for v in related_variants:
                group_matrix.append(self.representatives[v][left_border:right_border])
            group_matrix = np.asarray(group_matrix)

            # print pearsonr(group_matrix[0], group_matrix[1])


            distances_matrix = np.corrcoef(group_matrix)

            # print distances_matrix
            close_pairs = np.asarray(np.where(distances_matrix >= 0)).T
            print close_pairs





    def merge_variants(self):
        pass

    def get_overlap(self, unit):
        overlap_units = set()
        overlap_units.add(unit)
        target_units = [unit]

        while len(target_units) > 0 :
            cur_unit = target_units.pop(0)
            need_to_remove = []
            for u in self.units:
                if cur_unit[0] <= u[0] <= cur_unit[1] or \
                    cur_unit[0] <= u[1] <= cur_unit[1]:
                    overlap_units.add(u)
                    target_units.append(u)
                    need_to_remove.append(u)
            for u in need_to_remove:
                self.units.remove(u)
        return overlap_units



class variant():
    """
    Store the inforamtion of a reference variant:
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

    def __init__(self, chromosome, start, end, signals, members, seed, step=10):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.cutoff = np.max(signals)*0.1
        self.convex_cutoff = np.max(signals)/3
        self.members = members
        self.seed = seed

        units_by_cutoff= callunitbycutoff(self.chromosome, self.start, self.end, self.signals,
                                          cutoff=self.cutoff, step=self.step)
        # print [(peak[0], peak[1], peak[2]) for peak in peaks_by_cutoff]
        # print peaks_by_cutoff
        merged_units = merge_unit(self.chromosome, self.start, self.end, len(self.signals)*self.step/40.0,
                                  units_by_cutoff, self.signals, self.step)

        self.units = []

        for unit in merged_units:
            cur_chr, cur_start, cur_end, cur_signals = unit
            cur_units = callunit(cur_chr, cur_start, cur_end, cur_signals, self.step, self.convex_cutoff)
            self.units += cur_units

        # self.peaks = peaks_by_cutoff


class unit():
    def __init__(self, chromosome, start, end, signals, step=10):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.width = end - start
        if len(self.signals) != 0:
            self.height = np.max(self.signals)
        else:
            self.height = None
        # print self.chromosome, self.start, self.end, self.height


def callunit(chromosome, start, end, signals, step, convex_cutoff):
    units = []

    frontiers = []

    frontiers.append((chromosome, start, end, signals))

    while len(frontiers) > 0:
        cur_chromosome, cur_start, cur_end, cur_signals = frontiers.pop(0)
        split_index = splitable(cur_chromosome, cur_start, cur_end, cur_signals, convex_cutoff)
        if split_index is None:
            units.append((cur_chromosome, cur_start, cur_end, cur_signals))
        else:
            left = (cur_chromosome, cur_start, cur_start+split_index*step, cur_signals[:split_index])
            right = (cur_chromosome, cur_start+split_index*step, cur_end, cur_signals[split_index:])
            frontiers.append(left)
            frontiers.append(right)

    units = sorted(units, key=lambda x: x[1])
    units_obj = []
    for p in units:
        chromosome, start, end, signals = p
        u_obj = unit(chromosome, start, end, signals)
        units_obj.append(u_obj)
    return units_obj

def splitable(chromosome, start, end, signals, convex_cutoff):
    local_min = np.r_[True, signals[1:] < signals[:-1]] & np.r_[signals[:-1] < signals[1:], True]
    local_max = np.r_[True, signals[1:] > signals[:-1]] & np.r_[signals[:-1] > signals[1:], True]

    min_index = np.where(local_min)[0]
    max_index = np.where(local_max)[0]

    if len(max_index)==0 or len(min_index) == 0:
        return None

    left = max_index[0]
    right = max_index[-1]

    split_region = signals[left:right]
    if len(split_region) == 0:
        return None

    split_local_mins = np.r_[True, split_region[1:] < split_region[:-1]] & np.r_[split_region[:-1] < split_region[1:], True]

    split_local_mins_index = np.where(split_local_mins)[0]

    splitable_indexes = []

    for min in split_local_mins_index:
        left_part = signals[:left+min]
        right_part = signals[left+min:]

        left_unit_height = np.max(left_part)
        right_unit_height = np.max(right_part)
        min_height = split_region[min]

        # print left_peak_height, right_peak_height, min_height
        if left_unit_height > 2 * min_height and right_unit_height > 2 * min_height and \
                (left_unit_height >= convex_cutoff and right_unit_height >= convex_cutoff):
            print "split"
            splitable_indexes.append(left+min)
    if len(splitable_indexes) == 0:
        return None
    else:
        return int(np.mean(splitable_indexes))

def callunitbycutoff(chromosome, start, end, signals, cutoff, step=10):
    units = []

    units_index = np.where(signals >= cutoff)[0]
    if len(units_index) == 0:
        return units
    elif len(units_index) == 1:
        index = units_index[0]
        peak = (chromosome, start+index*step, start+(index+1)*step, np.asarray([signals[index]]))
        units.append(peak)
        return units

    cur_start = units_index[0]
    prev = units_index[0]
    for i in range(1, len(units_index)):
        cur_index = units_index[i]
        if cur_index - prev == 1 and i != len(units_index)-1:
            prev = cur_index
        elif cur_index - prev == 1 and i == len(units_index)-1:
            peak = (chromosome, start+cur_start*step, start+cur_index*step, signals[cur_start:cur_index])
            units.append(peak)
        else:
            peak = (chromosome, start+cur_start*step, start+prev*step, signals[cur_start:prev])
            units.append(peak)
            prev = cur_index
            cur_start = cur_index
    return units

def merge_unit(chromosome, start, end, cutoff, units, signals, step=10):
    result_units = []
    merged_unit = units[0]
    for i in range(1, len(units)):
        cur_unit = units[i]
        distance = cur_unit[1] - merged_unit[2]
        if distance > cutoff:
            result_units.append(merged_unit)
            merged_unit = cur_unit
        else:
            merged_unit = (chromosome, merged_unit[1], cur_unit[2],
                           signals[(merged_unit[1]-start)/10:(cur_unit[2]-start)/10])
    result_units.append(merged_unit)
    return result_units




