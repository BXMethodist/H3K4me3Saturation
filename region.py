
import numpy as np
import math
from scipy.stats import pearsonr
from clusterUtils import remove_duplicate


class Region():
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
        self.seeds = seeds
        self.labels = labels
        self.step = step
        self.variants = []
        self.units = set()

        self.create_variants(variants_members)
        self.merge_units_across_variants()
        self.merge_variants()

        self.representatives = [variant.signals for variant in self.variants]
        self.seeds = [variant.seed for variant in self.variants]
        self.labels = [variant.labels for variant in self.variants]

    def create_variants(self, variants_members):
        for i in range(len(self.representatives)):
            cur_variant = Variant(self.chromosome, self.start, self.end, self.representatives[i],
                                  variants_members[self.labels[i]], self.seeds[i], self.labels[i])
            self.variants.append(cur_variant)

            for unit in cur_variant.units:
                self.units.add((unit.start, unit.end, i, unit.splitted))
        return

    def merge_units_across_variants(self):
        group_units_by_overlap = []
        while len(self.units) > 0:
            unit = self.units.pop()
            cur_overlap = self.get_overlap(unit)
            group_units_by_overlap.append(cur_overlap)

        for group in group_units_by_overlap:
            # print [(g[0]-self.start, g[1]-self.start, g[2]) for g in group]
            related_variants = set([x[2] for x in group])
            if len(related_variants) == 1:
                continue
            left_border = (min([x[0] for x in group])-self.start)/self.step
            right_border = (max([x[1] for x in group])-self.start)/self.step

            if (right_border - left_border)*self.step > (self.end-self.start)/4.0:
                continue

            merge_boo = True
            for i in range(len(self.variants)):
                units_in_cur_group = [u for u in group if u[2] == i]
                if not all([u[3] for u in units_in_cur_group]):
                    merge_boo = False
                    break
                for j in range(i, len(self.variants)):
                    units_in_compare_group = [u for u in group if u[2] == j]
                    overlap_region = overlap_length(units_in_cur_group, units_in_compare_group)

                    if (total_length(units_in_cur_group) != 0 and total_length(units_in_compare_group) != 0) and \
                            (overlap_region/total_length(units_in_cur_group) < 0.6 or
                                         overlap_region/total_length(units_in_compare_group) < 0.6):
                        merge_boo = False

            if merge_boo:
                left_start = left_border * self.step + self.start
                right_end = right_border * self.step + self.start
                for v in range(len(self.variants)):
                    affected_units = [unit[0:2] for unit in group if unit[2] == v]
                    cur_variant = self.variants[v]
                    cur_variant.units = [unit for unit in cur_variant.units
                                         if (unit.start, unit.end) not in affected_units]


                    new_unit = Unit(self.chromosome, left_start, right_end,
                                    cur_variant.signals[left_border:right_border])
                    cur_variant.units.append(new_unit)
                    cur_variant.units = sorted(cur_variant.units, key=lambda x: x.start)
        return

    def merge_variants(self):
        need_to_merged_pair = set()
        for i in range(len(self.variants)):
            cur_variant = self.variants[i]
            for left_variant in self.variants[i+1:]:
                if cur_variant.similar(left_variant) and np.corrcoef(cur_variant.signals, left_variant.signals)[0][1]>=0.6:
                    need_to_merged_pair.add((cur_variant, left_variant))
        has_been_removed = set()

        need_to_merged_pair = group_pair(list(need_to_merged_pair))

        while len(need_to_merged_pair)>0:
            cur_pair = need_to_merged_pair.pop()
            for variant in cur_pair:
                if variant not in has_been_removed:
                    self.variants.remove(variant)
                    has_been_removed.add(variant)
            cur_pair = list(cur_pair)
            self.variants.append(self.combine_variants(cur_pair))
        return

    def combine_variants(self, group):
        members = group[0].members
        seed = group[0].seed
        labels = list(group[0].labels)
        for i in range(len(group)):
            members = np.concatenate((members, group[i].members), axis=0)
            seed += group[i].seed
            labels += list(group[i].labels)
        rep = np.mean(members, axis=0)
        seed = seed/len(group)
        return Variant(self.chromosome, self.start, self.end, rep, members, seed, np.asarray(labels), step=10)


    def get_overlap(self, unit):
        overlap_units = set()
        overlap_units.add(unit)
        need_to_remove = []
        left_border = unit[0]
        right_border = unit[1]
        for u in self.units:
            if left_border < u[0] < right_border or \
                left_border < u[1] < right_border or \
                u[0] <  left_border < u[1] or \
                u[0] < right_border < u[1]:
                overlap_units.add(u)
                need_to_remove.append(u)
                left_border = min(u[0], left_border)
                right_border = max(u[0], right_border)

        for u in self.units:
            if left_border < u[0] < right_border or \
                left_border < u[1] < right_border:
                overlap_units.add(u)
                if u not in need_to_remove:
                    need_to_remove.append(u)

        for u in need_to_remove:
            self.units.remove(u)
        return overlap_units



class Variant():
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

    def __init__(self, chromosome, start, end, signals, members, seed, labels, step=10):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.cutoff = np.max(signals)*0.1
        self.convex_cutoff = np.max(signals)/3
        self.members = members
        self.seed = seed
        self.labels = labels

        units_by_cutoff= callunitbycutoff(self.chromosome, self.start, self.end, self.signals,
                                          cutoff=self.cutoff, step=self.step)
        # print [(peak[0], peak[1], peak[2]) for peak in peaks_by_cutoff]
        # print peaks_by_cutoff
        merged_units = merge_unit(self.chromosome, self.start, self.end, len(self.signals)*self.step/40.0,
                                  units_by_cutoff, self.signals, self.step)

        merged_units = sorted(merged_units, key=lambda x: x[1])
        self.units = []

        for unit in merged_units:
            cur_chr, cur_start, cur_end, cur_signals = unit
            cur_units = callunit(cur_chr, cur_start, cur_end, cur_signals, self.step, self.convex_cutoff)
            self.units += cur_units

    def similar(self, other):
        if len(self.units) != len(other.units):
            return False
        for i in range(len(self.units)):
            unit = self.units[i]
            other_unit = other.units[i]
            if abs(unit.start - other_unit.start) <= 240 and abs(unit.end - other_unit.end) <= 240:
                continue
            else:
                return False
        return True



class Unit():
    def __init__(self, chromosome, start, end, signals, splitted=False, step=10):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.signals = signals
        self.step = step
        self.width = end - start
        self.splitted = splitted
        if self.width <= 0:
            print self.chromosome, self.start, self.end, "width is 0, why?"
        self.height = np.max(self.signals)


def callunit(chromosome, start, end, signals, step, convex_cutoff):
    units = []

    frontiers = []

    frontiers.append((chromosome, start, end, signals, False))

    while len(frontiers) > 0:
        cur_chromosome, cur_start, cur_end, cur_signals, splitted = frontiers.pop(0)
        split_index = splitable(cur_chromosome, cur_start, cur_end, cur_signals, convex_cutoff)
        if split_index is None:
            units.append((cur_chromosome, cur_start, cur_end, cur_signals, splitted))
        else:
            left = (cur_chromosome, cur_start, cur_start+split_index*step, cur_signals[:split_index], True)
            right = (cur_chromosome, cur_start+split_index*step, cur_end, cur_signals[split_index:], True)
            frontiers.append(left)
            frontiers.append(right)

    units = sorted(units, key=lambda x: x[1])
    units_obj = []
    for p in units:
        chromosome, start, end, signals, splitted = p
        u_obj = Unit(chromosome, start, end, signals, splitted=splitted)
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
            # print "split"
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
            if cur_start != prev:
                peak = (chromosome, start+cur_start*step, start+cur_index*step, signals[cur_start:cur_index])
                units.append(peak)
        else:
            if cur_start != prev:
                peak = (chromosome, start+cur_start*step, start+prev*step, signals[cur_start:prev])
                units.append(peak)
                prev = cur_index
                cur_start = cur_index
    # print [x[0:3] for x in units]
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

def overlap_length(units1, units2):
    sum = 0
    for unit1 in units1:
        for unit2 in units2:
            sum += overlap(unit1, unit2)
    return sum*1.0

def total_length(units):
    sum = 0

    for unit in units:
        sum += unit[1] - unit[0]
    return sum*1.0

def overlap(unit1, unit2):
    max_pos = min(unit1[1], unit2[1])
    min_pos = max(unit1[0], unit2[0])
    if min_pos >= max_pos:
        return 0*1.0
    else:
        return (max_pos - min_pos)*1.0

def group_pair(pairs):
    group_number = {}
    new_groups = set()
    for i in range(len(pairs)):
        cur_pair = pairs[i]
        if cur_pair[0] not in group_number and cur_pair[1] not in group_number:
            group_number[cur_pair[0]] = i
            group_number[cur_pair[1]] = i
        elif cur_pair[0] in group_number:
            group_number[cur_pair[1]] = group_number[cur_pair[0]]
        elif cur_pair[1] in group_number:
            group_number[cur_pair[0]] = group_number[cur_pair[1]]
    for i in range(len(pairs)):
        new_pairs = set()
        for key, value in group_number.items():
            if value == i:
                new_pairs.add(key)
        if len(new_pairs) > 0:
            new_groups.add(tuple(new_pairs))
    return new_groups



