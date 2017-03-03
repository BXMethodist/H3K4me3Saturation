
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

        split_indexes = {}
        for i in range(len(self.variants)):
            variant = self.variants[i]
            split_indexes[i] = variant.split_units()
            print split_indexes, i, "first split indexes"


        split_indexes = self.merge_split_index(split_indexes)

        print split_indexes

        for i in range(len(self.variants)):
            split_index = split_indexes[i]
            variant = self.variants[i]
            variant.split_on_common_index(split_index)

        self.representatives = np.asarray([variant.signals for variant in self.variants])
        self.seeds = np.asarray([variant.seed for variant in self.variants])
        self.labels = [variant.labels for variant in self.variants]

        # for variant in self.variants:
        #     print [(unit.start, unit.end) for unit in variant.units]

        # Does it worth to be plot and check
        self.plot = self.plotable()

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

        print group_units_by_overlap
        for group in group_units_by_overlap:
            # print [(g[0]-self.start, g[1]-self.start, g[2]) for g in group]
            related_variants = set([x[2] for x in group])
            if len(related_variants) == 1:
                continue
            left_border = (min([x[0] for x in group])-self.start)/self.step
            right_border = (max([x[1] for x in group])-self.start)/self.step

            # if (right_border - left_border)*self.step > (self.end-self.start)/3.0:
            #     continue

            merge_boo = True
            for i in range(len(self.variants)):
                units_in_cur_group = [u for u in group if u[2] == i]
                if any([u[3] for u in units_in_cur_group]):
                    merge_boo = False
                    break
                for j in range(i, len(self.variants)):
                    units_in_compare_group = [u for u in group if u[2] == j]
                    overlap_region = overlap_length(units_in_cur_group, units_in_compare_group)

                    if (total_length(units_in_cur_group) != 0 and total_length(units_in_compare_group) != 0) and \
                            (overlap_region/total_length(units_in_cur_group) < 0.6 or
                                         overlap_region/total_length(units_in_compare_group) < 0.6):
                        print overlap_region/total_length(units_in_cur_group), \
                            overlap_region / total_length(units_in_compare_group)
                        merge_boo = False

            if merge_boo:
                left_start = left_border * self.step + self.start
                right_end = right_border * self.step + self.start
                print group, "merge in progress"
                for v in range(len(self.variants)):
                    affected_units = [unit[0:2] for unit in group if unit[2] == v]
                    # print affected_units, "affected"
                    cur_variant = self.variants[v]
                    cur_variant.units = [unit for unit in cur_variant.units
                                         if (unit.start, unit.end) not in affected_units]

                    new_unit = Unit(self.chromosome, left_start, right_end,
                                    cur_variant.signals[left_border:right_border])
                    cur_variant.units.append(new_unit)
                    cur_variant.units = sorted(cur_variant.units, key=lambda x: x.start)

            for i in range(len(self.variants)):
                cur_variant = self.variants[i]
                for j in range(len(cur_variant.units)):
                    cur_unit = cur_variant.units[j]
                    max_overlap = 0
                    max_index = None
                    for k in range(i+1, len(self.variants)):
                        target_variant = self.variants[k]
                        for m in range(len(target_variant.units)):
                            target_unit = target_variant.units[m]
                            overlap_region_length = overlap((cur_unit.start, cur_unit.end),
                                                            (target_unit.start, target_unit.end))
                            overlap_percentage = min(overlap_region_length/(target_unit.end-target_unit.start),
                                                     overlap_region_length/(cur_unit.end-cur_unit.start))
                            # print overlap_percentage, i, j, k, m
                            if overlap_percentage > max_overlap:
                                max_overlap = overlap_percentage
                                max_index = (k, m)
                    # print max_overlap, " max overlap"
                    if 1.0 > max_overlap > 0.8:
                        k, m = max_index
                        target_unit = self.variants[k].units[m]
                        merge_start = min(cur_unit.start, target_unit.start)
                        merge_end = max(cur_unit.end, target_unit.end)
                        cur_unit.start = merge_start
                        cur_unit.end = merge_end
                        target_unit.start = merge_start
                        target_unit.end = merge_end
                        cur_unit.signals = self.variants[j].signals[
                            int((merge_start - self.variants[j].start)/self.variants[j].step):
                            int((merge_end - self.variants[j].start) / self.variants[j].step)]
                        target_unit.signals = self.variants[k].signals[
                                           int((merge_start - self.variants[k].start) / self.variants[k].step):
                                           int((merge_end - self.variants[k].start) / self.variants[k].step)]
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
            print cur_pair
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
        need_to_remove = set()

        frontiers = [unit]
        while len(frontiers) > 0:
            cur_unit = frontiers.pop(0)
            left_border = cur_unit[0]
            right_border = cur_unit[1]
            for u in self.units:
                if left_border < u[0] < right_border or \
                    left_border < u[1] < right_border or \
                    u[0] <  left_border < u[1] or \
                    u[0] < right_border < u[1]:
                    if u not in overlap_units:
                        overlap_units.add(u)
                        frontiers.append(u)
                        need_to_remove.add(u)
        for u in need_to_remove:
            self.units.remove(u)
        return overlap_units

    def plotable(self):
        if len(self.variants) < 2:
            return False
        if self.end - self.start < 2000:
            return False
        if all(len(variant.units)==1 for variant in self.variants):
            return False
        return True

    def merge_split_index(self, splited_indexes):
        new_split_indexes = {}
        for i in range(len(splited_indexes.keys())):
            new_split_indexes[i] = set()
        for i in range(len(splited_indexes.keys())):
            cur_indexes = splited_indexes[i]
            for index in cur_indexes:
                related = set()
                for j in range(len(splited_indexes.keys())):
                    target_indexes = splited_indexes[j]
                    for target_index in target_indexes:
                        if abs(target_index - index) < 15:
                            related.add((target_index, j))
                if len(related) == 0:
                    new_split_indexes[i].add(int(index))
                else:
                    mean_index = int(np.mean([x[0] for x in related]))
                    related_variants = [x[1] for x in related]
                    for v in related_variants:
                        new_split_indexes[v].add(mean_index)
        return new_split_indexes


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

        units= callunitbycutoff(self.chromosome, self.start, self.end, self.signals,
                                          cutoff=self.cutoff, step=self.step)
        merged_units = merge_unit(self.chromosome, self.start, self.end, len(self.signals) * self.step / 20.0,
                                  units, self.signals, self.step)

        merged_units = sorted(merged_units, key=lambda x: x[1])
        self.units = []
        for unit in merged_units:
            cur_chr, cur_start, cur_end, cur_signals = unit
            self.units.append(Unit(cur_chr, cur_start, cur_end, cur_signals))

    def similar(self, other):
        if len(self.units) != len(other.units):
            return False
        for i in range(len(self.units)):
            unit = self.units[i]
            other_unit = other.units[i]
            if abs(unit.start - other_unit.start) <= 200 and abs(unit.end - other_unit.end) <= 200:
                continue
            else:
                return False
        return True

    def split_units(self):
        split_indexes = []
        frontiers = [(unit.chromosome, unit.start, unit.end, unit.signals, False) for unit in self.units]

        while len(frontiers) > 0:
            cur_chromosome, cur_start, cur_end, cur_signals, splitted = frontiers.pop(0)
            split_index = splitable(cur_chromosome, cur_start, cur_end, cur_signals, self.convex_cutoff)
            if split_index is None:
                continue
            else:
                left = (cur_chromosome, cur_start, cur_start + split_index * self.step, cur_signals[:split_index], True)
                right = (cur_chromosome, cur_start + split_index * self.step, cur_end, cur_signals[split_index:], True)
                frontiers.append(left)
                frontiers.append(right)
                split_indexes.append((cur_start - self.start)/self.step+split_index)

        split_indexes = sorted(split_indexes)
        return split_indexes

    def split_on_common_index(self, splitted_indexes):
        new_units = []
        split = False
        splitted_indexes = sorted(splitted_indexes)

        unit_split_map = {}

        for unit in self.units:
            related_index = set()
            for index in splitted_indexes:
                if unit.start < self.start+index*self.step < unit.end:
                    related_index.add(index)
            unit_split_map[unit] = related_index

        for key, value in unit_split_map.items():
            if len(value) == 0:
                new_units.append(key)
            else:
                cur_start = key.start
                value = sorted(list(value))
                for index in value:
                    left = Unit(self.chromosome, cur_start, self.start + index * self.step,
                                self.signals[int((cur_start-key.start)/self.step):index])
                    new_units.append(left)
                    cur_start = self.start + index * self.step
                right = Unit(self.chromosome, cur_start, key.end,
                             self.signals[value[-1]:int((key.end-self.start)/self.step)])
                new_units.append(right)

        self.units = new_units
        return


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

    split_indexes = None
    minimum_height = float('inf')

    for min in split_local_mins_index:
        left_part = signals[:left+min]
        right_part = signals[left+min:]

        left_unit_height = np.max(left_part)
        right_unit_height = np.max(right_part)
        min_height = split_region[min]

        # print left + min, min_height
        # print left_peak_height, right_peak_height, min_height
        if left_unit_height > 2 * min_height and right_unit_height > 2 * min_height and \
                (left_unit_height >= convex_cutoff and right_unit_height >= convex_cutoff):
            if split_indexes is None:
                split_indexes = left+min
                minimum_height = min_height
            elif min_height < minimum_height:
                split_indexes = left + min
                minimum_height = min_height
            print "split", split_indexes
    if split_indexes is None:
        return None
    else:
        return int(split_indexes)

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
                if (cur_index-cur_start)/(end-start) > 0.1 or \
                        np.max(signals[cur_start:cur_index]) > np.max(signals)*0.15:
                    units.append(peak)
        else:
            if cur_start != prev:
                peak = (chromosome, start+cur_start*step, start+prev*step, signals[cur_start:prev])
                if (cur_index - cur_start) / (end - start) > 0.1 or \
                                np.max(signals[cur_start:cur_index]) > np.max(signals) * 0.15:
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



