

# This class is used to create reference map, not the final version of referencemap for annotation.
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
    def __init__(self, chromosome, start, end, representatives,
                 variants_members, seeds, labels, sample_names,
                 step=10):
        self.chromosome = chromosome,
        self.start = start
        self.end = end
        self.representatives = representatives
        self.seeds = seeds
        self.labels = labels
        self.sample_names = sample_names
        self.step = step
        self.variants = []
        self.units = set()
        self.variants_members = variants_members

        self.create_variants(variants_members)
        self.merge_units_across_variants()
        self.merge_variants()

        split_indexes = {}
        for i in range(len(self.variants)):
            variant = self.variants[i]
            split_indexes[i] = variant.split_units()
            # print split_indexes, i, "first split indexes"


        split_indexes = self.merge_split_index(split_indexes)

        # print split_indexes

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
        self.transitions = self.type_transition()
        # print self.transitions
        # print self.chromosome, self.start, self.end

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

        #   group_units_by_overlap
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
                            (overlap_region) < 0.6 * total_length(units_in_cur_group or
                                         overlap_region < 0.6 * total_length(units_in_compare_group)):
                        # print overlap_region/total_length(units_in_cur_group), \
                        #     overlap_region / total_length(units_in_compare_group)
                        merge_boo = False

            if merge_boo:
                left_start = left_border * self.step + self.start
                right_end = right_border * self.step + self.start
                # print group, "merge in progress"
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
                        cur_unit.signals = self.variants[i].signals[
                            int((merge_start - self.variants[i].start)/self.variants[i].step):
                            int((merge_end - self.variants[i].start) / self.variants[i].step)]
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
            # print cur_pair
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
        # if self.end - self.start < 2000:
        #     return False
        # if all(len(variant.units)==1 for variant in self.variants):
        #     return False
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

    def type_transition(self):
        """
        This function is used to create map between each variant and store the transfer type between them.
        Available type includes: Broadness to Narrow (BN), Convex to Concave (CC), Shift (SH), Pattern (PT),
        and Shape or Other (SO).
        :return: a map, key is the (i, j), i, j are index of variant in self.variants. value are type, which is a string.
        """

        transition = {}
        if len(self.variants) == 1:
            return transition

        for i in range(len(self.variants)):
            for j in range(i+1, len(self.variants)):
                self.variants[i].units = sorted(self.variants[i].units, key=lambda x:x.start)
                self.variants[j].units = sorted(self.variants[j].units, key=lambda x: x.start)
                type = self.get_types(self.variants[i], self.variants[j])
                transition[(i, j)] = type
                transition[(j, i)] = type
        return transition

    def get_types(self, variant1, variant2):
        types = []
        if isShift(variant1, variant2):
            types.append('SH')
        if isBroadNarrow(variant1, variant2):
            types.append('BN')
        if isConvexConcave(variant1, variant2):
            types.append('CC')
        if len(types) == 0 and isPattern(variant1, variant2):
            types.append('PT')
        if len(types) == 0:
            types.append('SO')
        if ('SH' in types and 'CC' in types) or ('SH' in types and 'BN' in types):
            print variant1.chromosome, variant1.start, variant1.end
            print variant2.chromosome, variant2.start, variant2.end
        return types


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
        self.convex_cutoff = np.max(signals)*2.0/3.0
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

    def best_representative_sample(self):
        correlations = [np.corrcoef(self.signals, member)[0][1] for member in self.members]
        corr_ranks = [(correlations[i], i) for i in range(len(correlations)) ]
        corr_ranks = sorted(corr_ranks, key=lambda x:x[0], reverse=True)
        return [self.labels[rank[1]] for rank in corr_ranks[0:3]]


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
        if left_unit_height > 3 * min_height and right_unit_height > 3 * min_height and \
                (left_unit_height >= convex_cutoff and right_unit_height >= convex_cutoff):
            if split_indexes is None:
                split_indexes = left+min
                minimum_height = min_height
            elif min_height < minimum_height:
                split_indexes = left + min
                minimum_height = min_height
            # print "split", split_indexes
    if split_indexes is None:
        return None
    else:
        return int(split_indexes)

def callunitbycutoff(chromosome, start, end, signals, cutoff, step=10):
    units = []
    # print signals
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
                if (cur_index-cur_start)*step > 0.1 * (end-start) or \
                        np.max(signals[cur_start:cur_index]) > np.max(signals)*0.15:
                    units.append(peak)
        else:
            if cur_start != prev:
                peak = (chromosome, start+cur_start*step, start+prev*step, signals[cur_start:prev])
                if (cur_index - cur_start)*step > 0.1 * (end - start) or \
                                np.max(signals[cur_start:cur_index]) > np.max(signals) * 0.15:
                    units.append(peak)
                prev = cur_index
            cur_start = cur_index
    # print [x[0:3] for x in units]
    return units

def merge_unit(chromosome, start, end, cutoff, units, signals, step=10):
    result_units = []
    # print units
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

def isShift(variant1, variant2):
    """
    :param variant1: Class Variant obj
    :param variant2: Class Variant obj
    :return: boolean, whether these two 's relation are shift
    """
    # if len(variant1.units) != len(variant2.units):
    #     return False

    variant1_max_index = np.argmax(variant1.signals)
    variant2_max_index = np.argmax(variant2.signals)

    distance = abs(variant1_max_index-variant2_max_index)
    if variant1_max_index < variant2_max_index:
        rolling_variant1_signals = np.roll(variant1.signals, distance)
        rolling_variant1_signals = np.append(rolling_variant1_signals, rolling_variant1_signals[:distance])
        rolling_variant1_signals[:distance] = 0
        if np.corrcoef(rolling_variant1_signals, np.append(variant2.signals, [0]*distance))[0 ,1] > 0.95:
            return True
    else:
        rolling_variant2_signals = np.roll(variant2.signals, distance)
        rolling_variant2_signals = np.append(rolling_variant2_signals, rolling_variant2_signals[:distance])
        rolling_variant2_signals[:distance] = 0
        if np.corrcoef(rolling_variant2_signals, np.append(variant1.signals, [0]*distance))[0 ,1] > 0.95:
            return True
    return False

def isBroadNarrow(variant1, variant2):
    """
    logic, units length of v1 > 2* v2 and the bigger one need to be at least 0.6 of the length of the variant.
    the distance of the two variant submit need to be close enough.
    :param variant1: Class Variant obj
    :param variant2: Class Variant obj
    :return: boolean, whether these two's relation are Broad to Narrow
    """
    total_width = abs(variant1.end - variant1.start)
    variant1_max_index = variant1.start + np.argmax(variant1.signals)*variant1.step
    variant2_max_index = variant2.start + np.argmax(variant2.signals) * variant2.step

    v1_total_units_width = units_total_length(variant1.units, [x for x in range(len(variant1.units))])
    v2_total_units_width = units_total_length(variant2.units, [x for x in range(len(variant2.units))])

    if v1_total_units_width > 2*v2_total_units_width or v2_total_units_width > 2 * v1_total_units_width:
        pass
    else:
        return False

    if abs(variant1_max_index-variant2_max_index) > total_width*0.25:
        return False


    max_start = None
    max_end = None
    for unit in variant1.units:
        if max_start is None:
            max_start = unit.start
            max_end = unit.end
        elif unit.start == max_end:
            max_end = unit.end
        elif max_start <= variant1_max_index <= max_end:
            break
        else:
            max_start = unit.start
            max_end = unit.end
    if max_start <= variant1_max_index <= max_end:
        v1_max_width = max_end - max_start
    else:
        print variant1_max_index, max_end, max_start

    max_start = None
    max_end = None
    for unit in variant2.units:
        if max_start is None:
            max_start = unit.start
            max_end = unit.end
        elif unit.start == max_end:
            max_end = unit.end
        elif max_start <= variant2_max_index <= max_end:
            break
        else:
            max_start = unit.start
            max_end = unit.end
    if max_start <= variant2_max_index <= max_end:
        v2_max_width = max_end - max_start
    else:
        print variant2_max_index, max_end, max_start

    if v1_max_width * 4 <= v2_max_width and v2_max_width >= 0.6* total_width:
        return True
    if v2_max_width * 4 <= v1_max_width and v1_max_width >= 0.6* total_width:
        return True
    return False

def units_total_length(units, indexes):
    """
    :param units: variant unit object
    :return: int, total units length
    """
    sum = 0
    for i in indexes:
        unit = units[i]
        sum += abs(unit.end - unit.start)
    return sum

def isConvexConcave(variant1, variant2):
    """
    To check whether their relationship is convex to concave
    :param variant1: Variant Obj
    :param variant2: Variant Obj
    :return: Boolean
    """
    concave1 = isConcave(variant1)
    concave2 = isConcave(variant2)

    # print concave1, concave2

    if (concave1[0] and concave2[0]):
        return False
    elif (not concave1[0] and not concave2[0]):
        return False

    if (concave1[0]):
        left, right = concave1[1], concave1[2]
        mid = concave2[1]
    elif (concave2[0]):
        left, right = concave2[1], concave2[2]
        mid = concave1[1]

    left_submit = np.argmax(left.signals)*left.step + left.start
    right_submit = np.argmax(right.signals)*right.step + right.start
    mid_submit = np.argmax(mid.signals)*mid.step + mid.start

    if not left_submit < mid_submit < right_submit:
        return False
    break_point = (right.start + left.end)/2

    distance = min(abs(left_submit - mid_submit), abs(right_submit - mid_submit))

    if abs(mid_submit - break_point) < distance:
        return True
    else:
        return False

def isPattern(variant1, variant2):
    """
    diffirent number of nunits is a indicator for pattern change, or units have different locations
    :param variant1: Class Variant Obj
    :param variant2: Class Variant Obj
    :return: Boolean
    """
    if len(variant1.units) != len(variant2.units):
        return True

    if len(variant1.units) == 1:
        return False
    return True

def isConcave(variant):
    """
    :param variant: Class Variant obj
    :return: True or False
    """
    if len(variant.units) <= 1:
        return False, variant.units[0], None

    unit_indexes = [x for x in range(len(variant.units))]

    unit_indexes = sorted(unit_indexes, key=lambda x: variant.units[x].height, reverse=True)
    unit_indexes = unit_indexes[0:2]
    unit_indexes = sorted(unit_indexes)

    left, right = variant.units[unit_indexes[0]], variant.units[unit_indexes[1]]

    total_width = variant.end - variant.start
    left_submit = np.argmax(left.signals) * left.step + left.start
    right_submit = np.argmax(right.signals) * right.step + right.start

    if right_submit - left_submit > total_width *0.5:
        max_unit = left if left.height > right.height else right
        return False, max_unit, None

    if left.height > variant.convex_cutoff and \
                    right.height > variant.convex_cutoff and \
                                    (right.start - left.end) < min(right.end-right.start, left.end-left.start):
        return True, left, right
    else:
        max_unit = left if left.height > right.height else right
        return False, max_unit, None
