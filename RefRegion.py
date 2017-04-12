# This module takes the reference map objects from RegionCluster and transform it to a light weighted map for CallRegion.

import pickle
import pandas as pd
from collections import defaultdict
from clusterUtils import mean_normalization

class ReferenceRegion():
    def __init__(self, region):
        """
        :param region: region class from construction of reference map
        """
        self.chromosome = str(region.chromosome)
        self.chromosome = self.chromosome[self.chromosome.find("'")+1:self.chromosome.rfind("'")]
        self.start = region.start
        self.end = region.end
        self.member_size = len(region.variants_members)
        self.plot = region.plotable()
        self.id = '.'.join([self.chromosome[3:], str(self.start), str(self.end)])

        self.setVariants(region.variants)

        self.transitions = {}
        for key, value in region.transitions.iteritems():
            id1, id2 = key
            id1 = self.variants[id1].id
            id2 = self.variants[id2].id
            self.transitions[(id1, id2)] = value
        categories = []
        for value in self.transitions.values():
            categories += value
        self.categories = list(set(categories))

    def setVariants(self, variants):
        self.variants = []
        for i in range(len(variants)):
            variant = variants[i]
            self.variants.append(ReferenceVariant(variant, self, i))


class ReferenceVariant():
    def __init__(self, variant, parent, i):
        self.parent = parent
        self.chromosome = str(variant.chromosome)
        self.chromosome = self.chromosome[self.chromosome.find("'") + 1:self.chromosome.rfind("'")]
        self.start = variant.start
        self.end = variant.end
        self.representative = variant.signals
        self.id = '.'.join([self.chromosome[3:], str(self.start), str(self.end), str(i+1)])

        self.setUnits(variant.units)
        self.center = variant.center
        self.left_boundary = variant.left_boundary
        self.right_boundary = variant.right_boundary

    def setUnits(self, units):
        self.units = []
        for j in range(len(units)):
            unit = units[j]
            self.units.append(ReferenceUnit(unit, self, j))

class ReferenceUnit():
    def __init__(self, unit, parent, j):
        self.parent = parent
        self.chromosome = str(unit.chromosome)
        self.chromosome = self.chromosome[self.chromosome.find("'") + 1:self.chromosome.rfind("'")]
        self.start = unit.start
        self.end = unit.end
        self.id = parent.id + '.' + str(self.start) + '.' + str(self.end) + '.' + str(j+1)


def Annotation(path, output):
    """
    :param path: pickle file location from Region Cluster
    :return: a region annotation file
            colnames: chromosome, start, end, id
            id consist: region_id: chr_start_end, variant_id: chr_start_end_vid, unit_id:chr_start_end_uid
    """
    with open(path, 'rb') as f:
        regions = pickle.load(f)
    f.close()

    new_regions = []
    region_annotations = []
    variant_annotations = []
    units_annotations = []
    count = 0
    counts = defaultdict(int)

    stats = []

    for region in regions:
        region.representatives = mean_normalization(region.representatives)
        referenceRegion = ReferenceRegion(region)
        if referenceRegion.plot:
            count += 1
        counts[len(referenceRegion.variants)]+=1
        BN = False
        PT = False
        SO = False
        SH = False
        CC= False
        if 'BN' in referenceRegion.categories:
            BN = True
        if 'PT' in referenceRegion.categories:
            PT = True
        if 'SO' in referenceRegion.categories:
            SO = True
        if 'SH' in referenceRegion.categories:
            SH = True
        if 'CC' in referenceRegion.categories:
            CC = True

        cur_stat = [referenceRegion.id,
                    referenceRegion.end-referenceRegion.start,
                    len(referenceRegion.variants),
                    referenceRegion.member_size,
                    BN,
                    CC,
                    SH,
                    PT,
                    SO]
        stats.append(cur_stat)

        new_regions.append(referenceRegion)

        row = [referenceRegion.chromosome, str(referenceRegion.start), str(referenceRegion.end)]
        region_id = 'region_id:'+referenceRegion.id+';'
        row.append(region_id)
        region_annotations.append(row)

        for i in range(len(referenceRegion.variants)):
            cur_variant = referenceRegion.variants[i]
            row = [cur_variant.id, referenceRegion.id, cur_variant.chromosome, cur_variant.start, cur_variant.end]
            cur_changes = set()
            cur_changes_map = defaultdict(set)
            for key, value in region.transitions.items():
                if cur_variant.id in key:
                    cur_changes.add(value)
                    key1, key2 = key
                    if key1 == cur_variant.id:
                        cur_key = key2
                    elif key2 == cur_variant.id:
                        cur_key = key1
                    cur_changes_map[value].add(cur_key)
            if 'BN' in cur_changes:
                row += [cur_changes_map['BN']]
            else:
                row += [False]
            if 'CC' in cur_changes:
                row += [cur_changes_map['CC']]
            else:
                row += [False]
            if 'SH' in cur_changes:
                row += [cur_changes_map['SH']]
            else:
                row += [False]
            if 'PT' in cur_changes:
                row += [cur_changes_map['PT']]
            else:
                row += [False]
            if 'SO' in cur_changes:
                row += [cur_changes_map['SO']]
            else:
                row += [False]

            # row.append(cur_variant.representative)
            variant_annotations.append(row)

            for j in range(len(cur_variant.units)):
                cur_unit = cur_variant.units[j]
                row = [cur_unit.chromosome, str(cur_unit.start), str(cur_unit.end)]
                cur_unit_id = 'unit_id:'+cur_unit.id+';'
                row.append(cur_unit_id)
                units_annotations.append(row)

    region_df = pd.DataFrame(region_annotations)
    variant_df = pd.DataFrame(variant_annotations, columns=['variant_id', "region_id", 'chromosome', 'start', 'end',
                                                            'BroadtoNarrow', 'ConcavetoConvex', 'Shift', 'Pattern',
                                                            'Other'])
    variant_df = variant_df.set_index(['variant_id'])
    units_df = pd.DataFrame(units_annotations)
    stats_df = pd.DataFrame(stats, columns=['region_id', 'width', 'number of clusters', 'number of samples',
                                            'BroadtoNarrow', 'ConcavetoConvex', 'Shift', 'Pattern', 'Other'])
    stats_df = stats_df.set_index(['region_id'])

    region_df.to_csv(output+'region.tsv', sep='\t', index=False, header=False)
    variant_df.to_csv(output + 'variant.tsv', sep='\t', index=False, header=False)
    units_df.to_csv(output + 'units.tsv', sep='\t', index=False, header=False)
    stats_df.to_csv(output + 'stats.tsv', sep='\t')

    with open(output + '.pkl', 'wb') as f:
        pickle.dump(new_regions, f, pickle.HIGHEST_PROTOCOL)
    f.close()

    print "total ", count, 'region need to be plot'
    print counts
    return region_annotations, variant_annotations, units_annotations

# if __name__ == "__main__":
# import os
# pkls = os.listdir("./pkl_parts")
#
# regions = []
#
# for pkl in pkls:
#     with open("./pkl_parts/" + pkl, 'rb') as f:
#         region = pickle.load(f)
#     regions += region
#     f.close()
#
# import pickle
#
# with open('75refmap_combined_3kb_regions' + '.pkl', 'wb') as f:
#     pickle.dump(regions, f, pickle.HIGHEST_PROTOCOL)
#
# f.close()
#
# Annotation("75refmap_combined_3kb_regions.pkl", "./pkl/75_combined_3kb")