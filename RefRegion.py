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
        cur_stat = [referenceRegion.id,
                    referenceRegion.end-referenceRegion.start,
                    len(referenceRegion.variants),
                    referenceRegion.member_size]
        stats.append(cur_stat)

        new_regions.append(referenceRegion)

        row = [referenceRegion.chromosome, str(referenceRegion.start), str(referenceRegion.end)]
        region_id = 'region_id:'+referenceRegion.id+';'
        row.append(region_id)
        region_annotations.append(row)

        for i in range(len(referenceRegion.variants)):
            cur_variant = referenceRegion.variants[i]
            row = [cur_variant.chromosome, str(cur_variant.start), str(cur_variant.end)]
            cur_variant_id = 'variant_id:'+cur_variant.id+';'
            row.append(cur_variant_id)
            row.append(','.join(map(str, cur_variant.representative)))
            variant_annotations.append(row)

            for j in range(len(cur_variant.units)):
                cur_unit = cur_variant.units[j]
                row = [cur_unit.chromosome, str(cur_unit.start), str(cur_unit.end)]
                cur_unit_id = 'unit_id:'+cur_unit.id+';'
                row.append(cur_unit_id)
                units_annotations.append(row)

    region_df = pd.DataFrame(region_annotations)
    variant_df = pd.DataFrame(variant_annotations)
    units_df = pd.DataFrame(units_annotations)
    stats_df = pd.DataFrame(stats, columns=['region_id', 'width', 'number of clusters', 'number of samples'])
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

