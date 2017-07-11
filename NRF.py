import os, pandas as pd

def unique_counts(i, bowties, path):
    results = []
    for bowtie in bowties[i*100: (i+1)*100]:
        fname = bowtie[:-7]

        os.system('cp '+path+bowtie+' ./')

        bowtie_f = open(path + bowtie, 'r')

        if fname.find('_') != -1:
            fname = '_'.join(sorted(fname.split('_')))

        unique_reads = set()

        for line in bowtie_f:
            line = [x.strip() for x in line.split('\t')]
            index = None

            if line[1] == '-' or line[1] == '+':
                index = 1
            elif line[2] == '-' or line[1] == '+':
                index = 2
            if index is None:
                print 'something wrong!!!'
            else:
                strand = None
                if line[index] == "+":
                    strand = 1
                elif line[index] == '-':
                    strand = -1
                # if (strand, line[3].strip(), line[4].strip()) not in unique_reads:
                #     print (strand, line[3].strip(), line[4].strip())
                unique_reads.add((strand, line[index+1].strip(), line[index+2].strip()))
            # break
        bowtie_f.close()
        # break
        unique_reads_count = len(unique_reads)
        results.append((fname, unique_reads_count))

        print fname, unique_reads_count

        # break
    df = pd.DataFrame(results)
    # print df
    df.columns = ['sample_name', 'unique_reads']

    df.to_csv('NRF' + str(i) + '.csv', index=None)

# df = pd.read_csv('ENC_H3K4me3_sample_pairs.csv', index_col=0)
#
# df = df[pd.notnull(df['input'])]
#
# bowties = [index + '.bowtie' for index in df.index]

bowties = [x for x in os.listdir('/archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/') if x.endswith('.bowtie')]

unique_counts(1, bowties, path='/archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/')