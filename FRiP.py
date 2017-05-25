import pandas as pd, os

def FRiP(cutoff, path='./wigs',peak_path='/home/tmhbxx3/archive/KFH3K4me3/'):
    total_wig_signals = {}
    wigs = [x[:-4] for x in os.listdir(path) if x.endswith('.pkl') and x[:-4] != '']
    # print wigs
    for wig in wigs:
        # wig_pkl = load_obj('./wigs/'+wig)
        # total_signals=0
        # for key, value in wig_pkl.genome.iteritems():
        #     total_signals += np.sum(value.signals)
        # total_wig_signals[wig] = total_signals
        # print total_signals
        total_wig_signals[wig] = 5000000000

    FRIP10 = []

    peak_path = peak_path + str(cutoff) + 'cutoff/pooled/'
    for wig in wigs:
        peaks_xls10 = peak_path + 'archive_tmhkxc48_BroadH3K4me3_broadpeak201401_H3K4me3_dregion_pooled_' + wig + '.peaks.xls'
        df10 = pd.read_csv(peaks_xls10, sep='\t')
        sum10 = df10['total_signal'].sum()
        # print df10.columns
        FRIP10.append((wig, sum10 / total_wig_signals[wig], df10['width_above_cutoff'].sum(), df10.shape[0]))

        print 'complete', wig

    result_df10 = pd.DataFrame(FRIP10)

    result_df10.columns = ['sample_name', 'FRiP', 'total_width', 'peak_number']

    result_df10.to_csv('FRIP_336_'+str(cutoff)+'_cutoff.csv', index=None)
    return

def FRiP_correlation(cutoff, method='spearman', path='/Users/boxia/Desktop/reference_map/QC/FRIP/'):
    file_path = path + 'FRIP_336_'+str(cutoff)+'_cutoff.csv'
    df = pd.read_csv(file_path, index_col=0)
    result = df.corr(method=method)
    return [result.ix[0, 1], result.ix[0, 2], result.ix[1, 2]]

def FRiP_mean(cutoff, path='/Users/boxia/Desktop/reference_map/QC/FRIP/'):
    file_path = path + 'FRIP_336_'+str(cutoff)+'_cutoff.csv'
    df = pd.read_csv(file_path, index_col=0)
    return df['FRiP'].mean()

results = []
for cutoff in [3, 6] + range(10, 310, 10):
    cur_result = FRiP_mean(cutoff)
    results.append((cutoff, cur_result))

df = pd.DataFrame(results)
df.columns = ['cutoff', 'FRiP']

df.to_csv('FRiP_vs_cutoff.csv', index=None)

