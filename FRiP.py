import pandas as pd, os

def FRiP(cutoff, path='/home/tmhbxx3/archive/H3K4me3_Encode_wig_splits/code/wigs',peak_path='/home/tmhbxx3/scratch/ENC_H3K4me3/peaks/'):
    total_wig_signals = {}
    wigs = [x[:x.find('.')] for x in os.listdir(path) if x.endswith('.wig')]

    for wig in wigs:
        total_wig_signals[wig] = 5000000000

    FRIP = []

    peak_path = peak_path + str(cutoff) + '/pooled/'
    for wig in wigs:
        peaks_xls = peak_path + wig + '.bgsub.Fnor.peaks_' + str(cutoff) + '.xls'
        df = pd.read_csv(peaks_xls, sep='\t')
        sum_total_signal = df['total_signal'].sum()
        # print df10.columns
        FRIP.append((wig, sum_total_signal / total_wig_signals[wig], df['width_above_cutoff'].sum(), df.shape[0]))

        print 'complete', wig

    result_df = pd.DataFrame(FRIP)

    result_df.columns = ['sample_name', 'FRiP', 'total_width', 'peak_number']

    result_df.to_csv('./FRIP_csv/FRIP_ENC369_with_input_'+str(cutoff)+'_cutoff.csv', index=None)
    return result_df

def FRiP_correlation(cutoff, method='spearman', path='./FRIP_csv/'):
    file_path = path + 'FRIP_ENC369_with_input_' + str(cutoff) + '_cutoff.csv'
    df = pd.read_csv(file_path, index_col=0)
    result = df.corr(method=method)
    # print result
    # Frip_vs_total_width, Frip_vs_peaknumber, total_width_vs_peak_number
    return [result.ix[0, 1], result.ix[0, 2], result.ix[1, 2]]

def FRiP_mean(cutoff, path='./FRIP_csv/'):
    file_path = path + 'FRIP_ENC428_with_input_'+str(cutoff)+'_cutoff.csv'
    df = pd.read_csv(file_path, index_col=0)
    return df['FRiP'].mean()

if __name__ == "__main__":
    for cutoff in [3,6] + range(10,310,10):
        FRiP(cutoff)
    #
    # results = []
    # for cutoff in [3, 6] + range(10, 310, 10):
    #     cur_result = FRiP_mean(cutoff)
    #     results.append((cutoff, cur_result))
    #
    # df = pd.DataFrame(results)
    # df.columns = ['cutoff', 'FRiP']
    #
    # df.to_csv('FRiP_vs_cutoff.csv', index=None)
    results = []
    for cutoff in [3, 6] + range(10, 310, 10):
        results.append(FRiP_correlation(cutoff))
    df = pd.DataFrame(results)
    df.columns = ['FrIP_vs_total_width', 'FrIP_vs_peaknumber', 'total_width_vs_peak_number']
    df.to_csv('FRiP_correlation.csv', index=None)
