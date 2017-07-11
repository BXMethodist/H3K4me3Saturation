import os, pandas as pd
import Wig


if __name__ == "__main__":
    # wigPath = "/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/wigs/"
    #
    # wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]
    # wig_list = open('mapping75_NSC1.05_RSC0.8_.txt','r')
    # wigs = [x.strip() for x in wig_list.readlines()]
    # wig_list.close()
    # works = set()
    #
    # pair_df = pd.read_csv('ENC_H3K4me3_sample_pairs.csv', index_col=0)
    # pair_df = pair_df[pd.notnull(pair_df['input'])]
    #
    # archive_df = pd.read_csv('archived_sample_input_pair.csv', index_col=0)
    # archive_df = archive_df[pd.notnull(archive_df['input'])]
    #
    # candidates = list(pair_df.index) + list(archive_df.index)
    #
    # for wig_file in wigFiles:
    #     name = wig_file[:wig_file.find('.bgsub.Fnor.wig')]
    #     if name in wigs and name in candidates:
    #         # print wig_file
    #         works.add(name)
    #
    #         path = wigPath + wig_file
    #         # os.system('cp ' + path + ' ./')
    #     if name in wigs and name not in candidates:
    #         print name
    # for w in wigs:
    #     if w not in works:
    #         print w

    works = [x for x in os.listdir('./') if x.endswith(".wig")]

    for work in works[0:100]:
        wig = Wig.Wig(work)
        wig.save_split_chr(10000)

