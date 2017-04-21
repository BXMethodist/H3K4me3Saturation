# import os
#
# def runFastqdump(list):
#     f = open(list, 'r')
#     links = [x.strip() for x in f.readlines()]
#     f.close()
#
#     for link in links:
#         id = link[link.rfind('/')+1:]
#         print id
#         print link
#         pbs = open(id+".pbs", "w")
#         pbs.write("#!/bin/bash\n")
#         pbs.write("#PBS -r n\n")
#         pbs.write("#PBS -N "+id+'\n')
#         # pbs.write("#PBS -q mediummem\n")
#         pbs.write("#PBS -m e\n")
#         pbs.write("#PBS -M bxia@houstonmethodist.org\n")
#         pbs.write("#PBS -l walltime=24:00:00\n")
#         pbs.write("#PBS -l nodes=1:ppn=1\n")
#         # pbs.write("#PBS -l pmem=16000mb\n")
#         pbs.write("cd /home/tmhbxx3/archive/H3K4me3/GEO_with_input/input/FASTQ \n")
#         pbs.write("module load python/2.7.11\n")
#         pbs.write('wget '+link +'\n')
#         # pbs.write('wget '+ids[1] + '\n')
#         # pbs.write('gunzip ' + ids[0][ids[0].rfind('/')+1:]+'\n')
#         # pbs.write('gunzip ' + ids[1][ids[1].rfind('/') + 1:] + '\n')
#
#         # pbs.write('cat ' + id1 + '1.fastq ' + id2 + '.fastq >' + id1+'_'+id2 + '.fastq'+'\n')
#         #
#         # pbs.write('bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
#         #           + id1+'_'+id2 + ".fastq " + id1+'_'+id2 + ".bowtie\n")
#         pbs.close()
#         os.system('qsub '+id+".pbs")
#         # break
#     return
#
# list = './download.txt'
#
#
# runFastqdump(list)
# import pandas as pd
# df = pd.read_csv('75_combined_3kbhg38_RefSeq_allgene30003000with_cluster.tsv', index_col=0, sep='\t')
#
# df_result = df[df['TSS'].notnull()]
#
# f = open('clustersgenes.tsv', "w")
# for gene in df_result[df_result['number of clusters'] >2]['gene_name2'].unique():
#     f.write(gene+'\n')
# f.close()


# import os
#
# def runFastqdump(list, path, bowtie_path):
#     f = open(list, 'r')
#     ids = [x.strip() for x in f.readlines()]
#     f.close()
#
#     cmd = "/home/tmhbxx3/tools/sratoolkit/bin/fastq-dump -O /home/tmhbxx3/archive/H3K4me3/GEO_with_input/sample/FASTQ --split-3 "
#
#     bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
#
#     for name in ids:
#         if name in ids:
#             # print cmd + path +name
#             pbs = open(name+".pbs", "w")
#             pbs.write("#!/bin/bash\n")
#             pbs.write("#PBS -r n\n")
#             pbs.write("#PBS -N bowtie_"+name+'\n')
#             pbs.write("#PBS -q mediummem\n")
#             pbs.write("#PBS -m e\n")
#             pbs.write("#PBS -M bxia@houstonmethodist.org\n")
#             pbs.write("#PBS -l walltime=96:00:00\n")
#             pbs.write("#PBS -l nodes=1:ppn=8\n")
#             pbs.write("#PBS -l pmem=16000mb\n")
#             pbs.write("cd /home/tmhbxx3/archive/H3K4me3/GEO_with_input/bowtie\n")
#             pbs.write("module load python/2.7.11\n")
#             pbs.write("wget https://sra-download.ncbi.nlm.nih.gov/srapub/"+name)
#
#
#             pbs.write(cmd + path + name + '\n')
#             pbs.write(
#                 'cat ' + bowtie_path + name + '_1.fastq ' + bowtie_path + name + '_2.fastq >' + bowtie_path + name + '.fastq' + '\n')
#             pbs.write(bowtie_cmd + bowtie_path + name + ".fastq " + name + ".bowtie\n")
#             pbs.write("rm " + bowtie_path + name + ".fastq\n")
#             pbs.write("rm " + bowtie_path + name + "_1.fastq\n")
#             pbs.write("rm " + bowtie_path + name + "_2.fastq\n")
#             pbs.close()
#             os.system('qsub ' + name + ".pbs")
#             break
#         return
#
# list = '/home/tmhbxx3/archive/H3K4me3/GEO_with_input/pairSRR.txt'
# path = '/home/tmhbxx3/archive/H3K4me3/GEO_with_input/sample/SRA/'
# bowtie_path = '/home/tmhbxx3/archive/H3K4me3/GEO_with_input/sample/FASTQ/'
#
# runFastqdump(list, path, bowtie_path)
# import pandas as pd, numpy as np
#
# a = [[1,2,3],[4,5,6]]
#
# df = pd.DataFrame(a, columns=['A','B', 'C'])
# print df
# print df[['B']].mean(axis=1)

# print np.corrcoef(np.asarray([4,2]),np.asarray([-3,-1]))[0,1]
# f= open("hg19_chr_sizes.txt", "r")
#
# info = f.readlines()
# sizes = {}
# for line in info:
#     line = line.strip().split()
#     sizes[line[0]] = line[1]
#
# print sizes

a =[0,1,2,3,4,5]
b = 6
from bisect import bisect

print bisect(a, b)
print a

# if 1 <=2 <=3 <= 4:
#     print 'lalaa'