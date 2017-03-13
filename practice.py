# import os
#
# def runFastqdump(list):
#     f = open(list, 'r')
#     ids = [x.strip() for x in f.readlines()]
#     f.close()
#
#     for link in ids:
#         ids = link.split()
#         id1 = ids[0][ids[0].rfind('/')+1:-9]
#         id2 = ids[1][ids[1].rfind('/')+1:-9]
#         print id1, id2
#         pbs = open(id1+'_'+id2+".pbs", "w")
#         pbs.write("#!/bin/bash\n")
#         pbs.write("#PBS -r n\n")
#         pbs.write("#PBS -N "+id1+'_'+id2+'\n')
#         pbs.write("#PBS -q mediummem\n")
#         pbs.write("#PBS -m e\n")
#         pbs.write("#PBS -M bxia@houstonmethodist.org\n")
#         pbs.write("#PBS -l walltime=96:00:00\n")
#         pbs.write("#PBS -l nodes=1:ppn=8\n")
#         pbs.write("#PBS -l pmem=16000mb\n")
#         pbs.write("cd /home/tmhbxx3/archive/H3K4me3/ENCODE_with_input/FASTQ\n")
#         pbs.write("module load python/2.7.11\n")
#         pbs.write('wget '+ids[0] +'\n')
#         pbs.write('wget '+ids[1] + '\n')
#         pbs.write('gunzip ' + ids[0][ids[0].rfind('/')+1:]+'\n')
#         pbs.write('gunzip ' + ids[1][ids[1].rfind('/') + 1:] + '\n')
#
#         pbs.write('cat ' + id1 + '1.fastq ' + id2 + '.fastq >' + id1+'_'+id2 + '.fastq'+'\n')
#
#         pbs.write('bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
#                   + id1+'_'+id2 + ".fastq " + id1+'_'+id2 + ".bowtie\n")
#         pbs.close()
#         os.system('qsub '+id1+'_'+id2+".pbs")
#         break
#     return
#
# list = '/home/tmhbxx3/archive/H3K4me3/ENCODE_with_input/FASTQ/pair.txt'
#
#
# runFastqdump(list)
import pandas as pd
df = pd.read_csv('75_combined_3kbhg38_RefSeq_allgene30003000with_cluster.tsv', index_col=0, sep='\t')

df_result = df[df['TSS'].notnull()]

f = open('clustersgenes.tsv', "w")
for gene in df_result[df_result['number of clusters'] >2]['gene_name2'].unique():
    f.write(gene+'\n')
f.close()