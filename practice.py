import numpy as np
import pandas as pd


df = pd.read_csv('hg38_RefSeq_allgene.txt', sep='\t')

print df

df = df.set_index("#bin")

df.to_csv('hg38_RefSeq_allgene.txt', sep='\t', index=False)