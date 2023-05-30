#!/usr/bin/env python

import pandas as pd
import sys

df1 = pd.read_table(sys.argv[1], names=['READBASE', 'PRIMROSE'])
df2 = pd.read_table(sys.argv[2], names=['READBASE','TAG'])
df = df2.merge(df1,on='READBASE')
df[['PRIMROSE', 'TAG']].to_csv(f'{sys.argv[1]}.TAG.bed', sep="\t", header=False, index=False)
