#!/usr/bin/env python3
import pandas as pd
import sys

N=len(sys.argv)
df== pd.read_table(sys.argv[i], header=0)
for i in range(2,N-1):
    tmp = pd.read_table(sys.argv[i], header=0)
    df = df.merge(tmp, on='site', how='outer')
df.to_csv(sys.argv[N-1], sep='\t', header=True, index=False)
