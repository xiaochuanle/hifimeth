#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys

df1 = pd.read_table(sys.argv[1], header=0)
df2 = pd.read_table(sys.argv[2], header=0)

df = df1.merge(df2,on='site').dropna()

print('\t'.join([sys.argv[1], sys.argv[2], str(np.corrcoef(df[df.columns[-1]], df[df.columns[-2]])[0, 1])]))
