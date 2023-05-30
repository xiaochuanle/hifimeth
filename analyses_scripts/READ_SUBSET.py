#!/usr/bin/env python3

import pandas as pd
import sys

df1 = pd.read_table(sys.argv[1], names=['READ', 'col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7'])
df2 = pd.read_table(sys.argv[2], names=['READ'])
df1.loc[df1['READ'].isin(df2['READ'])].to_csv(f'{sys.argv[2]}.5mc-call.txt',sep="\t", index=False, header=False)
