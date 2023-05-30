import pandas as pd
import sys

PRE=sys.argv[1]
MinCPGCount = int(sys.argv[2])
MinDepth = int(sys.argv[3])
INPUT = sys.argv[4]

df = pd.read_table(INPUT, names=['BIN', 'CPG', 'READ', 'MET'])

BinCPGCount = df[['BIN','CPG']].drop_duplicates('CPG').groupby('BIN').count()

CPGCountFiltBin = BinCPGCount.loc[BinCPGCount.CPG>=MinCPGCount] # Filt by >=3 CPGs / 200 bp

ReadCPGCount = df[['BIN','CPG']].groupby('BIN').count()

BinCPGReadCov = pd.merge(CPGCountFiltBin, ReadCPGCount, left_index=True, right_index=True)

CovFilt = BinCPGReadCov.loc[BinCPGReadCov.CPG_y/BinCPGReadCov.CPG_x>=MinDepth].index.to_list()

CPGCtCovFilt = list(set(CovFilt).intersection(set(CPGCountFiltBin.index)))

CPGCtCovFiltDF = df.loc[df.BIN.isin(CPGCtCovFilt)]

AggreRatio= (CPGCtCovFiltDF[['BIN', 'MET']].groupby('BIN').sum() / CPGCtCovFiltDF[['BIN', 'MET']].groupby('BIN').count()).reset_index().rename(columns={'MET':'RATIO'})

BinReadCPGCt = CPGCtCovFiltDF[['BIN', 'READ', 'CPG']].groupby(['BIN','READ']).count().rename(columns={'CPG':'CPGCount'})

BinReadCPGRatio = (CPGCtCovFiltDF[['BIN', 'READ', 'MET']].groupby(['BIN','READ']).sum()/CPGCtCovFiltDF[['BIN', 'READ', 'MET']].groupby(['BIN','READ']).count()).rename(columns={'MET':'READ_RATIO'})

BinReadCPGCtRatio = pd.concat([BinReadCPGCt, BinReadCPGRatio],axis=1)

BinReadCPGCtRatioFilt = BinReadCPGCtRatio.loc[BinReadCPGCtRatio.CPGCount>=MinCPGCount]

BinReadRatios = BinReadCPGCtRatioFilt.reset_index().merge(AggreRatio, on='BIN')

BinReadRatios['ReadBinRatioDiff'] = BinReadRatios['READ_RATIO'] - BinReadRatios['RATIO']

BinDiffReadCt = BinReadRatios[['BIN', 'ReadBinRatioDiff']].loc[BinReadRatios['ReadBinRatioDiff'].abs()>=0.2 ].groupby('BIN').count().rename(columns={'ReadBinRatioDiff':'DiffReadCt'})

BinDiffReadCtFilt = BinDiffReadCt.loc[BinDiffReadCt.DiffReadCt>=4]

BinReadCt = BinReadRatios[['BIN', 'READ']].groupby('BIN').count().rename(columns={'READ':"ReadCt"}).reset_index()

BinReadDiffCt =  BinDiffReadCtFilt.reset_index().merge(BinReadCt,on='BIN')

# Bins with >=80% reads clustered as H or L
Bincluster = BinReadDiffCt.loc[BinReadDiffCt.DiffReadCt / BinReadDiffCt.ReadCt >= 0.8] 


# Statics for bins with clusters
BinReadRatiosClu = BinReadRatios.loc[(BinReadRatios.BIN.isin(Bincluster.BIN))][['BIN','ReadBinRatioDiff']]
BinReadRatiosClu.loc[BinReadRatiosClu.ReadBinRatioDiff>=0.2, 'Cluster'] = "H"
BinReadRatiosClu.loc[BinReadRatiosClu.ReadBinRatioDiff<=-0.2, 'Cluster'] = "L"
BinReadRatiosClu = BinReadRatiosClu.fillna('U')
CluReadSta = BinReadRatiosClu.groupby(['BIN','Cluster']).count().reset_index().pivot(index='BIN', columns='Cluster', values='ReadBinRatioDiff').fillna(0).reset_index()
CluReadSta = CluReadSta.merge(BinCPGCount.reset_index(), on='BIN')
CluReadSta = CluReadSta.merge(AggreRatio, on='BIN')
CluReadSta['H']=CluReadSta['H'].astype(int)
CluReadSta['L']=CluReadSta['L'].astype(int)
CluReadSta['U']=CluReadSta['U'].astype(int)
CluReadSta.to_csv(f'{PRE}.Bin_CluReadSta.tsv',sep="\t", header=True, index=False)

# Filt unclustered reads for each cluster_bin and assign H and L tag to each bin read
BinReadRatiosFilt = BinReadRatios.loc[(BinReadRatios.BIN.isin(Bincluster.BIN)) & ((BinReadRatios.ReadBinRatioDiff>=0.2) | (BinReadRatios.ReadBinRatioDiff<=-0.2))]
BinReadRatiosFilt.loc[BinReadRatiosFilt.ReadBinRatioDiff>=0.2, 'Cluster'] = "H"
BinReadRatiosFilt.loc[BinReadRatiosFilt.ReadBinRatioDiff<=-0.2, 'Cluster'] = "L"


CluCount = BinReadRatiosFilt[['READ','Cluster', 'ReadBinRatioDiff']].groupby(['READ','Cluster']).count().rename(columns={'ReadBinRatioDiff':'CluCount'}).reset_index()
LowMet = CluCount.loc[CluCount.Cluster=="L"][['READ', 'CluCount']].rename(columns={'CluCount':'LowMet'})
HighMet = CluCount.loc[CluCount.Cluster=="H"][['READ', 'CluCount']].rename(columns={'CluCount':'HighMet'})
ReadMet = LowMet.merge(HighMet, how='outer').fillna(0)
ReadMet['LowMet'] = ReadMet['LowMet'].astype(int)
ReadMet['HighMet'] = ReadMet['HighMet'].astype(int)

CluCPGCount = BinReadRatiosFilt[['READ', 'CPGCount', 'Cluster']].groupby(['Cluster', 'READ']).sum()

HCluCPGCount = CluCPGCount.loc['H'].rename(columns={'CPGCount':'H_CPGCount'})
LCluCPGCount = CluCPGCount.loc['L'].rename(columns={'CPGCount':'L_CPGCount'})

MCluCPGCount = HCluCPGCount.merge(LCluCPGCount, left_index=True, right_index=True, how='outer').fillna(0)
McluCPGCount = MCluCPGCount.reset_index().merge(ReadMet, on='READ')

McluCPGCount.loc[McluCPGCount.H_CPGCount > McluCPGCount.L_CPGCount, 'TAG'] = 'H'
McluCPGCount.loc[McluCPGCount.H_CPGCount < McluCPGCount.L_CPGCount, 'TAG'] = 'L'
McluCPGCount.loc[McluCPGCount.H_CPGCount < McluCPGCount.L_CPGCount, 'TAG'] = 'U'

McluCPGCount['H_CPGCount'] = McluCPGCount['H_CPGCount'].astype(int)
McluCPGCount['L_CPGCount'] = McluCPGCount['L_CPGCount'].astype(int)

McluCPGCount.to_csv(f'{PRE}.McluCPGCount.tsv', header=True, sep="\t", index=False)

