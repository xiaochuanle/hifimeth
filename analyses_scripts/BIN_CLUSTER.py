import pandas as pd
import sys
import argparse
from scipy.stats import fisher_exact

# Create a parser object
parser = argparse.ArgumentParser(description="Python script for clustering reads into HMG and LMG and identifying clustered bins")

# Add arguments to the parser
parser.add_argument("-P", "--prefix", type=str, default="SAMPLE", help="Prefix for the output files [default: SAMPLE]")
parser.add_argument("-c", "--CPG", type=int, default=5, help="Minimum number of CpGs in a bin for it to be considered [default: 5]")
parser.add_argument("-r", "--readCPG", type=int, default=5, help="Minimum number of CpGs in a read intersecting a bin for the read to be analyzed [default: 5]")
parser.add_argument("-d", "--depth", type=int, default=8, help="Minimum number of reads harboring enough CpGs (defined by -r) in a bin for the bin to be considered [default: 8] ")
parser.add_argument("-i", "--input", type=str, help="Input file containing four columns: bin_name CpG_coordinate read_id CpG_methylation_state")
parser.add_argument("-D", "--diff", type=float, default=0.20, help="Minimum methylation ratio diff to the bin methylation frequency for a read to be clustered to HMG or LMG [default: 0.20]")
parser.add_argument("-R", "--cluster_ratio", type=float, default=0.80, help="Minimum ratio of HMG&LMG grouped reads for identifying a clustered bin [default: 0.80]")

# Parse the arguments
args = parser.parse_args()

PRE = args.prefix
MinCPGBin = args.CPG
MinCPGRead = args.readCPG
MinDepth = args.depth
INPUT = args.input
MinRatioDiff = args.diff
ClusterRatio = args.cluster_ratio

def Fisher(df_row):
    return fisher_exact([[df_row.HMG_M, df_row.HMG_U],[df_row.LMG_M, df_row.LMG_U]])[-1]

df = pd.read_table(INPUT, names=['BIN', 'CPG', 'READ', 'MET'])

# Filter bins via minimum bin CpG count requirement
BinCPGCount = df[['BIN','CPG']].drop_duplicates('CPG').groupby('BIN').count()
CPGCountFiltBin = BinCPGCount.loc[BinCPGCount.CPG>=MinCPGBin] 

#
ReadCPGCount = df[['BIN','CPG']].groupby('BIN').count()
BinCPGReadCov = pd.merge(CPGCountFiltBin, ReadCPGCount, left_index=True, right_index=True)

CovFilt = BinCPGReadCov.loc[BinCPGReadCov.CPG_y/BinCPGReadCov.CPG_x>=MinDepth].index.to_list()

CPGCtCovFilt = list(set(CovFilt).intersection(set(CPGCountFiltBin.index)))

CPGCtCovFiltDF = df.loc[df.BIN.isin(CPGCtCovFilt)]

AggreRatio= (CPGCtCovFiltDF[['BIN', 'MET']].groupby('BIN').sum() / CPGCtCovFiltDF[['BIN', 'MET']].groupby('BIN').count()).reset_index().rename(columns={'MET':'RATIO'})

BinReadCPGCt = CPGCtCovFiltDF[['BIN', 'READ', 'CPG']].groupby(['BIN','READ']).count().rename(columns={'CPG':'CPGCount'})

BinReadCPGRatio = (CPGCtCovFiltDF[['BIN', 'READ', 'MET']].groupby(['BIN','READ']).sum()/CPGCtCovFiltDF[['BIN', 'READ', 'MET']].groupby(['BIN','READ']).count()).rename(columns={'MET':'READ_RATIO'})

BinReadCPGCtRatio = pd.concat([BinReadCPGCt, BinReadCPGRatio],axis=1)

BinReadCPGCtRatioFilt = BinReadCPGCtRatio.loc[BinReadCPGCtRatio.CPGCount>=MinCPGRead]

BinReadRatios = BinReadCPGCtRatioFilt.reset_index().merge(AggreRatio, on='BIN')

BinReadRatios['ReadBinRatioDiff'] = BinReadRatios['READ_RATIO'] - BinReadRatios['RATIO']

BinDiffReadCt = BinReadRatios[['BIN', 'ReadBinRatioDiff']].loc[BinReadRatios['ReadBinRatioDiff'].abs()>=MinRatioDiff ].groupby('BIN').count().rename(columns={'ReadBinRatioDiff':'DiffReadCt'})

HMG = BinReadRatios.loc[BinReadRatios.ReadBinRatioDiff>=MinRatioDiff]
HMG['HMG_M'] = round(HMG.CPGCount * HMG.READ_RATIO, 0)
HMG['HMG_M'] = HMG.HMG_M.astype(int)
HMG['HMG_U'] = HMG.CPGCount - HMG.HMG_M
LMG = BinReadRatios.loc[BinReadRatios.ReadBinRatioDiff<=-MinRatioDiff]
LMG['LMG_M'] = round(LMG.CPGCount * LMG.READ_RATIO, 0)
LMG['LMG_M'] = LMG.LMG_M.astype(int)
LMG['LMG_U'] = LMG.CPGCount - LMG.LMG_M

BinDiffReadCtFilt = BinDiffReadCt.loc[BinDiffReadCt.DiffReadCt>=MinDepth]

BinReadCt = BinReadRatios[['BIN', 'READ']].groupby('BIN').count().rename(columns={'READ':"ReadCt"}).reset_index()

BinReadDiffCt =  BinDiffReadCtFilt.reset_index().merge(BinReadCt,on='BIN')

# Bins with >=80% reads clustered as H or L
Bincluster = BinReadDiffCt.loc[BinReadDiffCt.DiffReadCt / BinReadDiffCt.ReadCt >= ClusterRatio] 

HMG_SUM = HMG.loc[HMG.BIN.isin(Bincluster.BIN)][['BIN', 'HMG_M', 'HMG_U']].groupby('BIN').sum()
LMG_SUM = LMG.loc[LMG.BIN.isin(Bincluster.BIN)][['BIN', 'LMG_M', 'LMG_U']].groupby('BIN').sum()
FisherTest = pd.merge(HMG_SUM, LMG_SUM, left_index=True, right_index=True)
FisherTest['Fisher_pvalue'] = FisherTest.apply(Fisher, axis=1)
FisherTest = FisherTest.loc[FisherTest['Fisher_pvalue']<0.001].copy()
Bincluster = Bincluster.loc[Bincluster.BIN.isin(FisherTest.index)]

# Statics for bins with clusters
BinReadRatiosClu = BinReadRatios.loc[(BinReadRatios.BIN.isin(Bincluster.BIN))][['BIN','ReadBinRatioDiff']]
BinReadRatiosClu.loc[BinReadRatiosClu.ReadBinRatioDiff>=MinRatioDiff, 'Cluster'] = "H"
BinReadRatiosClu.loc[BinReadRatiosClu.ReadBinRatioDiff<=-MinRatioDiff, 'Cluster'] = "L"
BinReadRatiosClu = BinReadRatiosClu.fillna('U')
CluReadSta = BinReadRatiosClu.groupby(['BIN','Cluster']).count().reset_index().pivot(index='BIN', columns='Cluster', values='ReadBinRatioDiff').fillna(0).reset_index()
CluReadSta = CluReadSta.merge(BinCPGCount.reset_index(), on='BIN')
CluReadSta = CluReadSta.merge(AggreRatio, on='BIN')
CluReadSta['H']=CluReadSta['H'].astype(int)
CluReadSta['L']=CluReadSta['L'].astype(int)
CluReadSta['U']=CluReadSta['U'].astype(int)
CluReadSta = CluReadSta.merge(FisherTest.reset_index(), on='BIN').copy()
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
