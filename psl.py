import pandas as pd

from Bio.SearchIO.BlatIO import _calc_millibad


PSL_HEADER = [
    "matches", # matches
    "mismatches", # misMatches
    "repmatches", # repMatches
    "ncount", # nCount
    "qnuminsert", # qNumInsert
    "qbaseinsert", # qBaseInsert
    "tnuminsert", # tNumInsert
    "tbaseinsert", # tBaseInsert
    "strand", # strand
    "qname", # qName
    "qsize", # qSize
    "qstart", # qStart
    "qend", # qEnd
    "tname", # tName
    "tsize", # tSize
    "tstart", # tStart
    "tend", # tEnd
    "blockcount", # blockCount
    "blocksizes", # blockSizes
    "qstarts", # qStarts
    "tstarts", # tStarts
]


def calc_identity_pct(row, is_protein=True):
    return 100.0 - _calc_millibad(row, is_protein)*0.1


def _get_identity_column(psl, is_protein=True):
    return psl.apply(calc_identity_pct, axis=1, is_protein=is_protein)


def parse_psl(psl_path, is_protein=True):
    df = pd.read_csv(psl_path, sep='\t', skiprows=5, names=PSL_HEADER)
    df.insert(21, "identity", _get_identity_column(df))

    # get only the best match for each query protein
    df_max = df.groupby('qname').apply(lambda x: x.loc[x['identity'].idxmax()])

    return df_max