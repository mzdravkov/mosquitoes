import pandas as pd

from Bio.SearchIO.BlatIO import _calc_score


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


def _get_score_column(psl, is_protein=False):
    return psl.apply(_calc_score, axis=1, is_protein=is_protein)


def parse_psl(psl_path):
    df = pd.read_csv(psl_path, sep='\t', skiprows=5, names=PSL_HEADER)
    df.insert(21, "score", _get_score_column(df))

    # get only the best match for each query protein
    df_max = df.groupby('qname').apply(lambda x: x.loc[x['score'].idxmax()])
    
    # change negative scores to 0
    df_max.loc[df_max['score'] < 0, 'score'] = 0

    return df_max