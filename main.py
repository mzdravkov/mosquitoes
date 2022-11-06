from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.SeqIO.FastaIO import SimpleFastaParser
from multiprocessing import Pool
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.api.taxonomy_api import TaxonomyApi
import tempfile



def get_taxonomy_leafs(taxonomy_response):
    """
    get_taxonomy_leafs will read the response from taxonomy API taxonomy_filtered_subtree
    and extract a set of the taxonomy ids of the subtree leafs (species).
    """
    edges = taxonomy_response['edges']
    queue = []
    # Add the root ('1') to start traversing the tree
    queue.append('1')
    leafs = set()
    while len(queue) > 0:
        node = queue.pop()
        children = edges[node].get('visible_children')
        if not children:
            leafs.add(node)
            continue
        for child in children:
            queue.append(str(child))
    return leafs

with ApiClient() as api_client:
    taxonomy_client = TaxonomyApi(api_client)
    try:
        response = taxonomy_client.taxonomy_filtered_subtree(['7164'])
        print(get_taxonomy_leafs(response))
    except Exception as e:
        print(f"Exception when calling TaxonomyApi: {e}\n")


def blast_protein(arg):
    id, protein, specie2 = arg
    with tempfile.NamedTemporaryFile() as f:
        f.write(bytes(protein, 'utf-8'))
        f.flush()
        cline = NcbiblastpCommandline(query=f.name, db=specie2 + '.faa',
                                    remote=False, num_alignments=1,
                                    word_size=6,
                                    matrix='BLOSUM62', evalue=10,
                                    gapopen=11, gapextend=1,
                                    outfmt='10')

        result = cline()[0]
        if not result:
            return id, None, None
        matched_protein_id, identity = result.split(',')[1:3]
        return id, matched_protein_id, identity

proteins = [
    'GCF_000005575.2_reduced', # Anopheles gambiae
    'GCF_016920715.1' # Anopheles arabiensis
    ]


def get_gene_correspondance_table(specie1, specie2):
    specie1_prot_fasta = specie1 + '.faa'
    with open(specie1_prot_fasta) as handle:
        proteins = []
        for header, protein in SimpleFastaParser(handle):
            id = header.split(' ')[0]
            proteins.append((id, protein, specie2))
            # matched_protein_id, identity = blast_protein(protein, specie2)
            # print(id, matched_protein_id, identity)
            # # break
        with Pool() as p:
            results = p.map(blast_protein, proteins)
        for id, matched_id, identity in results:
            print(id, matched_id, identity)

get_gene_correspondance_table(proteins[0], proteins[1])
