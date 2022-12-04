from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.SeqIO.FastaIO import SimpleFastaParser
from multiprocessing import Pool

# from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
# from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GenomeApi

# from ncbi.datasets.package import dataset

import tempfile

from taxonomy import get_species_under_taxon_subtree


def genome_assembly_metadata_fetch_has_error(response):
    if 'assemblies' in response:
        return False
    errors = (msg for msg in response['messages'] if 'error' in msg)
    for err in errors:
        return err['error']
    return False


def get_accessions(taxons):
    genome_client = GenomeApi()
    accessions = []
    for taxon in taxons:
        descriptors = genome_client.assembly_descriptors_by_taxon(taxon)
        # print(descriptors)
        err = genome_assembly_metadata_fetch_has_error(descriptors)
        if err:
            print('ERROR when getting taxon assembly: {} ({})'.format(taxon, err['message']))
            continue
        else:
            print('OK')
        assemblies = descriptors['assemblies']
        for assembly in assemblies:
            assembly_data = assembly['assembly']
            if assembly_data.get('assembly_category') == 'representative genome':
                accession = assembly_data['assembly_accession']
                organism = assembly_data['biosample']['description']['organism']
                accessions.append(accession)
                print('Taxon {} ({}) has assembly {}'.format(organism['tax_id'], organism['organism_name'], accession))
    return accessions

specie_taxons = get_species_under_taxon_subtree('7157')
print(specie_taxons)
print('Got {} taxons'.format(len(specie_taxons)))

# accessions = get_accessions(specie_taxons)


# zipfile_name = "gene_ds.zip"

# # download the data package through the api-client.
# with DatasetsApiClient() as api_client:
#     gene_api = DatasetsGeneApi(api_client)
#     try:
#         gene_dataset_download = gene_api.(
#             gene_ids,j
#             include_annotation_type=["FASTA_GENE", "FASTA_PROTEIN"],
#             _preload_content=False,
#         )

#         with open(zipfile_name, "wb") as f:
#             f.write(gene_dataset_download.data)
#     except DatasetsApiException as e:
#         sys.exit(f"Exception when calling GeneApi: {e}\n")

# # open the package zip archive so we can retrieve files from it
# package = dataset.GeneDataset(zipfile_name)
# # print the names and types of all files in the downloaded zip file
# print(package.get_catalog())

# # Use file types or names from the catalog to retrieve contents from specific files, e.g. protein fasta
# for protein_fasta, file_name in package.get_files_by_type("PROTEIN_FASTA"):
#     print(file_name, protein_fasta[:100])

# # get the data report and print the id and symbol for each downloaded gene
# for report in package.get_data_reports():
#     print(f"{report.gene_id}\t{report.symbol}")


def blast_protein(arg):
    """
    Takes a tuple of (protein_id, protein_seq, specie2_id) and finds the most
    similar protein in specie2, givinig the identity score between the two.
    protein_id - NCBI ref id of the protein in specie1
    protein_seq - the amino acid sequence of the protein in specie1
    specie2_id - the accession number of the specie2 genome assembly that will
                 be used to search for a correspondent protein
    Returns the tuple (protein_id, matched_protein_id, identity_score)
    """
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
            return id, None, 0.0
        matched_protein_id, identity = result.split(',')[1:3]
        return id, matched_protein_id, identity

proteins = [
    'GCF_000005575.2_reduced', # Anopheles gambiae
    'GCF_016920715.1' # Anopheles arabiensis
    ]


def get_gene_correspondance_table(specie1, specie2):
    correspondances = []
    specie1_prot_fasta = specie1 + '.faa'
    with open(specie1_prot_fasta) as handle:
        proteins = []
        for header, protein in SimpleFastaParser(handle):
            id = header.split(' ')[0]
            proteins.append((id, protein, specie2))
        with Pool() as p:
            results = p.map(blast_protein, proteins)
        for id, matched_id, identity in results:
            correspondances.append((id, matched_id, identity))
            print(id, matched_id, identity)
    avg_identity = sum(float(c[2]) for c in correspondances)/len(correspondances)
    print(avg_identity)

get_gene_correspondance_table(proteins[0], proteins[1])
