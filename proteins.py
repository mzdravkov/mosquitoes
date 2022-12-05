import logging
import sys
import tempfile

from multiprocessing import Pool
from os.path import join

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ncbi.datasets import GenomeApi
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.package import dataset

from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest
from ncbi.datasets.openapi.model.v1_annotation_for_assembly_type import V1AnnotationForAssemblyType

from storage import SEQUENCES_DIR
from storage import get_dataset_filename


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
            logging.info('%s %s %s', id, None, 0.0)
            return id, None, 0.0
        matched_protein_id, identity = result.split(',')[1:3]
        logging.info('%s %s %s', id, matched_protein_id, identity)
        return id, matched_protein_id, identity


def get_protein_correspondance_table(specie1, specie2):
    """
    Takes the acession ids of two species and returns a
    protein correspondance table and an average identity score.
    The correspondance table contains records of the following type:
    (protein_id, matched_protein_id, identity_score)
    """
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
    avg_identity = sum(float(c[2]) for c in correspondances)/len(correspondances)
    return correspondances, avg_identity


def download_dataset(accession):
    """
    Downloads the protein fasta file for the genome assembly with the given accession number.
    """
    zipfile_name = get_dataset_filename(accession)
    
    # Check if zipfile_name already exists, skipping the download if it does.
    try:
        with open(zipfile_name, 'r') as f:
            print('exists')
            logging.info('File %s already exists. Skipping.', zipfile_name)
            return
    except FileNotFoundError:
        pass

    # download the data package through the api-client.
    with DatasetsApiClient() as api_client:
        genome_api = GenomeApi(api_client)
        try:
            gene_dataset_download = genome_api.download_assembly_package(
                [accession],
                include_annotation_type=[V1AnnotationForAssemblyType("PROT_FASTA")],
                exclude_sequence=True,
                filename=zipfile_name,
                _preload_content=False,
            )

            with open(zipfile_name, "wb") as f:
                f.write(gene_dataset_download.data)
        except DatasetsApiException as e:
            sys.exit(f"Exception when calling GeneApi: {e}\n")


def get_protein_from_dataset(accession):
    """
    Extracts the protein fasta file from the dataset with the given accession number.
    """
    zipfile_name = get_dataset_filename(accession)

    # open the package zip archive so we can retrieve files from it
    package = dataset.GeneDataset(zipfile_name)

    # TODO: check why we can have more than one file and what to do with it
    return package.get_files_by_type("PROTEIN_FASTA")[0]