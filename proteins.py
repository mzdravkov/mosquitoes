import logging
import subprocess
import sys
import tempfile
import random

from multiprocessing import Pool
from os.path import join

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.SearchIO.BlatIO import BlatPslParser
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


def blast_protein_batch(args):
    """
    """
    prot_fasta_name, specie2 = args
    cmd = NcbiblastpCommandline(query=prot_fasta_name, db=specie2 + '.faa',
                                remote=False, num_alignments=1,
                                word_size=6,
                                matrix='BLOSUM62', evalue=10,
                                gapopen=11, gapextend=1,
                                outfmt='10')

    raw_result = cmd()[0]
    result = [line.split(',') for line in raw_result.split('\n')]
    return result


def get_protein_correspondence_table_blast(specie1, specie2, subset=False):
    """
    Takes the acession ids of two species and returns a
    protein correspondence table and an average identity score.
    The correspondence table contains records of the following type:
    (protein_id, matched_protein_id, identity_score)
    """
    correspondences = []
    specie1_prot_fasta = specie1 + '.faa'
    with open(specie1_prot_fasta) as handle:
        proteins = []
        for header, protein in SimpleFastaParser(handle):
            id = header.split(' ')[0]
            proteins.append((id, protein, specie2))
        if subset:
            proteins = random.sample(proteins, subset)
        with Pool() as p:
            results = p.map(blast_protein, proteins)
        for id, matched_id, identity in results:
            correspondences.append((id, matched_id, identity))
    avg_identity = sum(float(c[2]) for c in correspondences)/len(correspondences)
    return correspondences, avg_identity


def get_protein_correspondence_table_batched_blast(specie1, specie2, processes=16, subset=False):
    """
    Takes the acession ids of two species and returns a
    protein correspondence table and an average identity score.
    The correspondence table contains records of the following type:
    (protein_id, matched_protein_id, identity_score)
    """
    specie1_prot_fasta = specie1 + '.faa'
    
    part_files = [tempfile.NamedTemporaryFile() for i in range(processes)]

    try:
        with open(specie1_prot_fasta) as handle:
            proteins = []
            i = 0
            for header, protein in SimpleFastaParser(handle):
                part_files[i].write(bytes(header, 'utf-8'))
                part_files[i].write(bytes(protein, 'utf-8'))
                i += 1
                if i == processes:
                    i = 0
        
        for file in part_files:
            file.flush()

        with Pool() as p:
            results = p.map(blast_protein_batch, [(f.name, specie2) for f in part_files])

        correspondences = []
        total_sum = 0
        for batch in results:
            for line in batch:
                correspondences.append(line)
                total_sum += float(line[2])
        avg_identity = total_sum/len(correspondences)
        return correspondences, avg_identity
    finally:
        for f in part_files:
            f.close()


def blat_proteins(specie1, specie2, reverse=False):
    specie1_file = specie1 + '.faa'
    specie2_file = specie2 + '.faa'
    
    result_file = specie1 + '_' + specie2 + '.psl'

    # execute blat as a subprocess
    cmd = ['blat', '-prot', specie2_file, specie1_file, result_file]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # parse the output of blat and get a correspondence_table
    with open(result_file) as f:
        parser = BlatPslParser(f)
        correspondences = {}
        total_sum = 0
        for match in parser:
            ordered_hits = sorted(match.hsps, key=lambda hsp: hsp.ident_pct, reverse=True)
            best_match = None
            if len(ordered_hits) > 0:
                best_match = ordered_hits[0]
                # correspondences.append((best_match.query_id, best_match.hit_id, best_match.ident_pct))
                key = (best_match.query_id, best_match.hit_id)
                if reverse:
                    key = (best_match.hit_id, best_match.query_id)
                correspondences[key] = best_match.ident_pct
    return correspondences 


def add_missing_proteins_to_correspondences(correspondences, specie1, specie2):
    """
    Adds proteins that are not already present in the correspondences table.
    """
    all_specie1_proteins = {row[0] for row in correspondences}
    with open(specie1 + '.faa') as handle:
        for header, _ in SimpleFastaParser(handle):
            id = header.split(' ')[0]
            if id not in all_specie1_proteins:
                correspondences.append((id, None, 0))

    all_specie2_proteins = {row[1] for row in correspondences}
    with open(specie2 + '.faa') as handle:
        for header, _ in SimpleFastaParser(handle):
            id = header.split(' ')[0]
            if id not in all_specie2_proteins:
                correspondences.append((None, id, 0))


def get_protein_correspondence_table(specie1, specie2):
    """
    Takes the acession ids of two species and returns a
    protein correspondence table and an average identity score.
    The correspondence table contains records of the following type:
    (specie1_protein_id, specie2_protein_id, identity_score)
    """
    # Search specie1 proteins in specie2's proteome
    correspondences_forw = blat_proteins(specie1, specie2)
    # Search specie2 proteins in specie1's proteome
    correspondences_back = blat_proteins(specie2, specie1, reverse=True)
    
    correspondences = []
    
    total_sum = 0
    # Merge the two results into a single table
    all_pairs = set(correspondences_forw.keys()).union(correspondences_back.keys())
    for pair in all_pairs:
        if pair in correspondences_forw and pair in correspondences_back:
            # Should be the same in both places, but just in case
            identity = (correspondences_forw[pair] + correspondences_back[pair])/2
            correspondences.append((pair[0], pair[1], identity))
            total_sum += identity
        elif pair in correspondences_back:
            correspondences.append((pair[1], pair[0], correspondences_back[pair]))
            total_sum += correspondences_back[pair]
        else:
            correspondences.append((pair[0], pair[1], correspondences_forw[pair]))
            total_sum += correspondences_forw[pair]
            
    # if there are proteins in specie1 that don't have a match in specie2 or vice versa,
    # add them to the table with a score of 0.
    add_missing_proteins_to_correspondences(correspondences, specie1, specie2)
            
    avg_identity = total_sum/len(correspondences)
    
    return correspondences, avg_identity


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