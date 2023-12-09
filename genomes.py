import logging

from ncbi.datasets import GenomeApi

def genome_assembly_metadata_fetch_has_error(response):
    """
    Returns the error found in the response or False if there's no error.
    """
    if 'assemblies' not in response:
        errors = (msg for msg in response['messages'] if 'error' in msg)
        for err in errors:
            return err['error']
    return False


def get_best_available_assembly(assemblies):
    """
    Takes a list of genome assemblies and returns the best one available.
    """
    criteria = {
        'assembly_category': 'reference genome',
        'assembly_category': 'representative genome',
        'assembly_level': 'Chromosome',
    }

    for criterion, value in criteria.items():
        for assembly in assemblies:
            if assembly['assembly'].get(criterion) == value:
                return assembly
    return None


def get_accessions(taxons):
    """
    Takes an iterable of specie taxon ids and returns a dict with mappings
    between the taxons and the accession ids of their genome assemblies.
    """
    genome_client = GenomeApi()
    accessions = {}
    for taxon, _, _ in taxons:
        descriptors = genome_client.assembly_descriptors_by_taxon(str(taxon))
        err = genome_assembly_metadata_fetch_has_error(descriptors)
        if err:
            logging.error('ERROR when getting taxon assembly: {} ({})'.format(taxon, err['message']))
            continue
        assemblies = descriptors['assemblies']
        assembly = get_best_available_assembly(assemblies)
        if assembly:
            assembly_data = assembly['assembly']
            accession = assembly_data['assembly_accession']
            organism = assembly_data['biosample']['description']['organism']
            accessions[taxon] = accession
            logging.info('Taxon {} ({}) has assembly {}'.format(organism['tax_id'], organism['organism_name'], accession))
    return accessions