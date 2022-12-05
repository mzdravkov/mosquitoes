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


def get_accessions(taxons):
    """
    Takes an iterable of specie taxon ids and returns a list of the
    accession ids of their genome assemblies.
    """
    genome_client = GenomeApi()
    accessions = []
    for taxon, taxon_name, taxon_rank in taxons:
        descriptors = genome_client.assembly_descriptors_by_taxon(taxon)
        err = genome_assembly_metadata_fetch_has_error(descriptors)
        if err:
            logging.error('ERROR when getting taxon assembly: {} ({})'.format(taxon, err['message']))
            continue
        assemblies = descriptors['assemblies']
        for assembly in assemblies:
            assembly_data = assembly['assembly']
            if assembly_data.get('assembly_category') == 'representative genome':
                accession = assembly_data['assembly_accession']
                organism = assembly_data['biosample']['description']['organism']
                accessions.append(accession)
                logging.info('Taxon {} ({}) has assembly {}'.format(organism['tax_id'], organism['organism_name'], accession))
    return accessions