from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.api.taxonomy_api import TaxonomyApi
from ncbi.datasets import GenomeApi


TAXON_RANKS_UNDER_SPECIES = {
    'SPECIES',
    'SUBSPECIES',
    'FORMA',
    'VARIETAS',
    'STRAIN',
    'SUBVARIETY'
}


def filter_species_and_subranks(taxonomy_response):
    """
    filter_species_and_subranks will read the response from taxonomy API genome_tax_tree
    and extract information about the species (and other ranks under species) in the subtree.
    """
    # Add the root node
    queue = [taxonomy_response]
    leafs = set()
    while len(queue) > 0:
        node = queue.pop()
        if node['rank'].to_str() in TAXON_RANKS_UNDER_SPECIES:
            leaf = (node['tax_id'], node['title'], node['rank'].to_str())
            leafs.add(leaf)
        children = node.get('children')
        if children:
            for child in children:
                queue.append(child)
    return leafs


def get_species_in_taxon_subtree(taxon):
    """
    get_species_in_taxon_subtree gets a taxon identifier (as string), e.g. '7157' and returns
    a list of the species (and other ranks below species) in the subtree under the given taxon.
    Example:
    get_species_in_taxon_subtree('7157') will return a list of all species in the Culicidae family.
    Return value type: [(taxon_id, taxon_name, taxon_rank), ...]
    """
    genome_client = GenomeApi()
    response = genome_client.genome_tax_tree(taxon)
    return filter_species_and_subranks(response)