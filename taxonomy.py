from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.api.taxonomy_api import TaxonomyApi

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


def get_species_under_taxon_subtree(taxon):
    """
    get_species_under_taxon_subtree gets a taxon identifier (as string), e.g. '7157' and returns
    a set with the taxonomy identifiers of its subtree leafs (i.e. species under the taxon).
    Example:
    get_species_under_taxon_subtree('7157') will return a list of all species in the Culicidae family.
    """
    with ApiClient() as api_client:
        taxonomy_client = TaxonomyApi(api_client)
        try:
            response = taxonomy_client.taxonomy_filtered_subtree([taxon])
            print(response)
            return get_taxonomy_leafs(response)
        except Exception as e:
            print(f"Exception when calling TaxonomyApi: {e}")