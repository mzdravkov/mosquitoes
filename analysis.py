from collections import defaultdict

import networkx as nx

from storage import read_correspondences


# TODO: Get anthropophily values for all species

# Score indicating the degree of anthropophily of a species
# in the range of [-1, 1], where -1 indicates a species that
# is non-anthropophilic and 1 indicates a species that is
# very anthropophilic.
SPECIES_ANTHROPOPHILY = {
'GCA_000441895.2': -0.918, # Anopheles sinensis [Ree, Han-Il, et al. 2001][PMC2712014]
'GCF_000005575.2': 0.52, # Anopheles gambiae str. PEST [takken (PMC4381365)]
'GCF_006496715.1': 0.96, # Aedes albopictus [Alongkot Ponlawat 2005]
'GCF_013141755.1': -0.982, # Anopheles stephensi [Thomas, S., Ravishankaran 2017]
'GCF_013758885.1': -0.9, # Anopheles albimanus [Bruce-Chwatt 1966 (PMC2476083)]
'GCF_015732765.1': -0.33, # Culex quinquefasciatus [takken (PMC4381365)]
'GCF_016801865.1': -0.286, # Culex pipiens pallens [Joaquín Muñoz 2011]
'GCF_016920715.1': -0.09, # Anopheles arabiensis [tekken (PMC4381365)]
'GCF_017562075.2': -0.76, # Anopheles merus [Pamela C Kipyab 2013]
'GCF_943734665.1': 0.0, # Anopheles aquasalis
'GCF_943734685.1': 0.96, # Anopheles coluzzii [Martin C. Akogbéto 2018]
'GCF_943734695.1': 0.0, # Anopheles maculipalpis
'GCF_943734745.1': 0.65, # Anopheles darlingi [Marta Moreno 2017]
'GCF_943734845.2': 0.21, # Anopheles funestus [tekken (PMC4381365)]
}


def graph(correspondence_arrs):
    edges = set()
    index = {}
    reverse_index = {}
    i = 0
    for pair, correspondences in correspondence_arrs.items():
        for c in correspondences:
            if c[0] not in index:
                index[c[0]] = i
                reverse_index[i] = c[0]
                i += 1
            if c[1] not in index:
                index[c[1]] = i
                reverse_index[i] = c[1]
                i += 1
            e1 = index[c[0]]
            e2 = index[c[1]]
            edge = tuple(sorted([e1, e2]))
            # if edge not in edges:
            #     weight = 1.0 / (float(c[2]) + 1)
            #     edges[edge] = weight
            edges.add(edge)
        print('.', end='', flush=True)

    print('\n')

    g = nx.Graph()
    g.add_nodes_from(index.values())
    g.add_edges_from(edges)
    g.remove_node(index[''])
    # g.add_weighted_edges_from((e[0], e[1], w) for e, w in edges.items())
    return g, index, reverse_index


def get_correspondences_graph(correspondences):
    g, index, rindex = graph(correspondences)
    return g, index, rindex


def get_neighbors(g, index, rindex, gene, depth=1):
    vertices = {gene}
    for i in range(depth):
        new = set()
        for v in vertices:
            neighbors = [rindex[n] for n in nx.neighbors(g, index[gene]) if rindex[n] != '']
            new.update(neighbors)
        size_old = len(vertices)
        vertices.update(new)
        size_new = len(vertices)
        if size_old == size_new:
            return vertices, i
    return vertices, depth


def get_communities(graph, rindex):
    communities = nx.algorithms.community.louvain_communities(graph, resolution=100)
    # map ids to proteins accessions
    return [{rindex[p] for p in c} for c in communities]


def print_histogram(depths):
    for depth in range(1, max(depths.keys()) + 1):
        print(f'{depth}: {depths[depth]}')


def get_community_membership_index(communities):
    community_index = {}
    for i in range(len(communities)):
        for gene in communities[i]:
            community_index[gene] = i
    return community_index


def calculate_community_relevance(genes_similarity, species_similarity, anthropophily1, anthropophily2):
    # similar anthropophily:
    # rank = score(gene1, gene2)/score(specie1, specie2)
    # dissimilar anthropophily:
    # rank = score(specie1, specie2)/score(gene1, gene2)
    # => rank = (score(gene1, gene2)/score(specie2, specie1))^similarity(anthropophily1, anthropophily2)

    # TODO: think about finding a better way to measure anthropophily similarity
    anthropophily_similarity = 1 - abs(anthropophily1 - anthropophily2)

    return (genes_similarity/(species_similarity + 1)) ** anthropophily_similarity


def histogram_buckets(values, num_buckets):
    min_value = min(values)
    interval = max(values) - min_value
    step = interval / num_buckets

    hist = defaultdict(lambda: 0)
    for v in values:
        for i in range(num_buckets):
            bucket = i * step + min_value
            if v <= bucket:
                hist[bucket] += 1
                break
    return hist


if __name__ == '__main__':
    correspondences = read_correspondences()
    print('getting correspondence graph')
    g, index, rindex = get_correspondences_graph(correspondences)

    print('nodes in the graph', len(g.nodes))

    # groups = {}
    # depths = defaultdict(lambda: 0)
    # print('\ngetting groups')
    # i = 0
    # for gene in index:
    #     if gene == '':
    #         continue
    #     if i % 1000 == 0:
    #         print('.', end='', flush=True)
    #     groups[gene], depth = get_neighbors(g, index, rindex, gene, depth=10)
    #     depths[depth] += 1
    #     if depth == 10:
    #         print(f'{gene} -> {depth}')
    #     i += 1
    # print('\ndepths')
    # print_histogram(depths)
    # print('\ngroup sizes')
    # group_sizes = defaultdict(lambda: 0)
    # for key, group in groups.items():
    #     group_sizes[len(group)] += 1
    # print_histogram(group_sizes)
    # print('total members of groups', sum(group_sizes.values()))
    # print('number of groups', len(groups))
    # print('average group size', sum(len(g) for g in groups.values())/len(groups))
    # print('max group size', max(len(g) for g in groups.values()))

    # communities = nx.algorithms.community.louvain_communities(g, resolution=100)
    print('getting communities')

    communities = get_communities(g, rindex)

    print('number of communities', len(communities))
    print('average size of community', sum(len(c) for c in communities)/len(communities))
    print('max community size', max(len(c) for c in communities))

    com_hist = defaultdict(lambda: 0)
    for c in communities:
        com_hist[len(c)] += 1
    print_histogram(com_hist)

    community_index = get_community_membership_index(communities)

    community_relevancies = defaultdict(lambda: {})

    print('processing correspondences to find community relevance scores')

    for (genome1, genome2), correspondences in correspondences.items():
        species1_anthropophily = SPECIES_ANTHROPOPHILY[genome1]
        species2_anthropophily = SPECIES_ANTHROPOPHILY[genome2]

        print('.', end='', flush=True)
        total_species_score = sum(float(score) for _, _, score in correspondences)
        avg_species_score = total_species_score / len(correspondences)

        for correspondence in correspondences:
            gene1, gene2, score_str = correspondence
            score = float(score_str)
            score = 0.0 if score < 0 else score

            relevance = calculate_community_relevance(
                score,
                avg_species_score,
                species1_anthropophily,
                species2_anthropophily
            )

            if gene1 != '':
                community1 = community_index[gene1]
                community_relevancies[community1][(genome1, genome2)] = relevance

            if gene2 != '':
                community2 = community_index[gene2]
                if community1 != community2:
                    community_relevancies[community2][(genome1, genome2)] = relevance

    community_scores = {}
    for community, relevancies in community_relevancies.items():
        community_scores[community] = sum(relevancies.values())/len(relevancies)

    print('\nHistogram of relevance buckets')

    print(f'min: ', min(community_scores.values()))
    buckets = histogram_buckets(community_scores.values(), 10)
    for bucket, count in sorted(buckets.items(), key=lambda x: x[0]):
        print(f'{bucket}: {count}')

    print('Top 20 communities:')
    for community, score in sorted(community_scores.items(), key=lambda x: x[1], reverse=True)[:20]:
        print(f'{community}: {score}')
        print(f'\t{communities[community]}')


