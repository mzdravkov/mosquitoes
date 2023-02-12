import statistics

from collections import defaultdict
from io import StringIO
import subprocess
import tempfile

import networkx as nx
from proteins import get_protein_sequence

from storage import read_correspondences
from utils import print_table

from Bio  import AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqIO.FastaIO import SimpleFastaParser

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
'GCF_016801865.2': -0.286, # Culex pipiens pallens [Joaquín Muñoz 2011]
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


def calculate_relevance(genes_similarity, species_similarity, anthropophily1, anthropophily2):
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


def top_by_relevance_community(n):
    correspondences = read_correspondences()
    print('getting correspondence graph')
    g, index, rindex = get_correspondences_graph(correspondences)

    print('nodes in the graph', len(g.nodes))

    print('getting communities')

    communities = get_communities(g, rindex)

    print('number of communities', len(communities))
    print('average community size', sum(len(c) for c in communities)/len(communities))
    print('median community size', statistics.median(len(c) for c in communities))
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
            # negative and zero scores are breaking the relevance score calculation
            # because they may result in a complex number or division by zero, so we
            # just set them to 1
            score = 1 if score <= 0 else score

            relevance = calculate_relevance(
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
        # We calculate the relevance by using the signal-to-noise ratio.
        # Essentially, we want to find communities that have a high average
        # score, but this is not due to outliers that skew the average.
        mean = sum(relevancies.values())/len(relevancies)
        variance = sum((r - mean)**2 for r in relevancies.values())/len(relevancies)
        community_scores[community] = mean/variance

    print('\nHistogram of relevance buckets')

    print(f'min: ', min(community_scores.values()))
    buckets = histogram_buckets(community_scores.values(), 10)
    for bucket, count in sorted(buckets.items(), key=lambda x: x[0]):
        print(f'{bucket}: {count}')

    print('Top {n} communities:')
    for community, score in sorted(community_scores.items(), key=lambda x: x[1], reverse=True)[:n]:
        print(f'{community}: {score}')
        for protein in communities[community]:
            print(f'\t- {protein}')


def check_species_conflict(gene, genome, protein_species):
    if gene in protein_species and protein_species[gene] != genome:
        raise RuntimeError(
            "Gene {} in two species: {} and {}".format(
                gene, protein_species[gene], genome)
        )


def top_by_relevance(correspondences, n, test_species=None):
    protein_species = {}
    homologs = defaultdict(lambda: {})
    species = set()

    for (genome1, genome2), correspondences in correspondences.items():
        if genome1 == test_species or genome2 == test_species:
            continue

        species1_anthropophily = SPECIES_ANTHROPOPHILY[genome1]
        species2_anthropophily = SPECIES_ANTHROPOPHILY[genome2]

        print('.', end='', flush=True)
        total_species_score = sum(float(score) for _, _, score in correspondences)
        avg_species_score = total_species_score / len(correspondences)

        for correspondence in correspondences:
            gene1, gene2, score_str = correspondence

            if gene1:
                check_species_conflict(gene1, genome1, protein_species)
                protein_species[gene1] = genome1
                species.add(genome1)

                if gene1 not in homologs:
                    homologs[gene1] = {}

            if gene2:
                check_species_conflict(gene2, genome2, protein_species)
                protein_species[gene2] = genome2
                species.add(genome2)
            
                if gene2 not in homologs:
                    homologs[gene2] = {}

            if gene1 == '' or gene2 == '':
                continue

            score = float(score_str)
            # negative and zero scores are breaking the relevance score calculation
            # because they may result in a complex number or division by zero, so we
            # just set them to 1
            score = 1 if score <= 0 else score

            relevance = calculate_relevance(
                score,
                avg_species_score,
                species1_anthropophily,
                species2_anthropophily
            )
            
            homologs[gene1][gene2] = relevance
            
    print('\nnumber of proteins', len(homologs))
    print('average homologous set size', sum(len(c) for c in homologs)/len(homologs))
    print('median homologous set size', statistics.median(len(c) for c in homologs))
    print('max homologous set size', max(len(c) for c in homologs))
            
    protein_relevancies = {}
    
    for gene, relevancies in homologs.items():
        if len(relevancies) > 0:
            mean = sum(relevancies.values())/len(relevancies)
            protein_relevancies[gene] = mean
        else:
            # TODO: this is a bit random, but I'm not sure how to handle those cases
            protein_relevancies[gene] = 0.0
            
    num_species = len(species)

    # species_protein_ranking = defaultdict(lambda: {})
    # max_indices = defaultdict(lambda: 0)
    # for protein, relevance in sorted(protein_relevancies.items(), key=lambda x: x[1], reverse=True):
    #     species = protein_species[protein]
    #     species_protein_ranking[species][protein] = max_indices[species]
    #     max_indices[species] += 1

    # avg_ranks = {}
    # for protein in protein_relevancies:
    #     prot_homologs = homologs[protein]
    #     ranks = [species_protein_ranking[protein_species[h]][h] for h in prot_homologs if h in species_protein_ranking[protein_species[h]]]
    #     if len(ranks) > 0:
    #         mean = sum(ranks)/len(ranks)
    #         avg_ranks[protein] = mean
    
    # i = 0
    # data = []
    # for protein, avg_rank in sorted(avg_ranks.items(), key=lambda x: x[1]):
    #     if len(homologs[protein]) < num_species/2:
    #         continue
        
    #     if i > n:
    #         break

    #     data.append((protein, protein_relevancies[protein], avg_rank, len(homologs[protein])))
    #     # for homolog, pairwise_relevance in homologs[protein].items():
    #     #     data.append((protein, relevance, homolog, protein_relevancies.get(homolog, None), pairwise_relevance))
    #     i += 1
    # print_table(data, columns=('protein', 'relevance', 'avg_rank', '#homologs'))

    # i = 0
    # for protein, avg_rank in sorted(avg_ranks.items(), key=lambda x: x[1]):
    #     if len(homologs[protein]) < num_species/2:
    #         continue
        
    #     if i > n:
    #         break

    #     print(protein)
    #     for homolog in homologs[protein]:
    #         print('\t- {} {} {}'.format(homolog, homologs[protein][homolog], protein_relevancies.get(homolog, None)))

    #     i += 1

    # ranks = defaultdict(lambda: {})
    # for species, ranking in species_protein_ranking.items():
    #     for protein, _ in top_proteins:
    #         homologs_in_species = [h for h in homologs[protein] if protein_species[h] == species]
    #         if len(homologs_in_species) > 0:
    #             homolog_in_species = homologs_in_species[0]
    #             index = ranking.index(homolog_in_species)
    #             ranks[protein][species] = index
                
    # print(ranks)
            
    ordered_proteins = sorted(protein_relevancies.items(), key=lambda x: x[1], reverse=True)

    top_proteins = []
    i = 0
    data = []
    for protein, relevance in ordered_proteins:
        if len(homologs[protein]) < num_species/2:
            continue
        
        if i > n:
            break
        
        top_proteins.append(protein)

        for homolog, pairwise_relevance in homologs[protein].items():
            data.append((protein, relevance, homolog, protein_relevancies.get(homolog, None), pairwise_relevance))

        i += 1
    print_table(data, columns=('protein', 'relevance', 'homolog', 'homolog_relevance', 'pairwise_relevance'))
    

    top_proteins_homologs = defaultdict(lambda: [])
    for protein in top_proteins:
        top_proteins_homologs[protein].append((protein, protein_species[protein]))
        for homolog in homologs[protein]:
            top_proteins_homologs[protein].append((homolog, protein_species[homolog]))

    return top_proteins, top_proteins_homologs


def get_consensus_seq(homologs, species):
    sequences = {}
    for specie in species:
        accession = next((h for h, s in homologs if s == specie), None)
        if accession:
            sequence = get_protein_sequence(accession, specie)
            if sequence:
                sequences[accession] = sequence
            else:
                raise RuntimeError("Can't find protein sequence {} in genome {}".format(accession, specie))
            
    if len(sequences) == 0:
        return None
    elif len(sequences) == 1:
        return list(sequences.values())[0]
    
    with tempfile.NamedTemporaryFile() as f:
        for accession, sequence in sequences.items():
            f.write(bytes(">{}\n".format(accession), 'utf-8'))
            f.write(bytes("{}\n".format(sequence), 'utf-8'))
        f.flush()
        
        cmd = ['clustalo', '--seqtype', 'Protein', '-i', f.name]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        alignment_fasta = str(result.stdout, 'utf-8')
        alignment = AlignIO.read(StringIO(alignment_fasta), "fasta")
        align_summary = SummaryInfo(alignment)
        return align_summary.dumb_consensus()


def test_proteins(proteins, homologs, test_species=None):
    anthropophily_groups = {
            'LOW': {s for s in SPECIES_ANTHROPOPHILY
                    if s != test_species
                    and SPECIES_ANTHROPOPHILY[s] <= -0.333},
            'AMBIVALENT': {s for s in SPECIES_ANTHROPOPHILY
                           if s != test_species
                           and SPECIES_ANTHROPOPHILY[s] > -0.333
                           and SPECIES_ANTHROPOPHILY[s] <= 0.333},
            'HIGH': {s for s in SPECIES_ANTHROPOPHILY
                     if s != test_species
                     and SPECIES_ANTHROPOPHILY[s] > 0.333},
    }
    
    group_scores = defaultdict(lambda: [])
    
    for protein in proteins:
        test_homolog = next((h for h, s in homologs[protein] if s == test_species), None)
        if not test_homolog:
            print('No homolog for protein {} in test specie {}'.format(protein, test_species))
            continue
        
        test_homolog_seq = get_protein_sequence(test_homolog, test_species)
        for anthropophily, species in anthropophily_groups.items():
            consensus = get_consensus_seq(homologs[protein], species)
            if not consensus:
                print('No consensus sequence for protein {} and test specie {}'.format(protein, test_species))
                continue
            aligner = PairwiseAligner(scoring='blastp')
            score = aligner.score(test_homolog_seq, consensus)
            group_scores[anthropophily].append(score)
            print(protein, anthropophily, test_homolog, score)
            
    if len(group_scores) == 0:
        return None
            
    group_avg_scores = [(group, sum(scores)/len(scores)) for group, scores in group_scores.items()]
    group = sorted(group_avg_scores, key=lambda x: x[1], reverse=True)[0]
    expected_group = "LOW"
    if SPECIES_ANTHROPOPHILY[test_species] > -0.333 and SPECIES_ANTHROPOPHILY[test_species] <= 0.333:
        expected_group = "AMBIVALENT"
    elif SPECIES_ANTHROPOPHILY[test_species] > 0.333:
        expected_group = "HIGH"
    return (group[0], expected_group)


def analyse_protein(accession):
    all_correspondences = read_correspondences()
    print('getting correspondence graph')
    g, index, rindex = get_correspondences_graph(all_correspondences)

    print('nodes in the graph', len(g.nodes))

    print('getting communities')

    communities = get_communities(g, rindex)
    community_index = get_community_membership_index(communities)
    community = communities[community_index[accession]]
    
    table_data = []
    for (genome1, genome2), correspondences in all_correspondences.items():
        for correspondence in correspondences:
            if correspondence[0] in community and correspondence[1] in community:
                other = correspondence[0] if correspondence[0] != accession else correspondence[1]
                score = float(correspondence[2])
                avg_species_score = sum(float(c[2]) for c in correspondences)/len(correspondences)
                relevance = calculate_relevance(
                    score,
                    avg_species_score,
                    SPECIES_ANTHROPOPHILY[genome1],
                    SPECIES_ANTHROPOPHILY[genome2]
                )
                table_data.append((
                    genome1,
                    SPECIES_ANTHROPOPHILY[genome1],
genome2,
                    SPECIES_ANTHROPOPHILY[genome2],
                    avg_species_score,
                    other,
                    score,
                    relevance))
                break
    print_table(table_data, columns=['S1', 'S1 anthropophily', 'S2', 'S2 anthropophily', 'avg species score', 'corresp. protein', 'protein score', 'relevance'])
    mean = sum(row[-1] for row in table_data)/len(table_data)
    variance = sum((row[-1] - mean)**2 for row in table_data)/len(table_data)
    community_score = mean/variance
    print('mean:', mean)
    print('variance:', variance)
    print('community score:', community_score)
