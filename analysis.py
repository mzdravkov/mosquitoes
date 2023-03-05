import random
import statistics

from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from io import StringIO
import subprocess
import tempfile
import pprint

from sklearn.metrics import silhouette_score
import networkx as nx
from proteins import get_protein_sequence

from storage import read_correspondences
from utils import print_table

from Bio  import AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align.AlignInfo import SummaryInfo

# TODO: Get anthropophily values for all species

# Score indicating the degree of anthropophily of a species
# in the range of [-1, 1], where -1 indicates a species that
# is non-anthropophilic and 1 indicates a species that is
# very anthropophilic.
# SPECIES_ANTHROPOPHILY = {
# 'GCA_000441895.2': -0.918, # Anopheles sinensis [Ree, Han-Il, et al. 2001][PMC2712014]
# 'GCF_000005575.2': 0.52, # Anopheles gambiae str. PEST [takken (PMC4381365)]
# 'GCF_006496715.1': 0.96, # Aedes albopictus [Alongkot Ponlawat 2005]
# 'GCF_013141755.1': -0.982, # Anopheles stephensi [Thomas, S., Ravishankaran 2017]
# 'GCF_013758885.1': -0.9, # Anopheles albimanus [Bruce-Chwatt 1966 (PMC2476083)]
# 'GCF_015732765.1': -0.33, # Culex quinquefasciatus [takken (PMC4381365)]
# 'GCF_016801865.2': -0.286, # Culex pipiens pallens [Joaquín Muñoz 2011]
# 'GCF_016920715.1': -0.09, # Anopheles arabiensis [tekken (PMC4381365)]
# 'GCF_017562075.2': -0.76, # Anopheles merus [Pamela C Kipyab 2013]
# 'GCF_943734635.1': 1.0, # Anopheles cruzii [Kirchgatter 2014, Santos 2019]
# 'GCF_943734655.1': , # Sabethes cyaneus
# 'GCF_943734665.1': 0.0, # Anopheles aquasalis
# 'GCF_943734685.1': 0.96, # Anopheles coluzzii [Martin C. Akogbéto 2018]
# 'GCF_943734695.1': 0.0, # Anopheles maculipalpis
# 'GCF_943734725.1': , # Anopheles marshallii
# 'GCF_943734745.1': 0.65, # Anopheles darlingi [Marta Moreno 2017]
# 'GCF_943734755.1': 1.0, # Anopheles moucheti [Sinka 2010]
# 'GCF_943734845.2': 0.21, # Anopheles funestus [tekken (PMC4381365)]
# 'GCF_943737925.1': , # Anopheles nili
# }

SPECIES_ANTHROPOPHILY = {
'GCA_000441895.2': -1.0, # Anopheles sinensis [Ree, Han-Il, et al. 2001][PMC2712014]
'GCF_000005575.2': 1.0, # Anopheles gambiae str. PEST [takken (PMC4381365)]
'GCF_006496715.1': 1.0, # Aedes albopictus [Alongkot Ponlawat 2005]
'GCF_013141755.1': -1.0, # Anopheles stephensi [Thomas, S., Ravishankaran 2017]
'GCF_013758885.1': -1.0, # Anopheles albimanus [Bruce-Chwatt 1966 (PMC2476083)]
'GCF_015732765.1': -1.0, # Culex quinquefasciatus [takken (PMC4381365)]
'GCF_016801865.2': -1.0, # Culex pipiens pallens [Joaquín Muñoz 2011]
'GCF_016920715.1': -1.0, # Anopheles arabiensis [tekken (PMC4381365)]
'GCF_017562075.2': -1.0, # Anopheles merus [Pamela C Kipyab 2013]
'GCF_943734635.1': 1.0, # Anopheles cruzii [Kirchgatter 2014, Santos 2019]
'GCF_943734655.1': 1.0, # Sabethes cyaneus [Leticia Smith 2023]
'GCF_943734665.1': 1.0, # Anopheles aquasalis
'GCF_943734685.1': 1.0, # Anopheles coluzzii [Martin C. Akogbéto 2018]
'GCF_943734695.1': 1.0, # Anopheles maculipalpis
'GCF_943734725.1': 1.0, # Anopheles marshallii [Boris Makanga 2016]
'GCF_943734745.1': 1.0, # Anopheles darlingi [Marta Moreno 2017]
'GCF_943734755.1': 1.0, # Anopheles moucheti [Sinka 2010]
'GCF_943734845.2': 1.0, # Anopheles funestus [tekken (PMC4381365)]
'GCF_943737925.1': 1.0, # Anopheles nili [Antonio-Nkondjio 2013 (Anopheles mosquitoes book)]
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
    # g.remove_node(index[''])
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

    # return (genes_similarity/(species_similarity + 1)) ** anthropophily_similarity

    weight = 1 - species_similarity/100.0 if anthropophily_similarity == 1 else species_similarity/100.0
    return weight * anthropophily_similarity * (genes_similarity - species_similarity)


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
    """
    Makes sure that the gene is not already present in a
    genome different from the specified one.
    """
    if gene in protein_species and protein_species[gene] != genome:
        raise RuntimeError(
            "Gene {} in two species: {} and {}".format(
                gene, protein_species[gene], genome)
        )


@lru_cache(maxsize=None)
def top_by_relevance(n, test_species=None):
    """
    Finds a list of the N best protein candidates
    based on their relevance scores.
    If given a test_species it will ignore this species from the
    analysis so that it can be used to test the results.
    """
    correspondences = read_correspondences()
    # Map of type: protein => its_specie
    protein_species = {}
    # Map of type: protein => {homolog: alignment_score}
    homologs = defaultdict(lambda: {})
    # Set of all species (genome accessions)
    species = set()

    null_protein = 0

    # Calculate the relevance scores between each protein and every one of its homologs
    # (note: we only work with the best homolog to the protein from every species)
    for (genome1, genome2), correspondences in correspondences.items():
        ignore_relevancies = False
        if genome1 == test_species or genome2 == test_species:
            ignore_relevancies = True

        species1_anthropophily = SPECIES_ANTHROPOPHILY[genome1]
        species2_anthropophily = SPECIES_ANTHROPOPHILY[genome2]

        print('.', end='', flush=True)
        total_species_score = sum(float(score) for _, _, score in correspondences)
        avg_species_score = total_species_score / len(correspondences)

        for correspondence in correspondences:
            gene1, gene2, score_str = correspondence

            if gene2 == '':
                gene2 = 'null_' + str(null_protein)
                null_protein += 1

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

            # if gene1 == '' or gene2 == '':
            #     continue

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

            # If one of the species is a test specie
            # we should ignore the relevances from
            # comparisons against it.
            if ignore_relevancies:
                homologs[gene1][gene2] = None
                # homologs[gene2][gene1] = None
            else:
                homologs[gene1][gene2] = relevance
                # homologs[gene2][gene1] = relevance

    protein_relevancies = {}

    # Each protein is assigned the average of all the relevance
    # scores of itself and its homologs.
    for gene, relevancies in homologs.items():
        filtered_relevancies = [r for r in relevancies.values() if r is not None]
        if len(filtered_relevancies) > 0:
            mean = sum(filtered_relevancies)/len(filtered_relevancies)
            protein_relevancies[gene] = mean
        else:
            # TODO: this is a bit random, but I'm not sure how to handle those cases
            protein_relevancies[gene] = 0.0

    # Make sure that a protein is truly relevant by taking not only its own
    # relevance score, but also the relevance score of its homologs.
    consistent_relevancies = {}
    for gene, gene_relevance in protein_relevancies.items():
        filtered_homologs = [homolog for homolog, score in homologs[gene].items() if score is not None and homolog != '']
        homolog_sum = sum(protein_relevancies[homolog] for homolog in filtered_homologs)
        relevance_sum = homolog_sum + gene_relevance
        consistent_relevance = relevance_sum/(len(filtered_homologs) + 1)
        consistent_relevancies[gene] = consistent_relevance

    num_species = len(species)

    ordered_proteins = sorted(consistent_relevancies.items(), key=lambda x: x[1], reverse=True)
    # ordered_proteins = sorted(protein_relevancies.items(), key=lambda x: x[1], reverse=True)

    # Get the N top proteins, ignoring those that are not sufficiently represented
    top_proteins = []
    included_as_homolog = set()
    i = 0
    data = []
    for protein, relevance in ordered_proteins:
        if i >= n:
            break

        non_zero_homologs = [h for h, score in homologs[protein].items() if not h.startswith('null') and score and score > 0]
        # If it is not present in at least a third of the
        # species, then we ignore it.
        if len(non_zero_homologs) + 1 < num_species/3:
            continue

        if protein in included_as_homolog:
            continue

        top_proteins.append(protein)
        included_as_homolog.update(homologs[protein].keys())

        for homolog, pairwise_relevance in homologs[protein].items():
            prot_species = protein_species.get(protein, '')
            homolog_species = protein_species.get(homolog, '')
            relevance_str = '{:.2f}'.format(relevance)
            # homolog_relevance = protein_relevancies.get(homolog)
            homolog_relevance = consistent_relevancies.get(homolog)
            homolog_relevance_str = '{:.2f}'.format(homolog_relevance) if homolog_relevance else ''
            pairwise_relevance_str = '{:.2f}'.format(pairwise_relevance) if pairwise_relevance else ''
            data.append((prot_species, protein, relevance_str, homolog_species, homolog, homolog_relevance_str, pairwise_relevance_str))

        i += 1

    print('')
    print_table(data, columns=('protein_species', 'protein', 'relevance', 'homolog_species', 'homolog', 'homolog_relevance', 'pairwise_relevance'))

    top_proteins_homologs = defaultdict(lambda: [])
    for protein in top_proteins:
        top_proteins_homologs[protein].append((protein, protein_species[protein]))
        for homolog in homologs[protein]:
            if homolog == '':
                continue
            top_proteins_homologs[protein].append((homolog, protein_species[homolog]))

    control_proteins = random.sample(homologs.keys(), n)
    control_homologs = defaultdict(lambda: [])
    for protein in control_proteins:
        control_homologs[protein].append((protein, protein_species[protein]))
        for homolog in homologs[protein]:
            if homolog == '':
                continue
            control_homologs[protein].append((homolog, protein_species[homolog]))

    return top_proteins, top_proteins_homologs, control_proteins, control_homologs


def get_homolog_seqs(homologs, species):
    """
    Get the sequences for the homologs that are in the given species.
    """
    sequences = {}
    for specie in species:
        accession = next((h for h, s in homologs if s == specie and not h.startswith('null')), None)
        if accession:
            sequence = get_protein_sequence(accession, specie)
            if sequence:
                sequences[accession] = sequence
            else:
                raise RuntimeError("Can't find protein sequence {} in genome {}".format(accession, specie))
    return sequences


def get_consensus_seq(homologs, species):
    """
    Constructs a consensus sequence from the subset
    of homologs found in a list of species.
    """
    sequences = get_homolog_seqs(homologs, species)

    if len(sequences) == 0:
        return None
    elif len(sequences) == 1:
        return list(sequences.values())[0]

    # Write the sequences to a fasta file and use Clustal Omega to
    # align them. Use the alignment to generate the consensus sequence.
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
        return align_summary.dumb_consensus(threshold=0.5)


def get_anthropophily_groups(test_species=None):
    return {
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


def get_expected_species_group(test_species):
    expected_group = "LOW"
    if SPECIES_ANTHROPOPHILY[test_species] > -0.333 and SPECIES_ANTHROPOPHILY[test_species] <= 0.333:
        expected_group = "AMBIVALENT"
    elif SPECIES_ANTHROPOPHILY[test_species] > 0.333:
        expected_group = "HIGH"
    return expected_group


def test_proteins(proteins, homologs, test_species=None):
    anthropophily_groups = get_anthropophily_groups(test_species)

    group_scores = defaultdict(lambda: 0)

    for protein in proteins:
        test_homolog = next((h for h, s in homologs[protein] if s == test_species and not h.startswith('null')), None)
        if not test_homolog:
            # If we don't have a homologous protein in the test specie we cannot compare it
            # to the consensus sequences of the protein in the different anthropophily groups.
            # Instead, we check if some of the groups also misses this protein and choose that
            # group as "more similar".
            group_occurances = []
            for anthropophily, species in anthropophily_groups.items():
                if len(species) == 0:
                    continue
                sequences = get_homolog_seqs(homologs[protein], species)
                group_occurances.append((anthropophily, len(sequences)))
            best_group = sorted(group_occurances, key=lambda x: x[1])[0][0]
            group_scores[best_group] += 1

            print('No homolog for protein {} in test specie {}'.format(protein, test_species))
            print('Using an alteranative strategy to pick best group: ', best_group)
            continue

        test_homolog_seq = get_protein_sequence(test_homolog, test_species)
        max_score = 0
        best_group = None
        for anthropophily, species in anthropophily_groups.items():
            consensus = get_consensus_seq(homologs[protein], species)
            if not consensus:
                print('No consensus sequence for protein {} and test specie {}'.format(protein, test_species))
                continue
            aligner = PairwiseAligner(scoring='blastp')
            try:
                score = aligner.score(test_homolog_seq, consensus)
            except:
                print(test_homolog_seq)
                print(consensus)
            # group_scores[anthropophily].append(score)
            if score > max_score:
                max_score = score
                best_group = anthropophily
            print(protein, anthropophily, test_homolog, score)
        if best_group:
            group_scores[best_group] += 1

    if len(group_scores) == 0:
        return None

    # group_avg_scores = [(group, sum(scores)/len(scores)) for group, scores in group_scores.items()]
    # group = sorted(group_avg_scores, key=lambda x: x[1], reverse=True)[0]
    group = sorted(group_scores.items(), key=lambda x: x[1], reverse=True)[0]
    expected_group = get_expected_species_group(test_species)
    return (test_species, group[0], expected_group)


def get_top_proteins_and_validate(proteins_num, window=None):
    proteins = set()
    results = []
    control_results = []

    i = 0
    for test_specie in SPECIES_ANTHROPOPHILY:
        print('Progress: {}/{}'.format(i, len(SPECIES_ANTHROPOPHILY)))
        i += 1
        top_proteins, homologs, control_proteins, control_homologs = top_by_relevance(proteins_num, test_species=test_specie)
        if window:
            top_proteins = top_proteins[window[0]:window[1]]
            control_proteins = control_proteins[window[0]:window[1]]
        proteins.update(top_proteins)
        result = test_proteins(top_proteins, homologs, test_specie)
        print("RESULT: ", result)
        if result:
            results.append(result)
        control = test_proteins(control_proteins, control_homologs, test_specie)
        if control:
            control_results.append(control)

    print('Proteins:')
    for protein in proteins:
        print(protein)

    print('Validation results:')
    print_table(results, columns=('Test Species', 'Predicted', 'Expected'))
    accuracy = len([x for x in results if x[1] == x[2]])/len(results)
    print('Prediction accuracy:', accuracy)
    print("Control: ", len([x for x in control_results if x[1] == x[2]])/len(control_results))
    return accuracy, proteins


def moving_window_validation(proteins_num, window_size, window_step):
    start = 0
    end = window_size
    accuracies = []
    collected_proteins = set()
    while end < proteins_num:
        print('Window: {}-{}'.format(start, end))
        window = (start, end)
        accuracy, proteins = get_top_proteins_and_validate(proteins_num, window)
        collected_proteins.update(proteins)
        accuracies.append((window, accuracy))
        start += window_step
        end += window_step
    print("Collected proteins:")
    for protein in proteins:
        print(protein)
    pprint.pprint(accuracies)
    return accuracies


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


def get_clustering_silhouette_score():
    correspondences = read_correspondences()

    avg_species_scores = {}
    for (genome1, genome2), correspondences in correspondences.items():
        total_species_score = sum(float(score) for _, _, score in correspondences)
        avg_species_score = total_species_score / len(correspondences)
        avg_species_scores[(genome1, genome2)] = avg_species_score

    anthropophily_groups = get_anthropophily_groups()

    labels = []
    # Get a list of labels
    for species in SPECIES_ANTHROPOPHILY:
        for group, group_species in anthropophily_groups.items():
            if species in group_species:
                labels.append(group)

    distances = []
    # Get a matrix with distances between pairs of species
    for species1 in SPECIES_ANTHROPOPHILY:
        row = []
        for species2 in SPECIES_ANTHROPOPHILY:
            if species2 == species2:
                row.append(0)
            else:
                score = avg_species_scores.get((species1, species2), avg_species_scores[(species2, species1)])
                distance = 1/score
                row.append(distance)
        distances.append(row)

    return silhouette_score(distances, labels, metric='precomputed')


def avg_score_for_phenotypic_groups():
    correspondences = read_correspondences()
    avg_species_scores = {}
    for (genome1, genome2), correspondences in correspondences.items():
        total_species_score = sum(float(score) for _, _, score in correspondences)
        avg_species_score = total_species_score / len(correspondences)
        avg_species_scores[(genome1, genome2)] = avg_species_score

    anthropophily_groups = get_anthropophily_groups()
    total_average = sum(avg_species_scores.values())/len(avg_species_scores.values())
    print('Average score between all pairs of species:', total_average)
    group_pairs = combinations_with_replacement({k for k, v in anthropophily_groups.items() if len(v) > 0}, 2)
    for (g1, g2) in group_pairs:
        group1 = anthropophily_groups[g1]
        group2 = anthropophily_groups[g2]
        scores = []
        for (specie1, specie2), avg_score in avg_species_scores.items():
            if (specie1 in group1 and specie2 in group2) or (specie1 in group2 and specie2 in group1):
                scores.append(avg_score)
        if len(scores) > 0:
            print('Average score of pairs between group {} and {}: {}'.format(g1, g2, sum(scores)/len(scores)))

    print('Silhuette score:', get_clustering_silhouette_score())
