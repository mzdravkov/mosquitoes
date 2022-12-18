from collections import defaultdict

import networkx as nx

from storage import read_correspondences


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
            edges.add(edge)
        print('.', end='', flush=True)
        
    g = nx.Graph()
    g.add_nodes_from(index.values())
    g.add_edges_from(edges)
    return g, index, reverse_index


def get_correspondences_graph():
    correspondences = read_correspondences()
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


def print_histogram(depths):
    for depth in range(1, max(depths.keys()) + 1):
        print(f'{depth}: {depths[depth]}')
    

if __name__ == '__main__':
    print('getting correspondence graph')
    g, index, rindex = get_correspondences_graph()
    groups = {}
    depths = defaultdict(lambda: 0)
    print('\ngetting groups')
    i = 0
    for gene in index:
        if gene == '':
            continue
        if i % 1000 == 0:
            print('.', end='', flush=True)
        groups[gene], depth = get_neighbors(g, index, rindex, gene, depth=10)
        depths[depth] += 1
        if depth == 10:
            print(f'{gene} -> {depth}')
        i += 1
    print('\ndepths')
    print_histogram(depths)
    print('\ngroup sizes')
    group_sizes = defaultdict(lambda: 0)
    for group in groups.values():
        group_sizes[len(group)] += 1
    print_histogram(group_sizes)
    print(sum(group_sizes.values()))