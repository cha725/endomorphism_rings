import networkx as nx

def out_directed_subgraph(g: nx.DiGraph, gens: list) -> nx.DiGraph:
    # Input: g is a DiGraph, gens is a subset of the nodes of the graph
    # Output: the full subgraph of g induced by all the paths leaving nodes in gens

    if not set(gens) <= set(g.nodes()):
        print('The generators are not nodes of the graph.')
        return False

    if not gens:
        print('This subgraph is trivial, it has no generators.')
        return nx.DiGraph()

    subgraph_nodes = gens # include the generators in the new vertices

    for target in g.nodes():
        for source in gens:
            if nx.has_path(g,source,target) and source != target:
                # nx says that a vertex has a path from itself to itself
                subgraph_nodes.append(target)
                break
                # if there is a path source to target,
                # then we include the target and do not need to check for anymore paths
                # so we break instead of checking more sources.

    return nx.induced_subgraph(g, subgraph_nodes)

def in_directed_subgraph(g: nx.DiGraph, gens: list) -> nx.DiGraph:
    # Input: g is a DiGraph, gens is a subset of the nodes of the graph
    # Output: the full subgraph of g induced by all the paths finishing at nodes in gens

    if not set(gens) <= set(g.nodes()):
        print('The generators are not nodes of the graph.')
        return False

    if not gens:
        print('This subgraph is trivial, it has no generators.')
        return nx.DiGraph()

    subgraph_nodes = gens  # include the generators in the new vertices

    for source in g.nodes():
        for target in gens:
            if nx.has_path(g, source, target) and source != target:
                # nx says that a vertex has a path from itself to itself
                subgraph_nodes.append(source)
                break
                # if there is a path source to target,
                # then we include the source and do not need to check for anymore paths
                # so we break instead of checking more sources.

    return nx.induced_subgraph(g, subgraph_nodes)

def sources(g:nx.DiGraph) -> list:
    # Input: a nx.DiGraph
    # Output: a list of the indices of source vertices
    return [x for x in g.nodes() if g.in_degree(x) == 0]

def sinks(g:nx.DiGraph) -> list:
    # Input: a nx.DiGraph
    # Output: a list of the indices of the sink vertices
    return [x for x in g.nodes() if g.out_degree(x) == 0]

def create_radical_layers(g:nx.DiGraph, start=None) -> dict:
    # Input: (1) a nx.DiGraph
    #        (2) a list of starting vertices
    #            (optional arg: the default is the list of sources of the DiGraph)
    # Output: a dictionary that encodes the length of the shortest path to each vertex from a starting vertex
    if start is None:
        start = sources(g)
    return dict(enumerate(nx.bfs_layers(g, start)))

def draw_layered_graph(g:nx.DiGraph, start=None):
    # Input: (1) a nx.DiGraph
    #        (2) a list of starting vertices
    #            (optional arg: the default is the list of sources of the DiGraph)
    # Output: a diagram of the DiGraph that is layered with the starting vertices at the top

    if start is None:
        start = sources(g)
    # If an optional arg is not given, initialise it to the sources of the graph

    G = nx.DiGraph

    for layer in create_radical_layers(g, start):
        for vertex in layer:
            G.add_node(vertex, layer=layer, label=g.nodes[vertex]['label'])

    for edge in g.edges:
        G.add_edge(*edge)
    return G


if __name__ != "__main__":
    pass
else:
    G = nx.DiGraph()
    G.add_edges_from([(0,1),(1,2),(2,3),(3,4),(4,5),(0,3),(2,1),(5,2)])
    print(in_directed_subgraph(G,[1]))
    print(out_directed_subgraph(G,[1]))
    print(draw_layered_graph(G))
