import networkx as nx
import matplotlib.pyplot as plt

# Module operations on graphs

def out_directed_subgraph(g: nx.MultiDiGraph, gens: list) -> nx.MultiDiGraph:
    # Input: g is a MultiDiGraph, gens is a subset of the nodes of the graph
    # Output: the full subgraph of g induced by all the paths leaving nodes in gens

    if not set(gens) <= set(g.nodes()):
        print('The generators are not nodes of the graph.')
        return nx.MultiDiGraph()

    subgraph_nodes = set(gens)  # Use a set to prevent duplicates

    for target in g.nodes():
        for source in gens:
            if nx.has_path(g, source, target) and source != target:
                subgraph_nodes.add(target)
                break  # No need to check other sources

    return g.subgraph(subgraph_nodes)


def in_directed_subgraph(g: nx.MultiDiGraph, gens: list) -> nx.MultiDiGraph:
    # Input: g is a MultiDiGraph, gens is a subset of the nodes of the graph
    # Output: the full subgraph of g induced by all the paths finishing at nodes in gens

    if not set(gens) <= set(g.nodes()):
        print('The generators are not nodes of the graph.')
        return nx.MultiDiGraph()

    subgraph_nodes = set(gens)  # Use a set to prevent duplicates
    for source in g.nodes():
        for target in gens:
            if nx.has_path(g, source, target) and source != target:
                subgraph_nodes.add(source)
                break  # No need to check other targets

    return g.subgraph(subgraph_nodes)


def sources(g: nx.MultiDiGraph) -> list:
    # Input: a nx.MultiDiGraph
    # Output: a list of the indices of source vertices
    return [x for x in g.nodes() if g.in_degree(x) == 0]


def sinks(g: nx.MultiDiGraph) -> list:
    # Input: a nx.MultiDiGraph
    # Output: a list of the indices of the sink vertices
    return [x for x in g.nodes() if g.out_degree(x) == 0]


def create_radical_layers(g: nx.MultiDiGraph, start=None) -> dict:
    # Input: (1) a nx.MultiDiGraph
    #        (2) a list of starting vertices
    #            (optional arg: defaults to the sources of the MultiDiGraph)
    # Output: a dictionary encoding the length of the shortest path to each vertex from a starting vertex
    if start is None:
        start = sources(g)

    layers = {i: list(layer) for i, layer in enumerate(nx.bfs_layers(g, start))}
    return layers

def relabelling_map_wrt_radical_layers(g:nx.MultiDiGraph, start=None) -> dict:
    # Input: (1) a nx.MultiDiGraph
    #        (2) a list of starting vertices
    #            (optional arg: defaults to the sources of the MultiDiGraph)
    # Output: a relabeling map of the vertices of the nx.MultiDigraph
    #         labels are taken in order of radical layers

    layers = create_radical_layers(g, start)

    compressed = []
    for layer in layers:
        compressed += layers[layer]

    relabel_map = {}
    for vertex in compressed:
        relabel_map[vertex] = compressed.index(vertex)

    return relabel_map

def relabelling_wrt_radical_layers(g:nx.MultiDiGraph, start=None) -> nx.MultiDiGraph:
    relabel_map = relabelling_map_wrt_radical_layers(g)
    g = nx.relabel_nodes(g, relabel_map)
    layers = create_radical_layers(g)
    
    return create_layered_graph(g)


def create_layered_graph(g: nx.MultiDiGraph, start=None) -> nx.MultiDiGraph:
    # Input: (1) a nx.MultiDiGraph
    #        (2) a list of starting vertices
    #            (optional arg: defaults to the sources of the MultiDiGraph)
    # Output: a relabelling of the MultiDiGraph that includes a layer for each vertex

    if start is None:
        start = sources(g)

    layer_g = nx.MultiDiGraph()  # Properly initialize the graph

    # Preserve node attributes
    for layer, vertices in create_radical_layers(g, start).items():
        # .items() means the for loop goes over both the keys and the values
        # a bit like enumerate for lists
        for vertex in vertices:
            if vertex in g.nodes:  # Only preserve existing attributes
                layer_g.add_node(vertex, layer=layer, **g.nodes[vertex])

    # Preserve edge attributes
    for u, v, key, data in g.edges(keys=True, data=True):
        layer_g.add_edge(u, v, key=key, **data)

    return layer_g

def draw(g:nx.MultiDiGraph): # draws the networkx graph nicely
    layer_g = create_layered_graph(g)
    pos = {}
    labels = {}
    for layer_id, layer in create_radical_layers(g).items():
        for vertex_id, vertex in enumerate(layer):
            pos[vertex] = (vertex_id, -layer_id)
            label = layer_g.nodes[vertex].get('label', None)
            if label is None:
                labels[vertex] = vertex
            else:
                labels[vertex] = vertex, label

    nx.draw(g, pos, labels=labels, with_labels=True, node_size=1000, node_color='lightblue', font_size=10,
            font_weight='bold', arrows=True)
    plt.show()

def draw_subgraph(g:nx.MultiDiGraph, vert):

    if not set(vert) <= set(g.nodes()):
        return

    node_colors = ['red' if node in vert else 'lightblue' for node in g.nodes()]
    layer_g = create_layered_graph(g)
    pos = {}
    labels = {}
    for layer_id, layer in create_radical_layers(g).items():
        for vertex_id, vertex in enumerate(layer):
            pos[vertex] = (vertex_id, -layer_id)
            label = layer_g.nodes[vertex].get('label', None)
            if label is None:
                labels[vertex] = vertex
            else:
                labels[vertex] = (vertex, label)

    nx.draw(g, pos, labels=labels, with_labels=True, node_size=1000, node_color=node_colors, font_size=10,
            font_weight='bold', arrows=True)
    plt.show()

def draw_out_directed_subgraph(g: nx.MultiDiGraph, gens: list):
    # Input: A MultiDiGraph and a list of generators
    # Output: Plot of the out directed subgraph of g generated by gens
    subgraph = out_directed_subgraph(g, gens)
    draw_subgraph(g, list(subgraph.nodes()))
    plt.show()

def draw_in_directed_subgraph(g: nx.MultiDiGraph, gens: list):
    # Input: A MultiDiGraph and a list of generators
    # Output: Plot of the in directed subgraph of g generated by gens
    subgraph = in_directed_subgraph(g, gens)
    draw_subgraph(g, list(subgraph.nodes()))
    plt.show()


if __name__ == "__main__":

    G = nx.MultiDiGraph()
    G.add_edges_from([
        (0, 1, {'label': "a"}), (1, 2, {'label': "b"}),
        (2, 3, {'label': "c"}), (3, 4, {'label': "d"}),
        (4, 5, {'label': "e"}), (0, 3, {'label': "f"}),
        (2, 1, {'label': "g"}), (5, 2, {'label': "h"})
    ])

    G.nodes[0]['label'] = "start"
    G.nodes[5]['label'] = "end"
    print(create_radical_layers(G))
    print("\nLayered Graph Representation:", create_layered_graph(G).nodes)

    h = nx.MultiDiGraph()
    h.add_edges_from([(1,2),(1,3),(4,2),(3,4)])
    print(create_radical_layers(h))
    print(relabelling_wrt_radical_layers(h).edges())

