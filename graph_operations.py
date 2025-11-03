import networkx as nx
#from itertools import combinations

#from networkx.algorithms.shortest_paths.unweighted import predecessor
#from networkx.classes import all_neighbors

def create_radical_layers(G :nx.DiGraph): # G is a graph in networkx
    # A vertex is in radical layer l if it is l steps away from a vertex that is a source
    sources = [x for x in G.nodes() if G.in_degree(x) == 0]
    if sources == []:
        print('This nx.DiGraph has no vertices that are sources.')
        return False
    radical_layers = [sources]
    test_vertices = sources
    vertices_left = list(G.nodes())
    while vertices_left != []: # This is to make sure each vertex is assigned a unique radical level
        next_neighbors = []
        for vertex in test_vertices:
            next_neighbors += [x for x in G.neighbors(vertex) if x in vertices_left]
            vertices_left.remove(vertex)
        if next_neighbors != []:
            radical_layers += [next_neighbors]
            test_vertices = next_neighbors
        else:
            return radical_layers

def find_top(G :nx.DiGraph):
    # The top of a DiGraph is the set of vertices that are sources
    return create_radical_layers(G)[0]

def find_socle(G :nx.DiGraph):
    # The socle of a DiGraph is the set of vertices that are sinks
    socle = [x for x in G.nodes() if G.out_degree(x) == 0]
    if socle == []:
        print('This nx.DiGraph has no vertices that are sinks.')
        return False
    else:
        return socle

def relabel_wrt_radical_layers(G :nx.DiGraph):
    # This changes the labelling of the vertices in the DiGraph so that the vertices in radical layer l are strictly smaller than those in radical layer l+1
    radical_layers = create_radical_layers(G)
    mapping = {}
    counter = 0
    for layer in radical_layers:
        for vertex in layer:
            mapping[vertex] = counter
            counter += 1
    new_edges = [(mapping[i], mapping[j]) for i, j in list(G.edges())]
    return new_edges

if __name__ == "__main__":
    G = nx.DiGraph()
    G.add_nodes_from(range(5))
    G.add_edges_from([(0,1),(0,4),(2,3),(3,0),(2,5),(5,3),(1,3),(4,1)])
    print('nodes ', list(G.nodes()))
    print('edges ', list(G.edges()))
    print('radical layers ', create_radical_layers(G))
    print('relabelling ', relabel_wrt_radical_layers(G))
    print('find top ', find_top(G))
    print('find socle ', find_socle(G))
