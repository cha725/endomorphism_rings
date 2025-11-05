from itertools import permutations

from modulediagram import ModuleDiagram, Arrow
import networkx as nx
from collections import deque
import matplotlib.pyplot as plt

def composition(m1 : dict, m2 : dict):
    return {k: m2[v] for k, v in m1.items() if v in m2}


class SubModuleDiagram:

    def __init__(self, super_module, sub_module, inclusion_map):
        self.super_module = super_module # a ModuleDiagram object
        self.sub_module = sub_module # a ModuleDiagram object
        self.inclusion_map = inclusion_map # a list of tuples





class AuslanderAlgebra:

    def __init__(self,*args):
        self.args = args
        self.dim = len(self.args)
        self.original_modules = self.args
        self.vertices = self.vertex_graph()
        self.vertices_explicit = self.vertex_graph().nodes(data=True)

    def vertex_graph(self):
        """
        The vertices of the Auslander algebra are all the indecomposable summands
        of the projective modules.
        """
        G = nx.MultiDiGraph()
        counter = 0
        for mod in self.args:
            for sub in mod.generate_all_submodules():
                if sub[1].is_isomorphic_to(mod):
                    G.add_node(counter, original = True, structure = sub, label = chr(counter+65))
                else:
                    G.add_node(counter, original = False, structure = sub, label = chr(counter + 65))
                counter += 1
        return G


    def all_edges_graph(self):
        G = self.vertex_graph()
        for node1, node2 in permutations(G.nodes(), 2):  # Iterate over all pairs
            struct1 = nx.get_node_attributes(G, "structure")[node1][1]
            struct2 = nx.get_node_attributes(G, "structure")[node2][1]
            for hom in struct1.homomorphism_group(struct2):
                G.add_edge(node1, node2, map=hom)
        return G

    def remove_redundant_edges(self):
        G = self.all_edges_graph
        edges_to_remove = []

        for node1, node2, edge_data in G.edges(data=True):
            direct_map = edge_data['map']
            print(direct_map)

            for path in nx.all_simple_paths(G, source=node1, target=node2, cutoff=3):  # Limit path length
                composed_map = path[0]  # Initialize with first map
                print(path)
                valid_composition = True

                for i in range(len(path) - 1):
                    u, v = path[i], path[i + 1]
                    edge_map = G[u][v][0]['map']  # Get first edge map



                    composed_map = composition(composed_map, edge_map)  # Apply composition

                    print(composed_map)
                    if composed_map is None:
                        valid_composition = False
                        break

                if valid_composition and composed_map == direct_map:
                    edges_to_remove.append((node1, node2))

        # Remove redundant edges
        G.remove_edges_from(edges_to_remove)
        return G


    # def find_walks(self, start):
    #     G = self.all_edges_graph  # Graph
    #     ind_arrow_maps = [x[2]['map'] for y,x in enumerate(G.edges(data=True))]
    #     queue = deque()  # BFS queue
    #     max_paths = []  # Stores maximal paths
    #     # Initialize the queue with paths starting from `start`
    #     for edge in G.edges(start, data=True):
    #         queue.append([edge])  # Each path starts as a list with one edge
    #     all_paths = []  # Stores all valid paths
    #     while queue:
    #         path = queue.popleft()  # Get the next path to explore
    #         last_node = path[-1][1]  # Last node in the path
    #         new_arrows = list(G.edges(last_node, data=True))
    #         if not new_arrows:
    #             max_paths.append(path)  # No further edges, store as maximal
    #             continue
    #         maximal = True
    #         for new_arrow in new_arrows:
    #             comp = composition(path[-1][2]['map'], new_arrow[2]['map'])  # Ensure valid composition
    #             if comp:
    #                 new_path = path + [new_arrow]  # Create a new path
    #                 queue.append(new_path)  # Add to queue for further exploration
    #                 maximal = False
    #         if maximal:
    #             max_paths.append(path)
    #         all_paths.append(path)
    #     return all_paths, max_paths

    # def find_walks(self, start, max_length=None):
    #     G = self.all_edges_graph
    #     queue = deque()
    #     max_paths = []
    #
    #     for edge in G.edges(start, data=True):
    #         queue.append([edge])
    #
    #     while queue:
    #         path = queue.popleft()
    #         last_node = path[-1][1]
    #
    #         if max_length and len(path) >= max_length:
    #             max_paths.append(path)
    #             continue
    #
    #         found_extension = False
    #         for new_arrow in G.edges(last_node, data=True):  # Use iterator for efficiency
    #             comp = composition(path[-1][2]['map'], new_arrow[2]['map'])
    #             if comp:
    #                 queue.append(path + [new_arrow])
    #                 found_extension = True
    #
    #         if not found_extension:
    #             max_paths.append(path)
    #     return max_paths

    def find_walks(self, start, max_length=None):
        G = self.all_edges_graph
        queue = deque()
        max_paths = []

        for edge in G.edges(start, data=True):
            queue.append((edge,))

        while queue:
            path = queue.popleft()
            last_node = path[-1][1]

            if max_length and len(path) >= max_length:
                max_paths.append(path)
                continue

            found_extension = False
            for new_arrow in G.edges(last_node, data=True):
                # Compute the composition of the last edge in path with the new edge
                comp = composition(path[-1][2]['map'], new_arrow[2]['map'])

                if comp:  # Stop if composition is zero
                    queue.append((*path, new_arrow))
                    found_extension = True

            if not found_extension:
                max_paths.append(path)

        return max_paths

    def inj_diagram(self):
        diagrams = []
        for vertex in self.vertices_explicit:
            G = nx.DiGraph()
            max_paths = self.find_walks(vertex[0])
            id_map = {}
            for factor in range(vertex[1]['structure'][1].num_vertices):
                id_map[factor] = factor
            G.add_node(0, data=(vertex, id_map))
            for path in max_paths:
                comp = id_map
                for arr in path:
                    comp = composition(comp, arr[2]['map'])
                    target = arr[1]
                    L = len(list(G.nodes()))
                    G.add_node(L, data=(self.vertices_explicit[target], comp))
                    G.add_edge(L-1,L, data=arr)
            diagrams.append(G)
        return diagrams

    # def draw_inj_diagrams(self): # draws the networkx graph nicely
    #     diagrams = self.inj_diagram()
    #
    #
    #
    #     nx.draw(self.graph, labels=labels, with_labels=True, node_size=1000, node_color='lightblue', font_size=10,
    #             font_weight='bold', arrows=True)
    #     plt.show()

# Example usage
if __name__ == "__main__": # can import it to another file and nothing below this line will run
    # this is where you put test cases so you can see what the code does

    # M = ModuleDiagram(['1','1','3','2'], [(0,1),(0,2),(2,3)])
    # L = ModuleDiagram(['2','1','1','3'], [(2,1),(2,3),(3,0)])
    # N = ModuleDiagram(['2', 'A', '1', '3'], [(2, 1), (2, 3), (3, 0)])
    # Q = ModuleDiagram(['3','1'],[(1,0)])

    #print(Q.radical_layers)
    #M.draw()
    #print(M.create_sub_module([3]))
    #print('hi', [M.create_sub_module([1])[i].graph.nodes(data=True) for i in (0, 1)])
    #print(L.create_sub_module([1])[2])
    #print('\n', 'module, quotient module ', [M.create_quotient_module([0])[i].graph.nodes(data=True) for i in (0, 1)])
    #print('\n', 'quotient map ', M.create_quotient_module([0])[2])
    #print('isomorphism ', M.create_sub_module([1])[1].is_isomorphic_to(M.create_quotient_module([0])[1]))
    #print(M.top)
    #print(M.socle)
    #M.draw()
    #L.draw()
    #N.draw()
    #print('isomorphism ', M.is_isomorphic_to(N))
    #print('isomorphism ', M.is_isomorphic_to(L))
    #print('nodes', list(M.graph.nodes(data=True)))
    #print('submodules', len(M.generate_all_sub_modules()))
    #print('submodules', [M.generate_all_sub_modules()[i][2] for i in range(len(M.generate_all_sub_modules()))])
    #print('quotient modules', [M.generate_all_quotient_modules()[i][2] for i in range(len(M.generate_all_quotient_modules()))])
    #print('original vertex labels ', M.vertex_labels)
    #print('radical layers ', M.radical_layers)
    #print('radical labels', M.radical_labels)
    #print('hom', M.homomorphism_group(M))


    MA = ModuleDiagram(
        arrows=[Arrow(0,1), Arrow(0,1), Arrow(2,3)],
        vertex_composition={0:1, 1:1, 2:3, 3:2})
    MB = ModuleDiagram(
        arrows=[Arrow(0,1), Arrow(0,2), Arrow(0,3), Arrow(1,4), Arrow(2,4), Arrow(2,5), Arrow(3,5)],
        vertex_composition={0:2, 1:2, 2:1, 3:3, 4:1, 5:3})
    MC = ModuleDiagram(
        arrows=[Arrow(0,1), Arrow(0,2), Arrow(0,3), Arrow(1,4)],
        vertex_composition={0:3, 1:1, 2:2, 3:3, 4:1})



    A = AuslanderAlgebra(MB)
    print(A.vertices)
    print(A.vertices_explicit)
    print(A.all_edges_graph)
    print('edges data', A.all_edges_graph.edges(data=True))
    G = A.remove_redundant_edges()
    #
    # # Add nodes and edges
    # G.add_nodes_from(A.all_edges_graph.nodes(data=True))
    # G.add_edges_from(A.all_edges_graph.edges)

    # Draw the graph
    pos = nx.spring_layout(G)  # Compute layout
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=2000, font_size=12)

    # Draw edge labels (showing the "map" function)
    edge_labels = {(u, v, k): d for u, v, k, d in G.edges(keys=True, data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)

    # Show the plot
    plt.show()

    # print('vertices', A.vertices_explicit)
    # for i in range(9):
    #     print('\n paths at ', i, A.find_walks(i,10) )
    #print(len(A.inj_diagram()))
    #print(A.inj_diagram())
    # for i,j in enumerate(A.inj_diagram()):
    #     print('\n  inj at ', i)
    #     print('vertices', A.inj_diagram()[i].nodes(data=True))
    #     print('edges', A.inj_diagram()[i].edges(data=True))
    # for G in A.inj_diagram():
    #     print(type(G))
    #     nx.draw(G)
    #     plt.show()
    #print(A.inj_diagram()[3].edges(data=True))
    #print(A.find_walks(3))
