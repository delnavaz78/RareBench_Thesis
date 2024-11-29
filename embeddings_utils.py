from PrimeKG_analysis import *
import sys
sys.path.append('..')
from PrimeKG_utils import *
import math

def dfs_to_graph_direct_acyclic(nodes, edges):
    graph = ig.Graph()

    graph.add_vertices(nodes['node_idx'])
    for attribute in nodes.columns:
        graph.vs[attribute] = nodes[attribute]
        
    edge_index = create_primekg_edge_index(edges)
    graph.add_edges([tuple(x) for x in edge_index])
    for attribute in edges.columns:
        graph.es[attribute] = edges[attribute]

    graph = graph.as_directed(mode='acyclic')

    return graph



def annotated_disease(G):
   
    # Identify all disease nodes
    disease_nodes = [node.index for node in G.vs if node['node_type'] == 'disease']
    
    # Calculate number of disease annotated to each node directly
    for node in G.vs:
        if node.index in disease_nodes:
            continue  # Skip calculating IC for disease nodes themselves
        
        node['annotated_disease'] = []
        # Get neighbors of the current node
        neighbors = G.neighbors(node.index, mode="all")
        
        # Identify neighbors that are diseases
        associated_diseases = [neighbor for neighbor in neighbors if neighbor in disease_nodes]
        for disease in associated_diseases:
            if disease not in node['annotated_disease']:
                node['annotated_disease'].append(disease)

    return G


def get_childeren(G, node_type):
    phenotypes = []
    phenotypes = [node for node in G.vs if node['node_type'] == node_type]

    for node in phenotypes:
        node["children"] = []

        # Populate the "children" attribute with the indices of successor nodes that are in phenotypes
        for succ in G.successors(node.index):  # Use node.index here
            if G.vs[succ] in phenotypes and G.vs[succ]['node_idx'] not in node['children']:  # Check if the successor is also in phenotypes
                node["children"].append(G.vs[succ]["node_idx"])

    return(G)

def get_parents(G, node_type):
    phenotypes = []
    phenotypes = [node for node in G.vs if node['node_type'] == node_type]

    for node in phenotypes:
        node["parents"] = []

        # Populate the "parents" attribute with the indices of successor nodes that are in phe
        for pre in G.predecessors(node.index):  # Use node.index here
            if G.vs[pre] in phenotypes and G.vs[pre]['node_idx'] not in node['parents']:  # Check if the successor is also in phenotypes
                node["parents"].append(G.vs[pre]["node_idx"])
                
    return(G)

from collections import deque

def calculate_nt(G, node_type):
    """
    Calculate n(t) values (disease counts including descendants) for each node
    in a directed acyclic graph, starting from leaf nodes.
    
    Parameters:
        G (igraph.Graph): A directed acyclic graph where nodes may be annotated with diseases.
        
    Returns:
        dict: A dictionary mapping each node_id to its n(t) value.
    """
    # initiate an array of phenotypes
    phenotypes = []
    phenotypes = [node for node in G.vs if node['node_type'] == node_type]

    # step1: find the diseases which are directly annotated to each phenotype node
    annotated_disease(G)

    # find children and parents of each phenotype node
    get_childeren(G, node_type)         
    get_parents(G, node_type)

    # Step 2: Identify leaf nodes with no child
    queue=deque(
        node for node in phenotypes
        if len(node['children']) == 0
    )
    
    # Initialize a visited set to avoid re-processing nodes
    visited = set(queue)

    # Step 3: Bottom-up aggregation of disease counts
    while queue:
        node = queue.pop()

        # Process each parent of the current node
        for parent_idx in node['parents']:
            parent = G.vs[parent_idx]  # Retrieve the parent node object using the index

            if parent not in queue:
                queue.appendleft(parent)            
            
            if parent not in visited:
                for disease in node["annotated_disease"]:
                    if disease not in parent['annotated_disease']:
                        parent['annotated_disease'].append(disease)
            
            # Add parent to the queue if all its children have been processed
            all_children_visited = all(child in visited for child in parent['children'])
            if all_children_visited and parent not in visited:
                visited.appendleft(parent)
    
    # Step 4: Map results to node_id for clarity
    nt_values = {node["node_id"]: len(node["annotated_disease"]) for node in phenotypes}
    
    return nt_values

def calculate_ic_values(G, nt_values):
    """
    Calculate the Information Content (IC) for each node associated with a disease in the graph.
    
    Parameters:
        G (igraph.Graph): A graph with nodes of various types and relationships.
        
    Returns:
        dict: A dictionary mapping each node_id to its IC value.
    """
    
    disease = []
    disease = [node for node in G.vs if node['node_type'] == "disease"]
    N = len(disease)
    
    ic_values = {}
    
    # Calculate IC for each node connected to a disease
    for key,value in nt_values.items():
        node_id = key 
        n_t = value

        if n_t > 0:
            ic_values[node_id] = -math.log(n_t / N)
        else:
            ic_values[node_id] = 0.0  # If no associated diseases, IC is set to 0
    
    return ic_values

def get_weights(IC_values, l, e):
    weights ={}
    for key, value in IC_values.items():
        weight = l * value + e
        weights[key] = weight

    return weights