import sys
sys.path.append('..')
import math
import networkx as nx
from collections import deque

def add_node_types(edges, nodes):
    edges_with_node_type = edges.merge(nodes[['name', 'type']], left_on='subject', right_on='name', how='left')
    edges_with_node_type = edges_with_node_type.rename(columns={'type': 'x_type'}).drop(columns=['name'])
    edges_with_node_type = edges_with_node_type.merge(nodes[['name', 'type']], left_on='object', right_on='name', how='left')
    edges_with_node_type = edges_with_node_type.rename(columns={'type': 'y_type'}).drop(columns=['name'])
    return edges_with_node_type



def annotated_disease(G):
    """
    Annotates each node in a NetworkX graph with a list of directly associated disease nodes.

    Parameters:
        G (nx.DiGraph): A NetworkX graph with nodes having a 'node_type' attribute.

    Returns:
        G (nx.DiGraph): The graph with an 'annotated_disease' attribute added to each node.
    """

    for node in G.nodes():
        if G.nodes[node].get("node_type") == "Disease":
            continue
        
        nx.set_node_attributes(G, {node: []}, "annotated_disease")
        successors = list(G.successors(node))

        associated_diseases = [
            succ for succ in successors 
            if G.nodes[succ].get("node_type") == "Disease"
        ]
        
        # Annotate the current node with its associated diseases
        G.nodes[node]["annotated_disease"] = associated_diseases

    return G



def get_children(G, node_type):
    """
    Annotates nodes of a specified type with their children in a NetworkX graph.

    Parameters:
        G (nx.DiGraph): A NetworkX directed graph.
        node_type (str): The type of nodes to process (e.g., 'Phenotype').

    Returns:
        G (nx.DiGraph): The graph with a 'children' attribute added to nodes.
    """
    # Get all nodes of the specified type
    phenotypes = [node for node, data in G.nodes(data=True) if data.get("node_type") == node_type]
    
    for node in phenotypes:
        # Initialize the "children" attribute
        G.nodes[node]["children"] = []
        
        # Populate the "children" attribute with predecessors of the current node that are also phenotypes
        for pre in G.predecessors(node):
            if G.nodes[pre].get("node_type") == node_type and pre not in G.nodes[node]["children"]:
                G.nodes[node]["children"].append(pre)
    
    return G


def get_parents(G, node_type):
    """
    Annotates nodes of a specified type with their parents in a NetworkX graph.

    Parameters:
        G (nx.DiGraph): A NetworkX directed graph.
        node_type (str): The type of nodes to process (e.g., 'Phenotype').

    Returns:
        G (nx.DiGraph): The graph with a 'parents' attribute added to nodes.
    """
    # Get all nodes of the specified type
    phenotypes = [node for node, data in G.nodes(data=True) if data.get("node_type") == node_type]
    
    for node in phenotypes:
        # Initialize the "parents" attribute
        G.nodes[node]["parents"] = []
        
        # Populate the "parents" attribute with successors of the current node that are also phenotypes
        for succ in G.successors(node):
            if G.nodes[succ].get("node_type") == node_type and succ not in G.nodes[node]["parents"]:
                G.nodes[node]["parents"].append(succ)
    
    return G


def calculate_nt(G, node_type):
    """
    Calculate n(t) values (disease counts including descendants) for each node
    in a directed acyclic graph, starting from leaf nodes.

    Parameters:
        G (networkx.DiGraph): A directed acyclic graph where nodes may be annotated with diseases.
        node_type (str): The type of nodes to process (e.g., 'Phenotype').

    Returns:
        dict: A dictionary mapping each node ID to its n(t) value.
    """
    # Step 1: Identify phenotype nodes
    phenotypes = [node for node, attrs in G.nodes(data=True) if attrs.get('node_type') == node_type]


    # Step 2: Identify leaf nodes with no children
    queue = deque(
        node for node in phenotypes if len(G.nodes[node]["children"]) == 0
    )

    # step3: Initialize a visited set to avoid re-processing nodes
    visited = set(queue)

    # Step 4: Bottom-up aggregation of disease counts
    while queue:
        node = queue.pop()

        # Process each parent of the current node
        for parent in G.nodes[node]["parents"]:
            if parent not in visited:
                for disease in G.nodes[node]["annotated_disease"]:
                    if disease not in G.nodes[parent]["annotated_disease"]:
                        G.nodes[parent]["annotated_disease"].append(disease)

                # Add parent to the queue if all its children have been processed
                all_children_visited = all(
                    child in visited for child in G.nodes[parent]["children"]
                )
                if all_children_visited:
                    queue.appendleft(parent)
                    visited.add(parent)

    # Step 5: Map results to node_id for clarity
    nt_values = {
        node: len(G.nodes[node]["annotated_disease"]) for node in phenotypes
    }

    return nt_values


def calculate_ic_values(G, nt_values):
    """
    Calculate the Information Content (IC) for each node associated with a disease in the graph.
    
    Parameters:
        G (networkx.DiGraph): A directed graph with nodes of various types and relationships.
        nt_values (dict): A dictionary mapping each node_id to its n(t) value.
        
    Returns:
        dict: A dictionary mapping each node_id to its IC value.
    """
    
    # Identify disease nodes
    disease_nodes = [node for node, attrs in G.nodes(data=True) if attrs.get('node_type') == "Disease"]
    N = len(disease_nodes)  # Total number of disease nodes

    ic_values = {}

    # Calculate IC for each node connected to a disease
    for node_id, n_t in nt_values.items():
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