from typing import List, Dict, Union, TypedDict, Optional
import argparse
import re
import json
from networkx.algorithms.isomorphism.isomorphvf2 import GraphMatcher
from yaml import load, Loader
import networkx as nx
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

class MoleculeGraph(TypedDict):
    adjacencyMatrix: List[List[int]]
    edgeLabels: List[List[Optional[str]]]
    nodeLabels: List[str]
    nodeIds: List[Optional[str]]


def parse_atom_definition(s: str) -> Dict[str, Union[int, str, List, None]]:
    # https://reactionmechanismgenerator.github.io/RMG-Py/reference/molecule/adjlist.html#rmgpy-molecule-adjlist
    # TODO: site, morphology
    pattern = r'^(\d+)(?:\s+\*(\d*)\b)?\s+([A-Z][a-z]?)\s+(u\d+)(?:\s+(p\d+))?(?:\s+(c[-+]?\d+))?(?:\s+(\{\d+,\s*.*\}(?:\s+\{\d+,\s*.*\})*))?$'
    match = re.match(pattern, s.strip())
    if not match:
        raise ValueError("Invalid atom definition: " + s)
    groups = match.groups()
    bonds = [tuple(x.strip() for x in bond.strip('{}').split(',')) for bond in re.split(r'}\s+\{', groups[6])] if groups[6] is not None else []
    bonds = [(int(x[0]), x[1]) for x in bonds]
    return {
        'number': int(groups[0]),
        'label': int(groups[1]) if groups[1] is not None else None,
        'element': groups[2],
        'unpaired': int(groups[3][1::]),
        'pairs': int(groups[4][1::]) if groups[4] is not None else None,
        'charge': int(groups[5][1::]) if groups[5] is not None else None,
        'bonds': bonds
    }


def adjlist_to_graph(adjlist: List[str]) -> MoleculeGraph:
    atoms = [parse_atom_definition(x) for x in adjlist]
    id_map = {int(x['number']): i for i, x in enumerate(atoms)}
    g: MoleculeGraph = {
        'adjacencyMatrix': [[0] * len(atoms) for _ in atoms],
        'edgeLabels': [[None] * len(atoms) for _ in atoms],
        'nodeLabels': [x['element'] for x in atoms],
        'nodeIds': [x['label'] for x in atoms]
    }
    for i, atom in enumerate(atoms):
        for bond in atom['bonds']:
            g['adjacencyMatrix'][i][id_map[bond[0]]] = 1
            g['edgeLabels'][i][id_map[bond[0]]] = bond[1]
    return g


def get_subgraph(graph: MoleculeGraph, indices: List[int]) -> MoleculeGraph:
    return {
        'adjacencyMatrix': [[graph['adjacencyMatrix'][i][j] for j in indices] for i in indices],
        'edgeLabels': [[graph['edgeLabels'][i][j] for j in indices] for i in indices],
        'nodeLabels': [graph['nodeLabels'][i] for i in indices],
        'nodeIds': [graph['nodeIds'][i] for i in indices]
    }


def get_component_subgraphs(graph: MoleculeGraph) -> List[MoleculeGraph]:
    _, component_labels = connected_components(csgraph=csr_array(graph['adjacencyMatrix']), directed=False,
                                               return_labels=True)
    component_indices = {k: [i for i, l in enumerate(component_labels) if l == k] for k in set(component_labels)}
    return [get_subgraph(graph, indices) for indices in component_indices.values()]


def get_nx_graph(graph: MoleculeGraph) -> nx.Graph:
    g = nx.Graph()
    for i, (label, node_id) in enumerate(zip(graph['nodeLabels'], graph['nodeIds'])):
        g.add_node(i, label=label, id=node_id)
    for i in range(len(graph['adjacencyMatrix'])):
        for j in range(i + 1, len(graph['adjacencyMatrix'])):
            if graph['adjacencyMatrix'][i][j] == 1:
                g.add_edge(i, j, label=graph['edgeLabels'][i][j])
    return g


def get_label_filtered_isomorphisms(g1: nx.Graph, g2: nx.Graph) -> List[Dict]:
    """
    Returns all isomorphisms between G1 and G2 that preserve node and edge labels.
    Assumes labels are stored under the 'label' attribute.
    """
    matcher = GraphMatcher(g1, g2)
    valid_isomorphisms = []
    for iso in matcher.isomorphisms_iter():
        node_match = all(
            g1.nodes[n1].get('label') == g2.nodes[n2].get('label')
            for n1, n2 in iso.items()
        )
        if not node_match:
            continue
        edge_match = True
        for u, v in g1.edges():
            u_mapped, v_mapped = iso[u], iso[v]
            if not g2.has_edge(u_mapped, v_mapped):
                edge_match = False
                break
            label1 = g1[u][v].get('label')
            label2 = g2[u_mapped][v_mapped].get('label')
            if label1 != label2:
                edge_match = False
                break
        if edge_match:
            valid_isomorphisms.append(iso)
    return valid_isomorphisms


def main():
    parser = argparse.ArgumentParser(description="Process two input files and write to an output file.")
    parser.add_argument('--species', required=True, help='Path to the input chemkin/species_dictionary.txt')
    parser.add_argument('--reactions', required=True, help='Path to the input chemkin/reaction_adjacency_lists.txt')
    parser.add_argument('--output', default='atom-atom-maps.json', help='Path to the output atom-atom-maps.json')
    args = parser.parse_args()
    # Final output object to be persisted as JSON
    output = {
        'species': {},
        'reactions': []
    }
    # Dictionary of all species with their name as key and graph as value
    species = {}
    with open(args.species, 'r') as f:
        lines = f.readlines()
        last_name = None
        last_rows = []
        for line in lines:
            if line.strip() == '':
                if last_name is not None:
                    g = adjlist_to_graph(last_rows)
                    species[last_name] = (g, get_nx_graph(g))
                    last_name = None
                    last_rows.clear()
                continue
            if last_name is None:
                last_name = line.strip()
            elif not line.startswith('multiplicity'):
                last_rows.append(line.strip())

    # Fast lookup for isomorphism checks to only check species of same size
    atom_count_species_names_map = {}
    for name in species:
        output['species'][name] = species[name][0]
        atom_count = len(species[name][0]['nodeLabels'])
        if atom_count not in atom_count_species_names_map:
            atom_count_species_names_map[atom_count] = [name]
        else:
            atom_count_species_names_map[atom_count].append(name)

    with open(args.reactions, 'r') as f:
        data = load(f, Loader=Loader)
        for reaction in data['reactions']:
            # reaction_parts = reaction['reaction'].split(' <=> ')
            # reactants = [x.strip() for x in reaction_parts[0].split(' + ')]
            # products = [x.strip() for x in reaction_parts[1].split(' + ')]
            # print(reaction['index'], reaction['reaction_family'], reactants, products)
            if 'reactant' in reaction and 'product' in reaction:
                # Create a graph from the reactants/products adjacency list
                reactants_graph = adjlist_to_graph([
                    x.strip() for x in reaction['reactant'].split('\n') if x.strip() != ''
                ])
                products_graph = adjlist_to_graph(
                    [x.strip() for x in reaction['product'].split('\n') if x.strip() != ''
                ])
                # Find connected components in the reactants/products graph and extract all reactant/product subgraphs
                reactants_subgraphs = get_component_subgraphs(reactants_graph)
                products_subgraphs = get_component_subgraphs(products_graph)
                # Convert the subgraphs to nx graphs for isomorphism checking
                reactants_nx_subgraphs = [get_nx_graph(g) for g in reactants_subgraphs]
                products_nx_subgraphs = [get_nx_graph(g) for g in products_subgraphs]
                # Match all reactant/product subgraphs to the species list via isomorphism checking
                reactant_species_mappings = [[] for _ in reactants_subgraphs]
                reactant_species_names = [None] * len(reactants_subgraphs)
                product_species_mappings = [[] for _ in products_subgraphs]
                product_species_names = [None] * len(products_subgraphs)
                for i, reactant in enumerate(reactants_nx_subgraphs):
                    atom_count = len(reactants_subgraphs[i]['nodeLabels'])
                    for name in atom_count_species_names_map[atom_count]:
                        isomorphisms = get_label_filtered_isomorphisms(reactant, species[name][1])
                        if len(isomorphisms) > 0:
                            if reactant_species_names[i] is not None and reactant_species_names[i] != name:
                                print('[WARN] Reactant matches multiple species', reactant_species_names[i], name,
                                      reactants_subgraphs[i])
                                continue
                            reactant_species_names[i] = name
                            reactant_species_mappings[i].extend(isomorphisms)
                for i, product in enumerate(products_nx_subgraphs):
                    atom_count = len(products_subgraphs[i]['nodeLabels'])
                    for name in atom_count_species_names_map[atom_count]:
                        isomorphisms = get_label_filtered_isomorphisms(product, species[name][1])
                        if len(isomorphisms) > 0:
                            if product_species_names[i] is not None and product_species_names[i] != name:
                                print('[WARN] Product matches multiple species', product_species_names[i], name,
                                      products_subgraphs[i])
                                continue
                            product_species_names[i] = name
                            product_species_mappings[i].extend(isomorphisms)
                # Build all reactant/product isomorphism combinations
                def recurse_species_mappings(species_mappings, index, path, results):
                    if index == len(species_mappings):
                        results.append(path)
                        return results
                    for i in range(len(species_mappings[index])):
                        recurse_species_mappings(species_mappings, index + 1, path + [i], results)
                    return results
                reactant_combinations = recurse_species_mappings(reactant_species_mappings, 0, [], [])
                product_combinations = recurse_species_mappings(product_species_mappings, 0, [], [])
                # Build all final mappings in the form (label, source_name, source_index, target_name, target_index)
                mappings = set()
                for reactant_combination in reactant_combinations:
                    for product_combination in product_combinations:
                        mapping = []
                        for i, isomorphism_index in enumerate(reactant_combination):
                            isomorphism = reactant_species_mappings[i][isomorphism_index]
                            for j, _id in enumerate(reactants_subgraphs[i]['nodeIds']):
                                if _id is not None:
                                    mapping.append([_id, reactant_species_names[i], isomorphism[j], None, None])
                        for i, isomorphism_index in enumerate(product_combination):
                            isomorphism = product_species_mappings[i][isomorphism_index]
                            for j, _id in enumerate(products_subgraphs[i]['nodeIds']):
                                if _id is not None:
                                    _id_mapping = [x for x in mapping if x[0] == _id][0]
                                    _id_mapping[3] = product_species_names[i]
                                    _id_mapping[4] = isomorphism[j]
                        mappings.add(frozenset([tuple(x) for x in mapping]))
                output['reactions'].append({
                    'index': reaction['index'],
                    'reaction': reaction['reaction'],
                    'reaction_family': reaction['reaction_family'],
                    'mappings': [list(x) for x in mappings]
                })

    with open(args.output, 'w', encoding='utf-8', newline='') as f:
        json.dump(output, f)

if __name__ == "__main__":
    main()