import argparse
import copy
import json
import re
from typing import List, Dict, Union
from yaml import load, Loader


def parse_atom_definition(s: str) -> Dict[str, Union[int, str, List, None]]:
    # https://reactionmechanismgenerator.github.io/RMG-Py/reference/molecule/adjlist.html#rmgpy-molecule-adjlist
    # TODO: site, morphology
    pattern = r'^(\d+)(?:\s+\*(\d*)\b)?\s+([A-Z][a-z]?)\s+(u\d+)(?:\s+(p\d+))?(?:\s+(c[-+]?\d+))?(?:\s+(\{\d+,\s*.*\}(?:\s+\{\d+,\s*.*\})*))?$'
    match = re.match(pattern, s.strip())
    if not match:
        raise ValueError("Invalid atom definition: " + s)
    groups = match.groups()
    bonds = [tuple(x.strip() for x in bond.strip('{}').split(',')) for bond in re.split(r'}\s+\{', groups[6])] if \
    groups[6] is not None else []
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


def atom_definition_to_string(atom) -> str:
    data = [str(atom['number'])]
    if atom['label'] is not None:
        data.append('*' + str(atom['label']))
    data.append(atom['element'])
    data.append('u' + str(atom['unpaired']))
    if atom['pairs'] is not None:
        data.append('p' + str(atom['pairs']))
    if atom['charge'] is not None:
        data.append('c' + str(atom['charge']))
    for bond in atom['bonds']:
        data.append('{%s,%s}' % (bond[0], bond[1]))
    return ' '.join(data)


def main():
    parser = argparse.ArgumentParser(description="Process a reaction adjacency list file and extract atom maps.")
    parser.add_argument('--input', required=True, help='Path to the input chemkin/reaction_adjacency_lists.txt')
    parser.add_argument('--output', default='atom-atom-maps.json', help='Path to the output atom-atom-maps.json')
    args = parser.parse_args()
    output = []
    with open(args.reactions, 'r') as f:
        data = load(f, Loader=Loader)
        for reaction in data['reactions']:
            if 'reactant' in reaction and 'product' in reaction:
                reactants_graph = [
                    parse_atom_definition(x.strip())
                    for x in reaction['reactant'].split('\n')
                    if x.strip() != ''
                ]
                products_graph = [
                    parse_atom_definition(x.strip())
                    for x in reaction['product'].split('\n')
                    if x.strip() != ''
                ]

                reactant_id_index_map = {
                    atom['label']: i
                    for i, atom in enumerate(reactants_graph)
                    if atom['label'] is not None
                }
                reactant_number_id_map = {
                    atom['number']: atom['label']
                    for atom in reactants_graph
                    if atom['label'] is not None
                }
                product_id_index_map = {
                    atom['label']: i
                    for i, atom in enumerate(products_graph)
                    if atom['label'] is not None
                }
                product_number_id_map = {
                    atom['number']: atom['label']
                    for atom in products_graph
                    if atom['label'] is not None
                }
                mapped_products_graph = copy.deepcopy(reactants_graph)
                for _id, index in reactant_id_index_map.items():
                    reactant_atom = mapped_products_graph[index]
                    product_atom = products_graph[product_id_index_map[_id]]
                    reactant_atom['unpaired'] = product_atom['unpaired']
                    reactant_atom['pairs'] = product_atom['pairs']
                    reactant_atom['charge'] = product_atom['charge']
                    bonds = [b for b in reactant_atom['bonds'] if b[0] not in reactant_number_id_map]
                    for bond in product_atom['bonds']:
                        if bond[0] in product_number_id_map:
                            bond_target_id = product_number_id_map[bond[0]]
                            bond_mapped_number = mapped_products_graph[reactant_id_index_map[bond_target_id]]['number']
                            bonds.append((bond_mapped_number, bond[1]))
                    reactant_atom['bonds'] = product_atom['bonds']
                max_id = max(reactant_id_index_map.keys())
                for i, atom in enumerate(mapped_products_graph):
                    if atom['label'] is None:
                        max_id += 1
                        reactants_graph[i]['label'] = max_id
                        atom['label'] = max_id
                output.append({
                    'index': reaction['index'],
                    'reaction': reaction['reaction'],
                    'reaction_family': reaction['reaction_family'],
                    'reactant': [atom_definition_to_string(x) for x in reactants_graph],
                    'product': [atom_definition_to_string(x) for x in mapped_products_graph]
                })

    with open(args.output, 'w', encoding='utf-8', newline='') as f:
        json.dump(output, f, indent=2)


if __name__ == "__main__":
    main()
