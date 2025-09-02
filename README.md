# Utils for RMG

## Atom-Atom Maps
RMG does not currently support output of labeled adjacency lists. However, this pull request
allows the modification of RMG source to generate them: https://github.com/ReactionMechanismGenerator/RMG-Py/pull/2018

Using the generated `chemkin/reaction_adjacency_lists.txt` as input the
`atom_atom_maps.py` script generates a JSON with all reactions and their atom maps.

> python atom_atom_maps.py --input "/path/to/chemkin/reaction_adjacency_lists.txt" --output "atom-atom-maps.json"

The output JSON is formatted as follows:
```json
[
  {
    "index": 78,
    "reaction": "H2(4) + CO2(11) <=> O=CO(56)",
    "reaction_family": "1,3_Insertion_CO2",
    "reactant": [
      "1 *2 O u0 p2 c0 {3,D}",
      "2 *5 O u0 p2 c0 {3,D}",
      "3 *1 C u0 p0 c0 {1,D} {2,D}",
      "4 *4 H u0 p0 c0 {5,S}",
      "5 *3 H u0 p0 c0 {4,S}"
    ],
    "product": [
      "1 *2 O u0 p2 c0 {3,S} {4,S}",
      "2 *5 O u0 p2 c0 {3,D}",
      "3 *1 C u0 p0 c0 {1,S} {2,D} {5,S}",
      "4 *4 H u0 p0 c0 {1,S}",
      "5 *3 H u0 p0 c0 {3,S}"
    ]
  }
]
```

For information regarding the adjacency list format see: https://reactionmechanismgenerator.github.io/RMG-Py/reference/molecule/adjlist.html#rmgpy-molecule-adjlist