# Utils for RMG

## Atom-Atom Maps
RMG does not currently support output of labeled adjacency lists. However, this pull request
allows the modification of RMG source to generate them: https://github.com/ReactionMechanismGenerator/RMG-Py/pull/2018

Using the generated `chemkin/species_dictionary.txt` and `chemkin/reaction_adjacency_lists.txt` as input the
`atom_atom_maps.py` script generates a JSON with all reactions and their atom maps.

> python atom_atom_maps.py --species "/path/to/chemkin/species_dictionary.txt" --reactions "/path/to/chemkin/reaction_adjacency_lists.txt" --output "atom-atom-maps.json"

The output JSON is formatted as follows:
```json
{
  "species": {
    "CH4(1)": {
      "adjacencyMatrix": [[0, 1, 1, 1, 1],
                          [1, 0, 0, 0, 0],
                          [1, 0, 0, 0, 0],
                          [1, 0, 0, 0, 0],
                          [1, 0, 0, 0, 0]
      ],
      "edgeLabels": [[null, "S", "S", "S", "S"],
                     ["S", null, null, null, null],
                     ["S", null, null, null, null],
                     ["S", null, null, null, null],
                     ["S", null, null, null, null]
      ],
      "nodeLabels": ["C", "H", "H", "H", "H"]
    }
  },
  "reactions": [
    {
      "index": 12,
      "reaction": "OX(6) + C.[Pt](22) <=> HOX(8) + CH3X(7)",
      "reaction_family": "Surface_Abstraction_vdW",
      "mappings": [
        [
          [5, "OX(6)", 1, "HOX(8)", 2],
          [3, "CH4(1)", 1, "HOX(8)", 1],
          [1, "vacantX(3)", 0, "CH3X(7)", 4],
          [4, "OX(6)", 0, "HOX(8)", 0],
          [2, "CH4(1)", 0, "CH3X(7)", 0]
        ],
        [
          [5, "OX(6)", 1, "HOX(8)", 2],
          [1, "vacantX(3)", 0, "CH3X(7)", 4],
          [4, "OX(6)", 0, "HOX(8)", 0],
          [3, "CH4(1)", 4, "HOX(8)", 1],
          [2, "CH4(1)", 0, "CH3X(7)", 0]
        ],
        [
          [5, "OX(6)", 1, "HOX(8)", 2],
          [3, "CH4(1)", 2, "HOX(8)", 1],
          [1, "vacantX(3)", 0, "CH3X(7)", 4],
          [4, "OX(6)", 0, "HOX(8)", 0],
          [2, "CH4(1)", 0, "CH3X(7)", 0]
        ],
        [
          [5, "OX(6)", 1, "HOX(8)", 2],
          [1, "vacantX(3)", 0, "CH3X(7)", 4],
          [4, "OX(6)", 0, "HOX(8)", 0],
          [3, "CH4(1)", 3, "HOX(8)", 1],
          [2, "CH4(1)", 0, "CH3X(7)", 0]
        ]
      ]
    }
  ]
}
```

Reaction mappings cover all 