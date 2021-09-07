from rcsbsearch import Attr

# Setup the query
prot_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("Protein")
dna_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("DNA")
rna_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("RNA")
hyb_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("NA-hybrid")

q = prot_q
r = list(set(q()))

with open("./pdb_query/protein.txt", "w") as f:
    for i in r:
        print(i, file = f)
        
q = dna_q
r = list(set(q()))

with open("./pdb_query/dna.txt", "w") as f:
    for i in r:
        print(i, file = f)
        
q = rna_q
r = list(set(q()))

with open("./pdb_query/rna.txt", "w") as f:
    for i in r:
        print(i, file = f)
        
q = hyb_q
r = list(set(q()))

with open("./pdb_query/hybrid.txt", "w") as f:
    for i in r:
        print(i, file = f)