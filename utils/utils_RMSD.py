from vmd import molecule, atomsel
from Bio.SeqUtils import seq1

def get_chain_1l_sequence(mol_id, chain):
    sel = atomsel(f"name CA and chain {chain}", mol_id)
    seq = sel.resname
    return seq1("".join(seq)), sel

####    ####    ####    ####    ####    ####    ####
def full_rmsd(pdb_1, chain_1, pdb_2, chain_2):
    """
    """
    mol_id_1 = molecule.load("pdb", f"./pdb/{pdb_1}.pdb")
    mol_id_2 = molecule.load("pdb", f"./pdb/{pdb_2}.pdb")
    
    seq_1, sel_1 = get_chain_1l_sequence(mol_id_1, chain_1)
    seq_2, sel_2 = get_chain_1l_sequence(mol_id_2, chain_2)
    
    if seq_1 == seq_2 and len(seq_1) == len(seq_2):
        fit_m = sel_1.fit(sel_2)
        sel_1.move(fit_m)
        rmsd = sel_1.rmsd(sel_2)
    else:
        rmsd = None
        
    molecule.delete(mol_id_1)
    molecule.delete(mol_id_2)
    
    return(rmsd)
####    ####    ####    ####    ####    ####    ####

