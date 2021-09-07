from glob import glob
from pdbfixer import pdbfixer
from simtk.openmm.app import PDBFile
from tqdm import tqdm

ch_list = glob("./pdb_ssDNAProtChains/*")

for chain_name in tqdm(ch_list):
    # print(chain_name)
    #Add missing atoms, remove water, ...
    try:
        fix = pdbfixer.PDBFixer(filename=chain_name)
        fix.findMissingResidues()
        fix.findMissingAtoms()
        fix.addMissingAtoms()
        fix.removeHeterogens(keepWater=False)
        fix.topology.createStandardBonds()
        PDBFile.writeFile(fix.topology, fix.positions, open(chain_name, 'w'), keepIds = True)
    except Exception as e:
        print("PDB fix err.", chain_name, e)
        continue
        
    # break