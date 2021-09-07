from utils.utils_bash import *
from vmd import molecule, atomsel
from glob import glob
import subprocess
from itertools import groupby, count
import os
import re
import itertools
from itertools import groupby, count

def get_str_bio_chains(mol_str, mol_bio):
    str_ch = sorted(uniq(atomsel("protein or nucleic", mol_str).chain))
    bio_ch = []
    for i in mol_bio:
        s = atomsel("protein or nucleic", i)
        bio_ch += uniq(s.chain)
    bio_ch = sorted(bio_ch)
    return str_ch, bio_ch
#     return str_ch == bio_ch

def count_frame(mol_str, mol_bio):
    str_numframe = molecule.numframes(mol_str)
    bio_numframe = [molecule.numframes(m) for m in mol_bio]
    return [str_numframe] + bio_numframe

def find_pair(file_path, path_to_find_pair = ""):
    bashCommand = f"{path_to_find_pair}find_pair {file_path}"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
#     print("o", "\n".join(str(output).split("\\n")))
    if error != None:
        print("e", error)
    return str(output), error

def get_conseutive_elements(data):
    c = count()
    val = list(list(g) for _, g in groupby(data, lambda x: x-next(c)))
    return val

def analyze_for_singleStrand(mol_str, min_na_len = 4):
    file_path = molecule.get_filenames(mol_str)[0]
    
    fp, fp_err = find_pair(file_path)
    
    nb_bp = int(re.search(" ([0-9]+)\ +# number", fp).group(1))
    try:
        nb_bp = int(re.search(" ([0-9]+)\ +# number", fp).group(1))
    except:
        return {}
    
    # Extract NA residues involved in pairing
    lst_bp = re.findall(">(.):\.*([\-0-9]+)_:", fp) + [(i[1], i[0]) for i in re.findall(":\.*([\-0-9]+)_:(.)<", fp)]
    
    # Extract NA residues in the structure
    sel = atomsel("nucleic and altloc 'A' '' ' ' and not name OP1 C5' OP2 C3' O5' C1' O4' O3' C4' C2' P", mol_str) ##!! ADD sidechain/base selection
    lst_na = sorted(uniq([(i[0],i[1]) for i in zip(sel.chain, map(str, sel.resid))])) 
#     print(lst_na)
#     print(lst_bp)

    # Find ssNA residues
    ssNA = set(lst_na)-set(lst_bp)
    ssNA_dict = {}
    for i in ssNA:
        try:
            ssNA_dict[i[0]].append(i[1])
        except:
            ssNA_dict[i[0]] = [i[1]]
    
    # Remove chains that are too short
    for i in list(ssNA_dict):
        j = ssNA_dict[i]
        ssNA_dict[i] = sorted([int(k) for k in j])
        if len(j) < min_na_len:
            del ssNA_dict[i]
    
    # Find consecutive residues
    for i in list(ssNA_dict):
        cons_elem = get_conseutive_elements(ssNA_dict[i])
        long_ce = []
        for e in cons_elem:
            if len(e) >= min_na_len:
                long_ce.append(e)
        if len(long_ce) != 0:
            ssNA_dict[i] = long_ce
        else:
            del ssNA_dict[i]
    
    return ssNA_dict

def analyse_main_str(mol_str):
    try:
        ss_analyze = analyze_for_singleStrand(mol_str)
    except :
        ss_analyze = {}
    if len(ss_analyze) != 0:
        add = [ss_analyze]
    else:
        add = [False]
    return add

def analyse_biol_ass(str_bio, add):
    if len(str_bio) != 0:
        for bio in str_bio:
            # Merge all BA frames
            bio_assembly = open(bio, "r").read().split("\n")
            bio_assembly = [line for line in bio_assembly if line.startswith("ATOM")]
            open("tmp_findpair.pdb", "w").write("\n".join(bio_assembly))
            mol_tmp = molecule.load("pdb", "tmp_findpair.pdb")
            # Analyse
            try:
                ss_analyze = analyze_for_singleStrand(mol_tmp)
            except :
                ss_analyze = {}
        
            # If the analysis does not report double strand,
            # Repport the structure file, not BA, with the ssNA
            if ss_analyze != {}:
                add.append(ss_analyze)
            else:
                add.append(False)
            os.remove("tmp_findpair.pdb")
#     else:
#         add.append({})
    return add

def ssDNA_str_biass(struc, biass):
    ssDNA_dict = {}
    # Take the intersection of the ssDNA residues in main structure and bioassemblies
    # If no bioassembly, take only from main structure
    if biass == None:
        ssDNA_dict = struc
        pass
    # if bioassembly, take the intersection
    else:
        # Merge Biass analysis
        biass = {k:v for element in biass if element != False for k,v in element.items()}
        # Find common ssDNA chains between str and BA
        common_chains = intersection = set(struc.keys()) & set(biass.keys())
        # Find common ssDNA residues between str and BA
        for c in common_chains:
            ssDNA_dict[c] = set(concat_list(struc[c])) & set(concat_list(biass[c]))
    
    return ssDNA_dict

def concat_list(list_of_list):
    return list(itertools.chain.from_iterable(list_of_list))

def interacting_chains(s, ssDNA_dict):
    mol = molecule.load("pdb", f"pdb_str/{s.lower()}.pdb")
    int_ssdna_chains = []
    int_ssdna_resid = []
    int_prot_chains = []
    
    for st, ssID in zip(ssDNA_dict.keys(), ssDNA_dict.values()):
        try:
            ssID = concat_list(ssID)
        except:
            pass
        ssID = sorted(ssID)
        ssID = "'" + "' '".join(map(str, ssID)) + "'"
        sel = f"protein and (within 5.0 of (chain '{st}' and resid {ssID}))"
        
        # Find protein's interacting residues and chains,
        s_prot = atomsel(sel, mol)
        # If none, continue
        if len(s_prot) == 0:
            continue
        # Find   ssDNA's Interacting residues and chains.
        sub_sel = []
        for c in uniq(s_prot.chain):
            sub_sel.append(f"chain '{c}' ")
            sub_sel[-1] = sub_sel[-1] + "and resid '"
            sub_sel[-1] = sub_sel[-1] + "' '".join(list(map(str, uniq(atomsel("chain " + c + " and " + sel).resid))))
            sub_sel[-1] = sub_sel[-1] + "'"
                    
        new_sel = f"(chain '{st}' and resid {ssID}) and (within 5.0 of ("
        new_sel = new_sel + "("+") or (".join(sub_sel)+")"
        new_sel = new_sel +"))"
        
        s_dna = atomsel(new_sel, mol)
                
        int_ssDNA = get_conseutive_elements(sorted(uniq(s_dna.resid)))
        int_ssDNA = [i for i in int_ssDNA if len(i) >= 4]
        if len(int_ssDNA) != 0:
            int_ssdna_resid.append([st, int_ssDNA])
            int_ssdna_chains.append(st)
            int_prot_chains += uniq(s_prot.chain)
    molecule.delete(mol)
    
    return uniq(int_ssdna_chains), uniq(int_prot_chains), int_ssdna_resid