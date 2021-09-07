from utils.utils_rcsb import *
from utils.utils_bash import *
from utils.utils_benchmark import *
from tqdm import tqdm
import re

from rcsbsearch import Attr

from vmd import molecule, atomsel
from pdbfixer import pdbfixer
from simtk.openmm.app import PDBFile

from glob import glob
import subprocess
from itertools import groupby, count, chain, combinations, permutations
import os, json
import requests, gzip
import numpy as np

from os.path import getmtime
import time

import subprocess
from multiprocessing.pool import ThreadPool

def get_sels(aln1, aln2, mol_1, mol_2, aln_file):
    cont = False
    sel_1, sel_2 = find_all(aln1, aln2, mol_1, mol_2)
    
    if (sel_1 is None) or (sel_2 is None):
        return None, None, True
    if len(sel_1) != len(sel_2):
        return None, None, True
    return sel_1, sel_2, False
                    
    # If failed, try end-to end alignment
    if (sel_1 is None) or (sel_2 is None):
        sel_1, sel_2 = find_ext(aln1, aln2, mol_1, mol_2)
    elif len(sel_1) != len(sel_2):
        sel_1, sel_2 = find_ext(aln1, aln2, mol_1, mol_2)
    if (sel_1 is None) or (sel_2 is None):
        print("", end = "", file = open(aln_file[0].replace("results_fatcat", "seq_align"), "w"))
        cont = True
    if len(sel_1) != len(sel_2):
        print(sel_1, len(sel_1) , i1)
        print(sel_2, len(sel_2) , i2)
        cont = True
    return sel_1, sel_2, cont

def get_sequence_alignment(fatcat_aln : str):
    aln = [j[14:] for i,j in enumerate(fatcat_aln.split("\n")[:-3]) if (len(j) > 0 and i >= 4)]
    ch1 = "".join([i for i in aln[0::3]])
    info = "".join([i for i in aln[1::3]])
    ch2 = "".join([i for i in aln[2::3]])
            
    all_c1, all_c2 = "",""
    prev_info, prev_c1, prev_c2 = "", "", ""
    print_info, print_c1, print_c2 = "", "", ""
    for c1, inf, c2 in zip(ch1,info,ch2):
        if inf == "1":
            print_info, print_c1, print_c2 = print_info+"1", print_c1+c1, print_c2+c2
                
        if prev_info == "1" and inf != "1":
            if len(print_c1) > 10 and print_c1 == print_c2:
                all_c1 = all_c1 + " " + print_c1
                all_c2 = all_c2 + " " + print_c2
            print_info, print_c1, print_c2 = "", "", ""
        prev_info, prev_c1, prev_c2 = inf, c1, c2
    
    if len(print_c1) > 5 and print_c1 == print_c2:
        all_c1 = all_c1 + " " + print_c1
        all_c2 = all_c2 + " " + print_c2
    
    return all_c1.replace("-", " ")[1:], all_c2.replace("-", " ")[1:]

def find_all(aln1, aln2, mol_1, mol_2):
    sel_1 = atomsel(f"name CA and sequence {aln1}", mol_1)
    sel_2 = atomsel(f"name CA and sequence {aln2}", mol_2)
    return sel_1, sel_2

def find_ext(aln1, aln2, mol_1, mol_2):
    
    s = 0
    loop_while = True
    while loop_while:
        res1 = [aln1[s:i] for i in range(s+1, len(aln1))]
        res2 = [aln2[s:i] for i in range(s+1, len(aln2))]
        
        prev = "",""
        for seq1, seq2 in zip(res1, res2):
            try:
                sel_1 = atomsel(f"name CA and sequence {seq1}", mol_1)
                sel_2 = atomsel(f"name CA and sequence {seq2}", mol_2)
            except:
                continue
            # Extend until no match
            if len(sel_1) == len(sel_2):
                prev = seq1, seq2
            # If no match, check if previous match is long enought
            elif len(prev[0]) > len(aln1)/2:
                loop_while = False
                break
            # Else check next starting point for sequence
            else:
                break
        
        s+=1
        if s > len(aln1) or len(prev) == 0:
            sel_1, sel_2 = None, None
            return sel_1, sel_2

    sel_1 = atomsel(f"name CA and sequence {prev[0]}", mol_1)
    sel_2 = atomsel(f"name CA and sequence {prev[1]}", mol_2)
    return sel_1, sel_2


df_bound_unbound = pd.read_csv("./csv_tables/bound_unbound.preFATCAT.csv", index_col = 0)

seq_done = set(glob("./seq_align/*"))
for i,(cl, b,ub, _,_,_,) in tqdm(df_bound_unbound.iterrows(), total = len(df_bound_unbound)):
    if i < 150:
        continue
    b = b.split(" ")
    ub = ub.split(" ")
    
    # Limit disk's IO
    ldict_ssDNA = {f"{j[:4].lower()}_{j[5:].upper()}" : [glob(f"./pdb_ssDNAProtChains/{j[:4].lower()}_{j[5:].upper()}*"), [molecule.load("pdb", k) for k in glob(f"./pdb_ssDNAProtChains/{j[:4].lower()}_{j[5:].upper()}*")]] for j in b}
    ldict_prots = {f"{j[:4]}_{j[5:].upper()}" : [glob(f"./pdb_chains/{j[:4]}_{j[5:]}*"), [molecule.load("pdb", k) for k in glob(f"./pdb_chains/{j[:4]}_{j[5:]}*")]] for j in map(str.lower, ub)}
    ldict = dict(ldict_ssDNA, **ldict_prots)
    
    for i1, i2 in tqdm(list(combinations(ldict.items(), 2))[:]):
#         if i1[0] != "1rcn_E" or i2[0] != "3jw1_A":
#             continue

#         if len(i1[1][0]) > 1 or len(i2[1][0]) > 1:
#             print(i1)
#             print(i2)
#             print("----")
#         continue
        
        # Do not repeat previous alignment
        if ((f"./seq_align/{i1[0].lower()}_{i2[0].lower()}.txt" in seq_done)
         or (f"./seq_align/{i2[0].lower()}_{i1[0].lower()}.txt" in seq_done)):
            continue
        
        aln_file = glob(f"./results_fatcat/{i1[0].lower()}_{i2[0].lower()}*") + glob(f"./results_fatcat/{i2[0].lower()}_{i1[0].lower()}*")

        #Check the number of returned alignment files
        if len(aln_file) > 1 or len(aln_file) == 0:
#             print("#err", i1, i2, aln_file)
            continue
            
        
        aln = open(aln_file[0]).read()
        if "error" in aln:
            continue
        
        aln1, aln2 = get_sequence_alignment(aln)
        
        # Search sor sequence alignement
        try:
            for str_1, mol_1 in zip(*i1[1]):
                for str_2, mol_2 in zip(*i2[1]):
                    # Try full seq alignement
                    sel_1, sel_2, cont = get_sels(aln1, aln2, mol_1, mol_2, aln_file)
                    if cont:
                        continue
                        
                    try:
                        # Align
                         # r0 = sel_1.rmsd(sel_2)
                        fit_m = sel_2.fit(sel_1)
                        atomsel("all", mol_2).move(fit_m)
                        molecule.write(mol_2, "pdb", str_2)
                         # r1 = sel_1.rmsd(sel_2)

                        # save Align result for later use
                        with open(aln_file[0].replace("results_fatcat", "seq_align"), "w") as f:
                            print(str_1, sel_1, file = f)
                            print(str_2, sel_2, file = f)
                            sel_1_resid = atomsel(str(sel_1), mol_1)
                            sel_2_resid = atomsel(str(sel_2), mol_2)
                            print(str_1, "name CA and resid '" + "' '".join(map(str, sel_1_resid.resid)) + "'", file = f)
                            print(str_2, "name CA and resid '" + "' '".join(map(str, sel_2_resid.resid)) + "'", file = f)
                    except Exception as e:
                        print("--", e)
                        print(i1, i2, aln_file)
                        print(len(sel_1), len(sel_2))
                        print("", end = "", file = open(aln_file[0].replace("results_fatcat", "seq_align"), "w"))
                        pass
                    
        except Exception as e:
            print(i1, i2, aln_file)
            print("**", e)
            print("", end = "", file = open(aln_file[0].replace("results_fatcat", "seq_align"), "w"))
            continue

#         break
#     break
for mol in molecule.listall():
    molecule.delete(mol)