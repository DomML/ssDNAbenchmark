# python testBioPythonAlign.py --bound $(cat b_ub_pairs.txt | cut -d" " -f1) --unbound $(cat b_ub_pairs.txt | cut -d" " -f2) --prefix "./pdb_dl/" --suffix "pdb" > ssDNA_RMSD_sasa.csv 2> /dev/null

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import seq1

from vmd import molecule, atomsel, evaltcl
import argparse, sys
from tqdm import tqdm

import traceback, inspect

# Typing hint
from typing import List
ResidList = List[int]

# CLI arguments parsing
parser = argparse.ArgumentParser(description='Align bound and unbound proteins')
parser.add_argument("--prefix", type=str, help="Folder name of pdb")
parser.add_argument("--suffix", type=str, help="pdb extension")
parser.add_argument("--bound", type=str, help="List of bound pdb file", nargs='+')
parser.add_argument("--unbound", type=str, help="List of unbound pdb file", nargs='+')
parser.add_argument("--bound_chains", type=str, help="List of bound pdb file", nargs='+')
parser.add_argument("--unbound_chains", type=str, help="List of unbound pdb file", nargs='+')
parser.add_argument("-t", help="Display progress",default=False, action='store_true')
parser.add_argument("-v", help="Display alignment",default=False, action='store_true')
parser.add_argument("-e", help="Output errors in stdout",default=False, action='store_true')
parser.add_argument("-p", help="Display protein's VMD commands",default=False, action='store_true')
parser.add_argument("-s", help="Display sasa's VMD commands",default=False, action='store_true')
parser.add_argument("-d", help="Display DNA's VMD commands",default=False, action='store_true')
args = parser.parse_args()

def no_tqdm(mylist):
    return mylist

def main():

    if (not args.t): tqdm = no_tqdm

    # Check if feed with the same number of bound and unbound strucutres
    if (len(args.bound)   != len(args.unbound)
     or len(args.bound)   != len(args.bound_chains) 
     or len(args.unbound) != len(args.unbound_chains)):
        print("len(bound) != len(unbound)")
        sys.exit()

    print("#b", "b_chain", "ub" , "ub_chain", "allProtRmsd", "DnaProtRmsd", "b_chain/ub_chain", "allProtRmsd", "DnaProtRmsd", "b_Sasas(0)", "b_Sasas(1)", "ub_Sasas(0)", "ub_Sasas(1)", "b_len", "ubl_len", sep=",")

    # Loop over bound/unbound pairs
    for b, ub, bc, ubc in tqdm(list(zip(args.bound,
                               args.unbound,
                               args.bound_chains,
                               args.unbound_chains))):
        # Load structures
        b_struc, ub_struc = loadBUbStr(b, ub)
        if b_struc == None: continue

        b_chains, ub_chains, b_len, ub_len = getChains(b_struc, ub_struc)

        for i_bc, j_bl in zip(b_chains, b_len):
            if i_bc == bc:
                b_chains = [bc]
                b_len = [j_bl]
                break
        for i_ubc, j_ubl in zip(ub_chains, ub_len):
            if i_ubc == ubc:
                ub_chains = [ubc]
                ub_len = [j_ubl]
                break

        if args.v:
            print(b_chains, ub_chains)

        if len(b_chains) == 0 or len(ub_chains) == 0:
            print(b, "or", ub, "contain no protein chain")
            continue

        try:
            loopChains(b_chains, ub_chains, b_struc, ub_struc, b_len, ub_len, b, ub)
        except KeyboardInterrupt:
            sys.exit()

def getChains(b_struc: int, ub_struc: int):
    # get CA of the protein
    b_sel  = atomsel("(protein and name CA)", b_struc)
    ub_sel = atomsel("(protein and name CA)", ub_struc)
    # List all protein chains
    b_chains = sorted(list(set(b_sel.chain)))
    ub_chains = sorted(list(set(ub_sel.chain)))

    bl = []
    for i in b_chains:
        s = atomsel(f"(protein and name CA) and chain {i}", b_struc)
        bl.append(len(s))
    ubl = []
    for i in ub_chains:
        s = atomsel(f"(protein and name CA) and chain {i}", ub_struc)
        ubl.append(len(s))

    return b_chains, ub_chains, bl, ubl

def loadBUbStr(b: str, ub: str):
    try:
        # Load paired bound / unbound structures
        b_struc  = molecule.load("pdb", args.prefix+ b+"."+args.suffix)
        ub_struc = molecule.load("pdb", args.prefix+ub+"."+args.suffix)

        if args.v:
            print(molecule.numframes(b_struc))
            print(molecule.numframes(ub_struc))

        return b_struc, ub_struc
    except KeyboardInterrupt:
        sys.exit()
    except :
        return None, None

def printRMSDinfo(b_struc, ub_struc,
                  b_ch, ub_ch, b_resid,
                  ub_resid, b_sel, ub_sel):
    print(molecule.name(b_struc), molecule.name(ub_struc))
    print(f"set b [ atomselect 0 \"protein and name CA and chain {b_ch} and altloc '' A and resid "+" ".join(map(str, b_resid))+"\"]")
    print(f"set ub [ atomselect 1 \"protein and name CA and chain {ub_ch} and altloc '' A and resid "+" ".join(map(str, ub_resid))+"\"]")
    print(f"set ab [ atomselect 0 \"protein and chain {ub_ch} and altloc '' A \"]")
    print("set f [measure fit $b $ub]")
    print("$ab move $f")
    print("measure rmsd $b $ub")
    print(f"vmd -m pdb_dl/{molecule.name(b_struc)} pdb_dl/{molecule.name(ub_struc)} -e generic_visu.vmd")
    print(len(b_sel), len(ub_sel))
    print([i for i in zip(b_sel.resname, b_sel.altloc, b_sel.resid)])
    print([i for i in zip(ub_sel.resname, ub_sel.altloc, ub_sel.resid)])

def protRmsd(b_struc: int, ub_struc: int, 
             b_resid: ResidList, ub_resid: ResidList,
             b_ch: str, ub_ch: str):
    """
    """
    b_sel  = atomsel(f"protein\
                   and name CA\
                   and chain {b_ch}\
                   and altloc '' A\
                   and resid "+"'"+"' '".join(map(str, b_resid))+"'" , b_struc)
    b_all  = atomsel(f"all" ,  b_struc)
    ub_sel = atomsel(f"protein \
                   and name CA \
                   and chain {ub_ch} \
                   and altloc '' A \
                   and resid "+"'"+"' '".join(map(str, ub_resid))+"'" , ub_struc)

    if args.p:
        printRMSDinfo(b_struc, ub_struc, b_ch, ub_ch, b_resid, ub_resid, b_sel, ub_sel)

        # print(set(b_sel.altloc), set(ub_sel.altloc))
    # Fit and move
    b_fit_ub_matrix = b_sel.fit(ub_sel)
    b_all.move(b_fit_ub_matrix)
    
    # b
    b_sel.write("pdb", f"../sub_str/{molecule.name(b_struc)}_{b_ch}_CA.pdb")
    s = atomsel(f"protein and chain {b_ch}", b_struc)
    s.write("pdb", f"../sub_str/{molecule.name(b_struc)}_{b_ch}_all.pdb")
    # Ub
    ub_sel.write("pdb", f"../sub_str/{molecule.name(ub_struc)}_{ub_ch}_CA.pdb")
    s = atomsel(f"protein and chain {ub_ch}", ub_struc)
    s.write("pdb", f"../sub_str/{molecule.name(ub_struc)}_{ub_ch}_all.pdb")
    
    return b_sel.rmsd(ub_sel), b_sel.resid, ub_sel.resid
    pass

def dnaRmsd(b_struc: int, ub_struc: int, 
             b_resid: ResidList, ub_resid: ResidList,
             b_ch: str, ub_ch: str):
    """
    """
    b_sel  = atomsel(f"protein\
                   and name CA\
                   and chain {b_ch}\
                   and altloc '' A\
                   and resid "+"'"+"' '".join(map(str, b_resid))+"' " +"\
                   and same resid as within 5 of nucleic" , b_struc)
    # Get the positions of the resid in their list,
    # to project on the unbound resid list
    b_residIndex = [b_resid.index(i) for i in b_sel.resid]
    ub_DNAresid = [ub_resid[i] for i in b_residIndex]
    b_all  = atomsel(f"all" ,  b_struc)
    ub_sel = atomsel(f"protein\
                   and name CA\
                   and altloc '' A\
                   and chain {ub_ch}\
                   and resid "+"'"+"' '".join(map(str, ub_DNAresid))+"'" , ub_struc)
    
    if args.d:
        printRMSDinfo(b_struc, ub_struc, b_ch, ub_ch, b_resid, ub_DNAresid, b_sel, ub_sel)

    # Fit and move
    b_fit_ub_matrix = b_sel.fit(ub_sel)
    b_all.move(b_fit_ub_matrix)
    return b_sel.rmsd(ub_sel), b_sel.resid, ub_sel.resid

def loopChains(b_chains: list, ub_chains: list,
               b_struc: int, ub_struc: int, 
               b_len: list, ub_len: list,
               b = None, ub = None):
    """
    """
    for b_chain, bl in zip( b_chains, b_len):
        for ub_chain, ubl in zip(ub_chains, ub_len):
            if ((b_chain == ub_chain) or
                (len(b_chains) == 1 and len(ub_chains) == 1) or
                # (len(b_chains) == 1 and len(ub_chains) > 1) or
                (len(b_chains) > 1 and len(ub_chains) == 1)
                ):
                dispRes = True
                b_resid, ub_resid = alignChains(b_struc, ub_struc, b_chain, ub_chain, v = args.v)
                try:
                    allProtRmsd, _, _ = protRmsd(b_struc, ub_struc, b_resid, ub_resid, b_chain, ub_chain)
                except KeyboardInterrupt:
                    sys.exit()
                except Exception as e:
                    print("###", ",",b, b_chain, ub , ub_chain, ",","###")
                    if args.e:
                        print(e)
                        traceback.print_exception(type(e), e, e.__traceback__)
                    dispRes = False
                try:
                    DnaProtRmsd, dna_bselRes, dna_ubselRes = dnaRmsd(b_struc, ub_struc, b_resid, ub_resid, b_chain, ub_chain)
                except KeyboardInterrupt:
                    sys.exit()
                except Exception as e:
                    print("###", ",",b, b_chain, ub , ub_chain, ",","###")
                    if args.e:
                        print(e)
                        traceback.print_exception(type(e), e, e.__traceback__)
                    dispRes = False
                try:
                    b_Sasas = sasa(b_struc, b_chain, dna_bselRes)
                    ub_Sasas = sasa(ub_struc, ub_chain, dna_ubselRes)
                except KeyboardInterrupt:
                    sys.exit()
                except Exception as e:
                    print("Sasa error", b, ub)
                    if args.e:
                        traceback.print_exception(type(e), e, e.__traceback__)

                if dispRes:
                    print(b, b_chain, ub , ub_chain, "%.2f" % allProtRmsd, "%.2f" % DnaProtRmsd, \
                       f"{b_chain}/{ub_chain}","%.2f" % allProtRmsd, "%.2f" % DnaProtRmsd, \
                       "{0:.2f}".format(b_Sasas[0]), "{0:.2f}".format(b_Sasas[1]), \
                       "{0:.2f}".format(ub_Sasas[0]), "{0:.2f}".format(ub_Sasas[1]), bl, ubl, sep=",")

def sasa(molId: int, chain: str, restrictSelRes: atomsel = None) -> float:
    """
    """
    s = atomsel(f'protein and chain {chain}', molId)
    r = atomsel(f"protein and chain {chain} and \
    resid "+"'"+"' '".join(map(str,restrictSelRes))+"'", molId)

    # rtcl = evaltcl(f'set a{molId} [ atomselect {molId} "protein and chain {chain}"]'+";"+ \
    #         f'set d{molId} [ atomselect {molId} "protein and chain {chain} and resid '+"'"+"' '".join(map(str,restrictSelRes))+"'\"]"+";"+ \
    #         f"puts \"rtcl [measure sasa 1.4 $a{molId}] [measure sasa 1.4 $a{molId} -restrict $d{molId}]\"")
    # print(rtcl)

    if args.s:
        print(f'set a{molId}{chain} [ atomselect {molId} "protein and chain {chain}"]')
        print(f'set d{molId}{chain} [ atomselect {molId} "protein and chain {chain} and resid '+"'"+"' '".join(map(str,restrictSelRes))+"'\"]")
        print(f"measure sasa 1.4 $a{molId}{chain} -points a{molId}{chain}p")
        print(f"measure sasa 1.4 $a{molId}{chain} -restrict $d{molId}{chain} -points a{molId}{chain}d")
        print(f'draw delete all')
        print(f'draw color 0')
        print(f'foreach pt $a{molId}{chain}p' + " {draw point $pt}")
        print(f'draw color 1')
        print(f'foreach pt $a{molId}{chain}d' + " {draw point $pt}")

    return s.sasa(srad = 1.4), s.sasa(srad = 1.4, restrict = r)

def getAlignedResid(resid: List[int], match: List[str]) -> List[int]:
    """
    Get resid of matching residues
    """
    seqId = resid
    seqPointer = 0
    completedResidList = []

    for m in match:
        if m == "-":
            completedResidList.append(None)
        else:
            completedResidList.append(seqId[seqPointer])
            seqPointer += 1

    if args.v:
        print("completedResidList")
        print(completedResidList)
    return completedResidList

def replaceAA(joined3letterSeq: str):
    """
    Replace unusual AA with their normal counterpart
    joined3letterSeq:void joined 3 letter sequence (ie "".join(...))
    """
    replacePairs = [["mse", "met"]]
    modSeq = joined3letterSeq.lower()
    for o, r in replacePairs:
        modSeq = modSeq.replace(o,r)
    return modSeq

def alignChains(b_struc: int, ub_struc: int,
                b_chain: str, ub_chain: str, v=False) \
                -> (ResidList,ResidList):
    """
    """
    b_chainSel = atomsel("(protein and name CA and chain "+b_chain+" and altloc '' A)", b_struc)
    ub_chainSel = atomsel("(protein and name CA and chain "+ub_chain+" and altloc '' A)", ub_struc)

    b_seq = (seq1(replaceAA("".join(b_chainSel.resname))))
    ub_seq = (seq1(replaceAA("".join(ub_chainSel.resname))))

    # Alignement en penalisant les ouvertures de gap
    alignments = pairwise2.align.globalms(b_seq, ub_seq, 2, -1, -.5, -.1)
    # Get alignement information
    fa = format_alignment(*alignments[0]).split("\n")
    match = fa[1]
    b_match = fa[0]
    ub_match = fa[2]

    # Complete the resid list to take into account gap
    b_residList  = getAlignedResid( b_chainSel.resid, b_match)
    ub_residList = getAlignedResid(ub_chainSel.resid, ub_match)

    b_resid  = [k for j, k in zip(match, b_residList)  if (j == "|" and k != None)]
    ub_resid = [k for j, k in zip(match, ub_residList) if (j == "|" and k != None)]

    if v:
        print("\n".join(fa)[:-1])
        print(len(b_resid), len(ub_resid))
        print([i for i in zip(b_resid, ub_resid)])
    return b_resid, ub_resid


if __name__ == '__main__':
    main()