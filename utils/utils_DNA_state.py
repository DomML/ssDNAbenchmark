import sys, numpy as np
from vmd import molecule, atomsel
from tqdm import tqdm
import itertools as itt

# def tqdm(l):
#     return l

# Load molecule
pref = sys.argv[1]
molNamelist = sys.argv[2:]

def nbInter(molID, allchains):
    nbIntChains = 0
    multipleIntChains = False
    # Check in inter-chain interactions
    for ch in allchains:
        selCh = atomsel("(nucleic \
                 and same chain as within 5 of protein) \
                 and (name N3) \
                 and (not resname DC DT) \
                 and (within 5 of (chain "+ch+" and name N3 and resname DC DT))", molID)
        nbIntChains = len(sorted(set([i for i in selCh.chain] + [ch])))

        if nbIntChains > 1:
            multipleIntChains = True
            break
    return (nbIntChains,multipleIntChains)

def getInteractingChains(molId):
    sel1 = atomsel("(nucleic \
                 and same chain as within 5 of protein) \
                 and (name N3) \
                 and (resname DC DT) \
                 and (within 5 of (name N3 and not resname DC DT))", molId)

    sel2 = atomsel("(nucleic \
                 and same chain as within 5 of protein) \
                 and (name N3) \
                 and (not resname DC DT) \
                 and (within 5 of (name N3 and resname DC DT))", molId)

    allchains = sorted(set(sel1.chain+sel2.chain))
    return allchains, sel1, sel2

def getInteractions(molID, ch):
    outInt = False
    intInt = False

    selCh = atomsel("(nucleic \
             and chain "+ch+" \
             and same chain as within 5 of protein) \
             and (name N3 N1) \
             and not (within 7 of (not chain "+ch+" and name N3  N1))", molID)
    selTmp = atomsel("(nucleic \
             and chain "+ch+" \
             and same chain as within 5 of protein) and (name N3)", molID)

    ####
    ssr = np.array(sorted(set(selCh.resid)))
    allsel = np.array(sorted(set(selTmp.resid)))
    ssDNAp = np.zeros(allsel.shape)
    ssDNAp = [1 if (id in ssr) else 0 for (i,id) in enumerate(allsel)]

    # Unique number of consecutive 1's and 0's
    consSSR  = np.array([sum(1 for _ in group) for _, group in itt.groupby(ssDNAp)])
    # Uniques consecutive 1's and 0's
    valuSSR  = np.array([i for i, group in itt.groupby(ssDNAp)])
    # Unique number of consecutive 1's and 0's longer than 3
    overThr  = np.array([1 if i>=3 else 0 for i in consSSR])
    # Unique number of consecutive 1's longer than 3
    SSRoverThr = valuSSR*overThr
    lenOverThr = consSSR*overThr
    ####

    # Check terminal 1's
    if SSRoverThr[0] == 1 or SSRoverThr[-1] == 1:
        outInt = True
    # Check inside 1's
    if 1 in SSRoverThr[1:-1]:
        intInt = True
    return outInt, intInt, lenOverThr

def getInteractionsRegions(molID, allchains):
    outInt = False
    intInt = False
    allLen = [0]
    for ch in allchains:
        outIntr, intIntr, lenOverThr = getInteractions(molID, ch)
        allLen+=list(lenOverThr)
        if outIntr: outInt = True
        if intIntr: intInt = True

    chsLen = list(set(sorted(allLen)))[1:]

    return outInt, intInt, chsLen
    # , maxLen


print("molName", "len(sel1.chain)", "nbIntChains", "outInt", "intInt", "longestSSDNAchain", sep=",")
for molName in tqdm(molNamelist):
    myMol = molecule.load("pdb", pref+"/"+molName+".pdb")
    # print(molName, file=sys.stderr)

    # Check inter-DNA interactions
    allchains, sel1, sel2 = getInteractingChains(myMol)
    nbIntChains,multipleIntChains = nbInter(myMol, allchains)
    # Locate ssDNA in interacting DNA strand
    outInt, intInt, chsLen = getInteractionsRegions(myMol, allchains)

    print(molName.split("/")[-1], len(sel1.chain), nbIntChains, \
               int(outInt), int(intInt), sep=",")