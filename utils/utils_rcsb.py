import pandas as pd
import urllib
import urllib.request
import os
import requests
from vmd import molecule

from rcsbsearch import Attr

# Download newest cluster files
def dl_cluster(co_value, target_dir = "clusters"):
    urllib.request.urlretrieve(f"ftp://resources.rcsb.org/sequence/clusters/bc-{co_value}.out", f"{target_dir}/bc-{co_value}.out")

def dl_pdb(pdb, target_dir = "pdb"):
    try:
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.pdb", f"{target_dir}/{pdb}.pdb")
    except:
        pass

def get_cluster_by_pdb(pdb_chain:str):
    """
    pdb : one "structure.chain" of the cluster for which we want all the structures/chains in
    return : 
        "structure.chain"  list
        cluster id
    """
    rcsb_url = f"https://www.rcsb.org/pdb/rest/sequenceCluster?cluster=100&structureId={pdb_chain}"
    # The query
    response = urllib.request.urlopen(rcsb_url)
    return response.read().decode("utf-8") 
    
def get_cutom_RCSB_report(fields:list, fields_type:list, pdb_id:list):
    """
    fields : list of string, representing the querried fields. Check 'https://www.rcsb.org/pdb/results/reportField.do'
    fields_type : list of type for each field
    pdb_id : list of string, representing the querried PDBs
    00
    return:
        pandas dataframe
    """
    sub_size = 100
    full_pd_q = pd.DataFrame()
    fields_q = ",".join(fields)
    
    for sub_id in range(0,len(pdb_id), sub_size):
        pdb_ids_q = ",".join(pdb_id[sub_id:sub_id+sub_size])
        target = "tmpq"

        rcsb_adress="https://www.rcsb.org/pdb/rest/customReport.xml?"
        rcsb_pdbid=f"pdbids={pdb_ids_q}"
        rcsb_cols=f"&customReportColumns={fields_q}&service=wsfile&format=csv"
        x = rcsb_adress+rcsb_pdbid+rcsb_cols

        urllib.request.urlretrieve(rcsb_adress+rcsb_pdbid+rcsb_cols, "tmpq")

        pd_q = pd.read_csv("tmpq", dtype = {i:j for i,j in zip(fields, fields_type)})
        full_pd_q = pd.concat([full_pd_q, pd_q])
        os.remove("tmpq")
    return full_pd_q


def get_RCSB_query(contains_Protein=True, contains_DNA=True, contains_RNA=True, contains_NA_Hybrid=True, contains_Other=True):
    """
    Contains_* : True or False
    """
    # Setup the query
    prot_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("DNA")
    dna_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("DNA")
    rna_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("RNA")
    hyb_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("NA-hybrid")
    other_q = Attr('entity_poly.rcsb_entity_polymer_type').exact_match("Other")

    #
    if (contains_Protein):
        query = prot_q
    else:
        query = ~prot_q
    if (contains_DNA):
        query = query & dna_q
    else:
        query = query & ~dna_q
    if (contains_RNA):
        query = query & rna_q
    else:
        query = query & ~rna_q
    if (contains_NA_Hybrid):
        query = query & hyb_q
    else:
        query = query & ~hyb_q
    if (contains_Other):
        query = query & other_q
    else:
        query = query & ~other_q
        
    results = set(query())
    return [i for i in results]


def OLD_get_RCSB_query(contains_Protein=True, contains_DNA=True, contains_RNA=True, contains_NA_Hybrid=True):
    """
    Contains_* : True for (Y)es, False for (N)o
    """
    if (contains_Protein):cp = "Y"
    else: cp = "N"
    if (contains_DNA):cd = "Y"
    else: cd = "N"
    if (contains_RNA):cr = "Y"
    else: cr = "N"
    if (contains_NA_Hybrid):ch = "Y"
    else: ch = "N"

    # Setup the query
    query_text = f"""
<orgPdbQuery>
<queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
<containsProtein>{cp}</containsProtein>
<containsDna>{cd}</containsDna>
<containsRna>{cr}</containsRna>
<containsHybrid>{ch}</containsHybrid>
</orgPdbQuery>
"""
    print(query_text)
    url = 'https://www.rcsb.org/pdb/rest/search'
    header = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    # The query
    response = requests.post(url, data=query_text, headers=header)

    if response.status_code == 200:
        return response.text.split('\n')[:-1]
    else:
        print("Failed to retrieve results")
        print(response.status_code)
        return []

def get_numframes(pdb_id):
    m = molecule.load("pdb", f"pdb/{pdb_id}.pdb")
    f = molecule.numframes(m)
    molecule.delete(m)
    return f
    
def get_taxon_dict(taxonID_file="./Taxonomy/taxonId_to_superkingdom.txt"):
    taxon_list = open(taxonID_file).read().split("\n")[:-1]
    taxon_dict = {i.split("\t")[0]:i.split("\t")[1] for i in taxon_list}
    return taxon_dict
    
def taxonID_to_superkingdom(taxonID, taxon_dict):
    """
    """
    try:
        t = taxon_dict[str(taxonID)]
    except Exception as e:
        print(e, taxonID)
        t=None
    return t