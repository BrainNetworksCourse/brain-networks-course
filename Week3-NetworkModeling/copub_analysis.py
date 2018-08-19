"""
given a set of base authors, find all coauthors and map co-publication

"""

import os,sys
import pandas
from Bio import Entrez
import itertools
import pickle
import time


sys.path.append('../utils')
from utils import get_joint_pubs

Entrez.email='joe@schmo.edu' # enter your email address here
retmax=200000

base_authors=['bassett-d','sporns-o','breakspear-m','milham-m']

# find all coauthors from list of base authors

def get_pmids_from_list(st):
    pmids=[]
    if not isinstance(st, list):
        st=[st]
    for search_term in st:
        handle = Entrez.esearch(db="pubmed", retmax=retmax,term=search_term)
        record = Entrez.read(handle)
        handle.close()
        pmids=pmids+[int(i) for i in record['IdList']]

    return(pmids)


# get records
def get_pubmed_records(pmids):
    handle = Entrez.efetch(db="pubmed", id=",".join(['%d'%i for i in pmids]), retmode="xml")
    records=Entrez.read(handle)
    return(records)

# get full author listing
def get_authors_from_pmids(records):
    authors=[]
    for i in records['PubmedArticle']:
        a=[]
        try:
            if 'AuthorList' in i['MedlineCitation']['Article']:
                for au in i['MedlineCitation']['Article']['AuthorList']:
                    ln=au['LastName'].replace(' ','')
                    initials=au['Initials']
                    a.append(ln+'-'+initials)
            authors.append(a)
        except KeyError:
            pass
    return(authors)



def get_copub_data(authors):
    copub={}
    all_author_names=[]
    for a in authors:
        for i in itertools.combinations(a,2):
            if i[0]==i[1]:
                continue
            for a in i:
                all_author_names.append(a)
            if not i in copub:
                copub[i]=1
            else:
                copub[i]+=1

    all_author_names=list(set(all_author_names))
    return(copub,all_author_names)


pmids=get_pmids_from_list('"resting state" AND "fMRI" AND "connectivity"')
records=get_pubmed_records(pmids)
authors=get_authors_from_pmids(records)
copub,all_author_names=get_copub_data(authors)
pickle.dump((copub,all_author_names),open('../data/pubmed/copub_data.pkl','wb'))
