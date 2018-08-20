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

journals_to_use=['Cereb. Cortex','Nat Hum Behav',
 'Gigascience', 'Cortex','Nat Neurosci',
 'BMC Neurosci','J. Neurophysiol.',
 'Elife','Neuron',
 'J Cogn Neurosci','J. Neurosci.','Neuroimage','Proc. Natl. Acad. Sci. U.S.A.',
 'Hum Brain Mapp','Biol. Psychiatry','Nat. Med.',
 'Netw Neurosci']
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
def get_authors_from_pmids(records,journal_filter=None):
    authors=[]
    orcid_dict={}
    journals=[]
    for i in records['PubmedArticle']:
        a=[]
        j=i['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        pmid=int(i['MedlineCitation']['PMID'])
        if journal_filter:
            if not j in journal_filter:
                continue
        if 'AuthorList' in i['MedlineCitation']['Article']:
            if not 'eng' in i['MedlineCitation']['Article']['Language']:
                continue
            journals.append(j)
            for au in i['MedlineCitation']['Article']['AuthorList']:
                try:
                    ln=au['LastName'].replace(' ','')
                    initials=au['Initials']
                    if len(initials)>0:
                        a.append(ln+'-'+initials)
                except KeyError:
                    pass
                if 'Identifier' in au:
                    if len(au['Identifier'])>0:
                        try:
                            orcid_dict[ln+'-'+initials]=au['Identifier'][0]
                        except KeyError:
                            pass
        if len(a)>0:
            authors.append(a)
    return(authors,orcid_dict,journals)



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
authors,orcid_dict,journals=get_authors_from_pmids(records,journals_to_use)

copub,all_author_names=get_copub_data(authors)
pickle.dump((copub,all_author_names,orcid_dict),open('../data/pubmed/copub_data.pkl','wb'))
