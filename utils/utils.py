"""
utils for brain networks course
"""

from Bio import Entrez

# get number of joint pubs for each pair
def get_joint_pubs(author_pair, email,retmax=200000,):
    """
    get all pubs for two authors
    author_pair: tuple containing two pubmed search terms
    email: email address for use by Entrez
    """
    Entrez.email=email # enter your email address here
    assert len(author_pair)==2
    search_term='%s AND %s'%(author_pair[0],author_pair[1])
    handle = Entrez.esearch(db="pubmed", retmax=retmax,term=search_term)
    record = Entrez.read(handle)
    handle.close()
    pmids=list(set([int(i) for i in record['IdList']]))
    return(len(pmids))


