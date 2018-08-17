"""
utils for brain networks course
"""

from Bio import Entrez
import igraph
import networkx as nx
import tempfile

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


# convert nx graph to igraph, using graphml as intermediary

def nx_to_igraph(G):
    temp_dir = tempfile.mkdtemp()
    tmpfile=os.path.join(temp_dir,'graph.graphml')
    nx.write_graphml(G,tmpfile)
    ig = igraph.read(tmpfile,format="graphml")
    os.remove(tmpfile)
    os.rmdir(temp_dir)
    return(ig)
