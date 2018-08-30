"""
utils for brain networks course
"""

from Bio import Entrez
import igraph
import networkx as nx
import tempfile,os,itertools,pickle
import numpy
import sklearn.preprocessing

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


def module_degree_zscore(G,partition,zscore=True):
    """
    given an adjacency matrix and partition,
    return module-degree z score for each node
    set zscore to false to get unnormalized version
    """
    partition=numpy.array(partition)
    unique_levels=numpy.unique(partition)
    adjmtx=nx.to_numpy_array(G)
    degree=numpy.sum(adjmtx,0)
    mdzs=numpy.zeros(len(degree))
    for l in unique_levels:
        a=adjmtx[partition==l,:][:,partition==l]
        d=numpy.sum(a,1)
        if zscore:
            mdzs[partition==l]=sklearn.preprocessing.scale(d)
        else:
            mdzs[partition==l]=d
    return(mdzs)

def participation_coefficient(G,partition):
    mod_degree=module_degree_zscore(G,partition,zscore=False)
    adjmtx=nx.to_numpy_array(G)
    degree=numpy.sum(adjmtx,0)

    pc=1 - (mod_degree/degree)**2
    return(pc)

def get_pubdata(researchers,datafile='../data/pubmed/pubdata.pkl',email='bill@gmail.com'):
    """
    get pubmed data for a set of researchers
    """

    if os.path.exists(datafile):
        numpubs=pickle.load(open(datafile,'rb'))
    else:
        numpubs={}


    # first get pubmed search terms from names
    for i in researchers:
        n_s=researchers[i]['name'].lower().split(' ')
        researchers[i]['pubmed_name']=n_s[1]+'-'+n_s[0][0]


    for i in itertools.combinations(list(researchers.keys()),2):
        if not i in numpubs:
            # use cached data if present
            numpubs[i]=get_joint_pubs((researchers[i[0]]['pubmed_name'],
                                       researchers[i[1]]['pubmed_name']),
                                       email)


    # save pub data
    pickle.dump(numpubs,open(datafile,'wb'))
    return(numpubs)


def mk_random_graph(G_init,verbose=False,maxiter=5):
    edgelist=numpy.random.randint(len(G_init.nodes),
                                  size=(len(G_init.edges),2))
    good_list=False
    iter=0
    while not good_list:
        if iter>maxiter:
            print('hit maxiter')
            return None
        if verbose:
            print(len(edgelist))
        edgelist=edgelist[edgelist[:,0]!=edgelist[:,1]]
        if verbose:
            print('self-edge removal',len(edgelist))
        edgelist=numpy.sort(edgelist,axis=1)
        edgelist=numpy.unique(edgelist,axis=0)
        if verbose:
            print('duplicate removal',len(edgelist))
        if len(edgelist)==len(G_init.edges):
            good_list=True
        else:
            iter+=1
            edgelist=numpy.vstack((edgelist,
                        numpy.random.randint(len(G_init.nodes),
                                size=(len(G_init.edges)-len(edgelist),2)) ))
    G_rand=nx.Graph()
    for i in range(edgelist.shape[0]):
        G_rand.add_edge(edgelist[i,0],edgelist[i,1])

    return(G_rand)

# implementation of Knuth's Algorithm U for finding all set partitions with a given number of blocks
# from https://codereview.stackexchange.com/questions/1526/finding-all-k-subset-partitions

def algorithm_u(ns, m):
    """
    implementation of Knuth's Algorithm U for finding all set partitions with a given number of blocks
    from https://codereview.stackexchange.com/questions/1526/finding-all-k-subset-partitions
    example: algorithm_u([1,2,3,4,5,6,7,8],2)
    """

    def visit(n, a):
        ps = [[] for i in range(m)]
        for j in range(n):
            ps[a[j + 1]].append(ns[j])
        return ps

    def f(mu, nu, sigma, n, a):
        if mu == 2:
            yield visit(n, a)
        else:
            for v in f(mu - 1, nu - 1, (mu + sigma) % 2, n, a):
                yield v
        if nu == mu + 1:
            a[mu] = mu - 1
            yield visit(n, a)
            while a[nu] > 0:
                a[nu] = a[nu] - 1
                yield visit(n, a)
        elif nu > mu + 1:
            if (mu + sigma) % 2 == 1:
                a[nu - 1] = mu - 1
            else:
                a[mu] = mu - 1
            if (a[nu] + sigma) % 2 == 1:
                for v in b(mu, nu - 1, 0, n, a):
                    yield v
            else:
                for v in f(mu, nu - 1, 0, n, a):
                    yield v
            while a[nu] > 0:
                a[nu] = a[nu] - 1
                if (a[nu] + sigma) % 2 == 1:
                    for v in b(mu, nu - 1, 0, n, a):
                        yield v
                else:
                    for v in f(mu, nu - 1, 0, n, a):
                        yield v

    def b(mu, nu, sigma, n, a):
        if nu == mu + 1:
            while a[nu] < mu - 1:
                yield visit(n, a)
                a[nu] = a[nu] + 1
            yield visit(n, a)
            a[mu] = 0
        elif nu > mu + 1:
            if (a[nu] + sigma) % 2 == 1:
                for v in f(mu, nu - 1, 0, n, a):
                    yield v
            else:
                for v in b(mu, nu - 1, 0, n, a):
                    yield v
            while a[nu] < mu - 1:
                a[nu] = a[nu] + 1
                if (a[nu] + sigma) % 2 == 1:
                    for v in f(mu, nu - 1, 0, n, a):
                        yield v
                else:
                    for v in b(mu, nu - 1, 0, n, a):
                        yield v
            if (mu + sigma) % 2 == 1:
                a[nu - 1] = 0
            else:
                a[mu] = 0
        if mu == 2:
            yield visit(n, a)
        else:
            for v in b(mu - 1, nu - 1, (mu + sigma) % 2, n, a):
                yield v

    n = len(ns)
    a = [0] * (n + 1)
    for j in range(1, m + 1):
        a[n - m + j] = j - 1
    return f(m, n, 0, n, a)
