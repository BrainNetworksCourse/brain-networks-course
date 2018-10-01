# get publication data for week 2 exzample
# Extracted from notebook

import itertools
import pickle
from Bio import Entrez
from brainnetworks.utils import get_joint_pubs
numpubs={}

researchers={'DB':{'name':'Danielle Bassett'},
            'AF':{'name':'Alex Fornito'},
            'MB':{'name':'Michael Breakspear'},
            'SP':{'name':'Steve Petersen'},
            'MC':{'name':'Michael Cole'},
            'JP':{'name':'Jonathan Power'},
            'DF':{'name':'Damien Fair'},
            'AZ':{'name':'Andrew Zalesky'},
            'LC':{'name':'Luca Cocchi'}}


# first get pubmed search terms from names
for i in researchers:
    n_s=researchers[i]['name'].lower().split(' ')
    researchers[i]['pubmed_name']=n_s[1]+'-'+n_s[0][0]


email='bill@gmail.com'  # email address for use by Entrez
for i in itertools.combinations(list(researchers.keys()),2):
    if not i in numpubs:
        # use cached data if present
        numpubs[i]=get_joint_pubs((researchers[i[0]]['pubmed_name'],
                                   researchers[i[1]]['pubmed_name']),
                                   email)
# save pub data
print(numpubs)
pickle.dump(numpubs,open('../data/pubmed/pubdata.pkl','wb'))
