"""
utility functions for analyzing core-nets data
"""

import pandas,numpy
import networkx as nx


def get_tract_graph(FLNe_threshold=10e-4,limit_to_target=False):
  """
    Load data from Markov et al. 2012 study
    limit_to_target: only include target regions
  """
  tracing_data = pandas.read_excel('core-nets/Cercor_2012 Table.xls')

  tract_dict={}
  source_regions=[]
  target_regions=[]
  for i in tracing_data.index:
      stpair=(str(tracing_data.loc[i]['SOURCE']),str(tracing_data.loc[i]['TARGET']))
      source_regions.append(stpair[0])
      target_regions.append(stpair[1])
      if not stpair in tract_dict:
          tract_dict[stpair]=[tracing_data.loc[i]['FLNe']]
      else:
          tract_dict[stpair].append(tracing_data.loc[i]['FLNe'])

  source_regions=list(set(source_regions))
  target_regions=list(set(target_regions))

  if limit_to_target:
      adjacency_mtx = pandas.DataFrame(numpy.zeros((len(target_regions),len(target_regions))),
                                index=target_regions,columns=target_regions)
  else:
      adjacency_mtx = pandas.DataFrame(numpy.zeros((len(source_regions),len(source_regions))),
                                index=source_regions,columns=source_regions)

  for r in adjacency_mtx.index:
      for c in adjacency_mtx.columns:
          if (r,c) in tract_dict:
              adjacency_mtx.loc[r][c] = numpy.mean(tract_dict[(r,c)])

  G = nx.DiGraph()

  for r in adjacency_mtx.index:
      for c in adjacency_mtx.columns:
          if adjacency_mtx.loc[r][c]>FLNe_threshold:
              G.add_edge(r, c, weight=adjacency_mtx.loc[r][c] )

  return(G,adjacency_mtx)
