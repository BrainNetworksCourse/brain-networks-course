"""
utility functions for analyzing core-nets data
"""

import pandas,numpy
import networkx as nx


def get_tract_graph(FLNe_threshold=10e-4):
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

  # confirm that we have the correct numbers of regions
  assert len(source_regions)==91
  assert len(target_regions)==29

  # make sure that all of the target regions are present in the source regions list
  assert len(set(source_regions).intersection(set(target_regions))) == len(target_regions)

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
