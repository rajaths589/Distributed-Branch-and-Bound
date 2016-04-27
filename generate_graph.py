import numpy as np
import networkx as nx

m = np.random.rand(8, 8)
m = (m + m.transpose())/2
m = m - np.diag(m.diagonal())
m = np.array(m*10000, dtype=np.int32)
m = np.array(m/100, dtype=np.int32)

#mask = np.random.randint(0,2,size=m.shape).astype(np.bool)
#r = m*mask
#r = ((r + r.transpose())/2).astype(np.int32)

G = nx.from_numpy_matrix(m)
G = G.to_directed()
nx.write_edgelist(G, path="tsp_testcase.txt", delimiter=" ")
