import networkx as nx
import pandas as pd

icp = pd.read_csv('./注释文件/ICP_rank.txt',
                  sep="\t")

G_weighted = nx.read_edgelist(
    './共表达网络/net.txt', nodetype=str,
    data=(('cor', float),), create_using=nx.Graph())
# G_weighted=nx.from_pandas_edgelist(df, 'lnc', 'gene_id', create_using=nx.Graph, edge_attr='cor')


pIDs = icp.loc[:, 'ENSG'].tolist()
network_nodes = G_weighted.nodes()

propagate_input = {}
for node in network_nodes:
    propagate_input[node] = 0
    if node in pIDs:
        propagate_input[node] = 1

weighted_personalized_pagerank = nx.pagerank(G_weighted, alpha=0.85, personalization=propagate_input, max_iter=100,
                                             tol=1e-06)

df_metrics = pd.DataFrame(dict(
    pagerank=weighted_personalized_pagerank
))
df_metrics.index.name = 'ENSG'
df_metrics.to_csv('./page-rank结果/result.txt',
                  sep="\t")
