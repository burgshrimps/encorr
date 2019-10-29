import networkx as nx
import matplotlib.pyplot as plt

def create_network_df(all_conn_dfs, pidx):
    G = nx.Graph()

    for conn_df_all_phases in all_conn_dfs:
        conn_df = conn_df_all_phases[pidx]
        conn_m = conn_df.to_numpy()
        for i in range(len(conn_m)):
            for j in range(len(conn_m[0])):
                if conn_m[i,j] != 0:
                    ecolor = 'r' if conn_m[i,j] > 0 else 'b'
                    G.add_edge(conn_df.index[i], conn_df.columns[j], weight=abs(conn_m[i,j]), color=ecolor)

    pos = nx.spring_layout(G)
    #pos = nx.circular_layout(G)

    edges = G.edges()
    ecolors = [G[u][v]['color'] for u,v in edges]
    weights = [G[u][v]['weight']*10 for u,v in edges]

    ncolor_mapper = {'tet1' : 'grey', 
                     'tet2' : 'red',
                     'tet3' : 'orange',
                     'tet4' : 'gold',
                     'tet5' : 'green',
                     'tet6' : 'lime',
                     'tet7' : 'cyan',
                     'tet8' : 'cornflowerblue',
                     'tet9' : 'violet',
                     'tet10' : 'pink',
                     'tet11' : 'peachpuff',
                     'tet12' : 'khaki',
                     'tet13' : 'lavender',
                     'tet14' : 'skyblue',
                     'tet15' : 'palegreen',
                     'tet16' : 'deeppink'}
    ncolors = []
    for node in G.nodes():
        ncolors.append(ncolor_mapper[node.split(':')[0]])

    nx.draw(G, pos, edges=edges, edge_color=ecolors, node_color=ncolors, width=weights, )

    labels = {}
    for node in G.nodes():
        labels[node] = '$' + node.split(':')[1] + '$'

    nx.draw_networkx_labels(G,pos,labels,font_size=12)

    plt.show()

