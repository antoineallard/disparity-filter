from networkx import DiGraph as _DiGraph
from networkx import Graph as _Graph
from networkx import is_directed as _is_directed
from networkx import isolates as _isolates
from numpy import absolute as _absolute
from numpy import nan as _nan
from numpy import power as _power
from numpy import sqrt as _sqrt
from pandas import DataFrame as _DataFrame
from pandas import concat as _concat



def compute_alpha(G, weight='weight'):
    ''' Compute significance scores (alpha) for weighted edges in G as defined in Serrano et al. 2009
        Args
            G: Weighted NetworkX graph
        References
            M. A. Serrano et al. (2009) Extracting the multiscale backbone of complex weighted networks. PNAS, 106:16, pp. 6483-6488.
    '''
    # Returns
    #     Weighted graph with a significance score (alpha) assigned to each edge
    if _is_directed(G):
        _compute_alpha_directed(G, weight)
        # return _compute_alpha_directed(G, weight)
    else:
        _compute_alpha_undirected(G, weight)
        # return _compute_alpha_undirected(G, weight)



def _compute_alpha_directed(G, weight):
    '''See docstring for `compute_alpha`.'''
    # N = DiGraph()
    # N.add_nodes_from(G.nodes(data=True))

    for u in G:

        k_out = G.out_degree(u)
        k_in = G.in_degree(u)

        if k_out > 1:
            sum_w_out = sum(_absolute(G[u][v][weight]) for v in G.successors(u))
            for v in G.successors(u):
                w = G[u][v][weight]
                p_ij_out = float(_absolute(w))/sum_w_out
                # alpha_ij_out = 1 - (k_out-1) * integrate.quad(lambda x: (1-x)**(k_out-2), 0, p_ij_out)[0] # pylint: disable=cell-var-from-loop
                alpha_ij_out = _power(1 - p_ij_out, k_out - 1)
                G.add_edge(u, v, alpha_out=float('{:.6f}'.format(alpha_ij_out)))
                # N.add_edge(u, v, alpha_out=float('{:.6f}'.format(alpha_ij_out)))
                # N[u][v].update(G[u][v])

        elif k_out == 1 and G.in_degree(list(G.successors(u))[0]) == 1:
            # we need to keep the connection as it is the only way to maintain the connectivity of the network
            # there is no need to do the same for the k_in, since the link is built already from the tail
            v = list(G.successors(u))[0]
            G.add_edge(u, v, alpha_out=float('{:.6f}'.format(0)), alpha_in=float('{:.6f}'.format(0)))
            # N.add_edge(u, v, alpha_out=0., alpha_in=0.)
            # N[u][v].update(G[u][v])

        elif k_out == 1:
            for v in G.successors(u):
                G.add_edge(u, v, alpha_out=1.)
                # N.add_edge(u, v, alpha_out=1.)
                # N[u][v].update(G[u][v])

        if k_in > 1:
            sum_w_in = sum(_absolute(G[v][u][weight]) for v in G.predecessors(u))
            for v in G.predecessors(u):
                w = G[v][u][weight]
                p_ij_in = float(_absolute(w))/sum_w_in
                # alpha_ij_in = 1 - (k_in-1) * integrate.quad(lambda x: (1-x)**(k_in-2), 0, p_ij_in)[0] # pylint: disable=cell-var-from-loop
                alpha_ij_in = _power(1 - p_ij_in, k_in - 1)
                G.add_edge(v, u, alpha_in=float('{:.6f}'.format(alpha_ij_in)))
                # N.add_edge(v, u, alpha_in=float('{:.6f}'.format(alpha_ij_in)))
                # N[v][u].update(G[v][u])

        elif k_in == 1:
            for v in G.predecessors(u):
                G.add_edge(v, u, alpha_in=float('{:.6f}'.format(1)))
                # N.add_edge(v, u, alpha_in=1.)
                # N[v][u].update(G[v][u])

    # return N



def _compute_alpha_undirected(G, weight):
    '''See docstring for `compute_alpha`.'''
    # B = Graph()
    # B.add_nodes_from(G.nodes(data=True))

    for u in G:

        k = len(G[u])

        if k > 1:
            sum_w = sum(_absolute(G[u][v][weight]) for v in G[u])
            for v in G[u]:
                w = G[u][v][weight]
                p_ij = float(_absolute(w))/sum_w
                # alpha_ij = 1 - (k-1) * integrate.quad(lambda x: (1-x)**(k-2), 0, p_ij)[0] # pylint: disable=cell-var-from-loop
                alpha_ij = _power(1 - p_ij, k - 1)
                try:
                    alpha = G[v][u]['alpha']
                    # alpha = B[v][u]['alpha']
                except KeyError:
                    alpha = 1
                G.add_edge(u, v, alpha=float('{:.6f}'.format(min([alpha, alpha_ij]))))
                # B.add_edge(u, v, alpha=float('{:.6f}'.format(min([alpha, alpha_ij]))))
                # B[u][v].update(G[u][v])

        elif k == 1 and G.degree(list(G[u])[0]) == 1:
            # we need to keep the connection as it is the only way to maintain the connectivity of the network
            v = list(G[u])[0]
            G.add_edge(u, v, alpha=float('{:.6f}'.format(0)))
            # B.add_edge(u, v, alpha=0.)
            # B[u][v].update(G[u][v])

        elif k == 1:
            for v in G[u]:
                try:
                    alpha = G[v][u]['alpha']
                except KeyError:
                    alpha = 1
                G.add_edge(u, v, alpha=float('{:.6f}'.format(min([alpha, 1]))))
                # B.add_edge(u, v, alpha=float('{:.6f}'.format(min([alpha, 1]))))
                # B[u][v].update(G[u][v])

    # return B



def filter_graph(G, alpha_t=None, cut_mode='or'):
    ''' Performs a cut of the graph previously filtered through the disparity_filter function.

        Args
        ----
        G: Weighted NetworkX graph

        alpha_t: double (default='0.4')
            The threshold for the alpha parameter that is used to select the surviving edges.
            It has to be a number between 0 and 1.

        cut_mode: string (default='or')
            Possible strings: 'or', 'and'.
            It applies only to directed graphs. It represents the logic operation to filter out edges
            that do not pass the threshold value, combining the alpha_in and alpha_out attributes
            resulting from the disparity_filter function.

        Returns
        -------
        B: Weighted NetworkX graph
            The resulting graph contains only edges that survived from the filtering with the alpha_t threshold

        References
        ---------
        .. M. A. Serrano et al. (2009) Extracting the multiscale backbone of complex weighted networks. PNAS, 106:16, pp. 6483-6488.
    '''
    if _is_directed(G):
        return _filter_graph_directed(G, alpha_t, cut_mode)
    else:
        return _filter_graph_undirected(G, alpha_t)



def _filter_graph_directed(G, alpha_t, cut_mode='or'):
    '''See the docstring for the `filter_graph` function.'''

    if alpha_t == None:
        alpha_t = G.graph['optimal_alpha']

    B = _DiGraph()
    B.add_nodes_from(G.nodes(data=True))

    for u, v, w in G.edges(data=True):

        alpha_in = w['alpha_in']
        alpha_out = w['alpha_out']

        if cut_mode == 'or':
            if alpha_in < alpha_t or alpha_out < alpha_t:
                B.add_edge(u, v)
                B[u][v].update(G[u][v])

        elif cut_mode == 'and':
            if alpha_in < alpha_t and alpha_out < alpha_t:
                B.add_edge(u, v)
                B[u][v].update(G[u][v])

    return B



def _filter_graph_undirected(G, alpha_t):
    '''See the docstring for the `filter_graph` function.'''
    B = _Graph()
    B.add_nodes_from(G.nodes(data=True))

    for u, v, w in G.edges(data=True):

        alpha = w['alpha']

        if alpha < alpha_t:
            B.add_edge(u, v)
            B[u][v].update(G[u][v])

    return B



def _remove_edge_and_count_vertices(x, g):
    g.remove_edge(x.v_source, x.v_target)
    g.remove_nodes_from(list(_isolates(g)))
    return g.number_of_nodes()



def find_optimal_alpha(graph, save_optimal_alpha_data=False):
    graph_copy = graph.copy()

    nb_vertices = graph.number_of_nodes()
    nb_edges = graph.number_of_edges()

    df = _DataFrame([[u, v, d['weight'], d['alpha']] for u, v, d in graph_copy.edges(data=True)],
                      columns=['v_source', 'v_target', 'weight', 'alpha'])

    df.sort_values(by=['alpha', 'weight'], ascending=False, inplace=True)

    df.reset_index(drop=True, inplace=True)

    df['Remaining number of edges'] = graph_copy.number_of_edges() - df.index.values - 1
    df['Remaining number of vertices'] = df.apply(_remove_edge_and_count_vertices, g=graph_copy, axis=1)
    df['Remaining fraction of edges'] = df['Remaining number of edges'] / nb_edges
    df['Remaining fraction of vertices'] = df['Remaining number of vertices'] / nb_vertices

    #df.drop_duplicates(subset=['alpha'], keep='last', inplace=True)

    last_row = {'v_source': _nan, 'v_target': _nan,
                'weight': _nan, 'alpha': _nan,
                'Remaining number of edges': nb_edges, 'Remaining number of vertices': nb_vertices,
                'Remaining fraction of edges': 1, 'Remaining fraction of vertices' :1}
    df = _concat([_DataFrame([last_row]), df])

    df['dist_to_diagonal'] = _absolute(df['Remaining fraction of vertices'] - df['Remaining fraction of edges']) / _sqrt(2)

    df.reset_index(drop=True, inplace=True)

    idx = df['dist_to_diagonal'].idxmax()
    graph.graph['optimal_alpha_data_index'] = idx
    alpha_t = df.iloc[idx]['alpha']
    graph.graph['optimal_alpha'] = alpha_t

    graph.graph['optimal_alpha_data'] = df

    if save_optimal_alpha_data:
        df.to_csv('finding_optimal_alpha.csv.zip', compression='zip')



def plot_optimal_alpha(graph):
    import seaborn as sns
    sns.set_theme(style="ticks", palette="deep")

    df = graph.graph['optimal_alpha_data']
    idx = graph.graph['optimal_alpha_data_index']

    x_elbow = df.iloc[idx]['Remaining fraction of edges']
    y_elbow = df.iloc[idx]['Remaining fraction of vertices']
    x_inter = (x_elbow + y_elbow) / 2
    y_inter = (x_elbow + y_elbow) / 2

    ax = sns.lineplot(x='Remaining fraction of edges', y='Remaining fraction of vertices', data=df)

    ax.plot([0, 1, None, x_inter, x_elbow],
            [0, 1, None, y_inter, y_elbow],
            linestyle='--', linewidth=0.5)
    ax.plot(x_elbow, y_elbow, ls='None', marker='o', markersize=3)

    ax.set_aspect('equal')
    sns.despine(ax=ax)
    ax.get_figure().savefig('finding_optimal_alpha.pdf')
