# Disparity filter for weighted graphs

Implementation of the disparity filter from [Serrano et al.](https://doi.org/10.1073/pnas.0808904106) to extract the backbone of a weighted graph.


## Usage

The original graph must be a `Graph` or `DiGraph` object from [NetworkX](https://networkx.org/) .

The optimal value for `alpha` corresponds to the threshold that minimizes the number of edges while maximizing the number of vertices. Use `plot_optimal_alpha()` to obtain an illustration.

```python
import disparity_filter_weighted_graphs as dfil
import networkx as nx

# Load the graph from a weighted edgelist.
graph = nx.read_weighted_edgelist('path_to_file', create_using=nx.Graph())

# Compute the 'alpha' value for each edge.
dfil.compute_alpha(graph)

# Find the optimal value for alpha. The dataframe used to find the optimal
#   value for alpha is saved to `finding_optimal_alpha.csv.zip`.
dfil.find_optimal_alpha(graph, save_optimal_alpha_data=True)

# Plot the position of the optimal value for alpha.
dfil.plot_optimal_alpha(graph)

# Create a filtered version of the original graph using the optimal value for alpha.
backbone = dfil.filter_graph(graph)
```


## Reference

Extracting the multiscale backbone of complex weighted networks<br>
M. Ángeles Serrano, Marián Boguñá, and Alessandro Vespignani<br>
[Proc. Natl. Acad. Sci. U.S.A. 106, 6483-6488 (2009)](https://doi.org/10.1073/pnas.0808904106)
