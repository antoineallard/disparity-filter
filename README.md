backbone_network
=======================
Takes a NetworkX graph and removes edges to create a "backbone" graph. For more information about this method see the referenced article.



Usage
-----
```
from backbone_network import get_graph_backbone
graph_backbone = get_graph_backbone(graph)
```

Dependencies
------------
* Numpy 1.8.0
* Scipy 0.11.0
* NetworkX 1.8.1

## References

M. A. Serrano et al. (2009) Extracting the Multiscale Backbone of Complex Weighted Networks. PNAS, 106:16, pp. 6483-6488.
