## Version 1.0 (currently under development)

- removed the module architecture to keep it simple
- now uses the analytic expression instead of numerical integration
- corrected a mistake where the second alpha_ij overwrites the first one in undirected graphs
- now adds the alpha values to the original graph as edge attributes instead of creating a new graph
- added functions to compute the optimal alpha and to plot an illustration of the calculation


## Version 0.9

Original notes by Malcolm van Raalte (see https://github.com/malcolmvr/backbone_network @ fd5c692ce758fc95cf128aaaaa65d1b228440665)

This module implements the disparity filter to compute a significance score of edge weights in networks.
Forked from: https://github.com/aekpalakorn/python-backbone-network/blob/master/backbone.py
With the following changes:
 - formatted to pylint standards
 - architected as a module with no code that runs on load
 - broke large functions into smaller ones
 - copy all nodes so that completely disconnected nodes aren't removed and so that node attributes are not removed
 - copy all the original edge attributes so that they are not removed
 - bug fix: changed G.in_degree(G.successors(u)[0]) to G.in_degree(list(G.successors(u))[0])
