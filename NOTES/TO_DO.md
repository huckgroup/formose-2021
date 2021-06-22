# To do list

## deal with NorthNet code

Not much of the code is used in these scripts, only these parts are imported:

~~~python
from NorthNet import Classes
from NorthNet import network_generation as n_gen
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops
~~~

The code can either be copied into the repository (any downsides?), or NorthNet must be made ready for 'release' (however, it is undergoing an overhaul, and needs some more work and testing before being released). `Classes` and `network_generation` are the important components.

## Finish creating the reaction expression heatmap (or similar)

Figure 3X needs a representation in which the features of the data can be easily discerned. Getting this to work is a combination of the right metric and the right colour map.
