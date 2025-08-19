# abwd_declust: Mainshock and Aftershock Classification via the Woodell & Abrahamson 2014 and 2018 Distance-Window

This repository contains the code to perform the declustering analysis of earthquake catalogues via the distance-window proposed by Woodell & Abrahamson 2014 & 2018 in:

> Wooddell KE, Abrahamson NA. Classification of Main Shocks and Aftershocks in the NGA-West2 Database. Earthquake Spectra. 2014;30(3):1257-1267. doi:10.1193/071913EQS208M

> PEER 2020/02 - Data Resources for NGA-Subduction Project. Yousef Bozorgnia (NGA-Sub PI),  Jonathan P. Stewart (Report Editor)

As the publication and report identified, the method has been applied to the NGA-WEST2 Ground-Motion Database. In essence, the authors proposed three (3) different distance-window metrics from which to differentiate between main-shock and aftershock-type events. The various distance metrics defining the distance window can be visually appreciated via the figure below. In essence, the method utilizes the traditional Time-Distance-Window formulation after Gardner & Knopoff, where the Time Window is kept as per the original method but the distance window is modified. This modification resulted in a more accurate classification of mainshock versus aftershock events since the distance metric now relied on the actual rupture plane and not a preset circular distance from the epicentre.

![Summary of Distance Metrics](https://github.com/AaronJR474/abwd_declust/blob/main/abwd_declust/example/abwd_distance.png)

## 1.0 Structure

The folder [abwd_declust](abwd_declust/) contains the function [abwd_declust_v2_1.py](abwd_declust_v2_1.py) which is used to carry out the declustering analysis.

The folder [example](abwd_declust/example) contains the Jupyter Notebook [NGMDB_ABWD_DECLUST_V2_2.ipynb](abwd_declust/example/NGMDB_ABWD_DECLUST_V2_2.ipynb) which essentially demonstrates the use of the function and subsequent post-processing of the results. Within this folder, is also the data of the New Zealand Ground-Motion Database (NZGMDB) which was used to run the example demonstrated in the Jupyter Notebook.

## 2.0 Rupture Area Creation and Plotting

The function relies on the rupture areas being packaged into shapely geometries as defined in the Jupyter Notebook example [NGMDB_ABWD_DECLUST_V2_2.ipynb](abwd_declust/example/NGMDB_ABWD_DECLUST_V2_2.ipynb). Thus prior to utilizing the function, as done in the example, rupture areas of each earthquake should be plotted for sense-checking. Failure to do so may lead to errors within the function and/or result in the deviance of the method. An example plot of rupture areas for the NZGMDB and subsequent subdivision of the planar geometries are presented in the figures below.

Finally, some rupture areas within earthquake catalogues often have custom geometry as a result of further studies and subsequent simulations. Therefore, as seen in the example, these specific rupture areas must be treated separately and then spliced into the final rupture area shapely geometry list. 

![Plot of All Rupture Araes and Centroids](https://github.com/AaronJR474/abwd_declust/blob/main/abwd_declust/example/rupture_areas_centroids_epicentres.png)

![Plot of Rupture Area Subdivisions](https://github.com/AaronJR474/abwd_declust/blob/main/abwd_declust/example/rupture_areas_sub.png)

## 3.0 Installation and Requirements

The code was created and developed on Windows within the PyCharm IDE using Python 3.12. In terms of installation, simply download the entire repository and follow the steps identified in the example Jupyter Notebook mentioned previously. Additionally, ensure that all packages are installed as identified in [requirements.txt](requirements.txt).

## 4.0 Limitations and Further Improvements

As with any code, limitations are always inherent along with the need for further improvements. Some of the key limitations and thus need for further improvements are as follows:

1. Custom Geometries for larger rupture areas require a "tighter" convex hull to more accurately capture the intended geometry. This can be done manually in GIS software such as ArcGIS and GlobalMapper, however, at the time of development the use of the _convexhull_ from the _Scipy_ module provided the best outcome. For the next update, it is the intention to add a tighter hull around the geometry and provide an updated example.
2. Geometries are packaged in _Shapely's_ format which unfortunately is a bit resource intensive i.e., one (1) computation for a chosen distance metric (e.g., crjb) takes roughly 1 hour to complete for approximately 5000 events. In the next update, alternative more efficient means of treating the geometry are to be implemented in order to reduce computation times.

## 5.0 References

This code was built utilizing resources from the [openquake](https://docs.openquake.org/oq-engine/3.1/openquake.hmtk.seismicity.declusterer.html) declustering package.
