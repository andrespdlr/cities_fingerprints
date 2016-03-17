# cities_fingerprints

A two-module python API to produce the analysis contained in A typology of Street Patterns by Rémi Louf, Marc Barthelemy (ATSP)

http://bit.ly/1QZ83IK.

Module 1 CFPP_cities_fingerprint_Module_GIS:

Performs the geographic operation and interacts with esri layers and shape files.

- It assumes the existence of polygon layers  
- It operates in three steps:
	- Step 0: Select the urban blocks from mexican geography and extract nodes and centroids
	- Step 1: From a polygon layer extracts nodes and centroids of the polygon
	- Step 2: From the original polygon layer extracts areas and from point layers
			        (centroids and nodes) extracts coordinates.

NOTE: This code can only be used in QGIS Python Console

Module 2 CFPP_cities_fingerprint_Module_DATA:

Performs the data analysis associated with the calculation and distribution of Φ

- It assumes the existence of CSV files produced by MODULE_GIS
- Produces the key variable from ATSP called Φ for each block

This repo contains:
- Module 1 script
- Module 2 script (.py and .ipynb)
- output folder with csv output files
- layers folder, with an instance (due to weight concerns) of three mexican states layer
- images folder, with examples and visualization of the outputs
![Alt text](https://github.com/andrespdlr/cities_fingerprints/blob/master/images/df_map.png "Instance of Map Mexico City")
- presentation with introduction, motivation and some results (cities_fingerprints.pdf)
