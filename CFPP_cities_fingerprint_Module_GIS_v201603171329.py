                                                  ##############################
                                                  #   **Cities Fingerprint**   #
                                                  #         MODULE_GIS         #                                                  
                                                  #                            #
                                                  #   		    codename: k1   #
                                                  #                            #
                                                  #   version = 201603171329   #
                                                  #    Issues known: none      #
                                                  #                            #
                                                  #   A script by AndresPDLR   #
                                                  ##############################


'''

This program extracs the fingerprint of a urban landscape according to the method
exposed in A typology of Street Patterns by Rémi Louf, Marc Barthelemy (ATSP) 

    http://bit.ly/1QZ83IK


Python Console - QGis Documentation: 

	http://bit.ly/1QJeMoz


- It assumes the existence of polygon layers  
- It operates in three steps:
	- <Step 0> Select the urban blocks from mexican geography and extract nodes and centroids
	- <Step 1> From a polygon layer extracts nodes and centroids of the polygon
	- <Step 2> From the original polygon layer extracts areas and from point layers
			   (centroids and nodes) extracts coordinates.


** NOTE: This code can only be used in QGIS Python Console **

'''


# Requirements -----------------------------------------------------------------
import datetime

import math as mt
from sys import argv
import csv

# GIS Analysis
import processing
import os
from PyQt4.QtCore import QVariant
from PyQt4.QtCore import *
#-------------------------------------------------------------------------------

''' The database with the layers to be used is at a mexican state (edo) level '''

edos = ["ags",  "bc", "bcs",  "camp", "coah", "col", "chis", "chih",
        "df" , "dgo", "gro",   "gto",  "hgo", "jal",  "mex", "mich",
        "mor", "nay",  "nl",   "oax",  "pue", "qro", "qroo",  "slp",
        "sin", "son", "tab", "tamps", "tlax", "ver",  "yuc",  "zac"] # List of 32 mexican states acronyms

# Functions --------------------------------------------------------------------

def extract_areas(edo):
    ''' Extracs the area of a polygon '''
    # Input: name of a state (edo)
    # Output: csv with all areas for each polygon in the layer
    areas = {}
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + "_manzanas_urbanas.shp"
    manzanas_urbanas_layer = iface.addVectorLayer(path, edo + "_manzanas_urbanas", "ogr")
    manzanas_urbanas_features = manzanas_urbanas_layer.getFeatures()
    for feature in manzanas_urbanas_features:
        manzana_ID = feature['CVEGEO']
        geometry = feature.geometry()
        m_area = geometry.area()
        if manzana_ID in areas.keys():
            areas[manzana_ID].append(m_area) 
        else:
            areas[manzana_ID] = m_area
    writer = csv.writer(open('C:/Users/Andres/k1_areas_{}.csv'.format(edo), 'wb'))
    for key, value in areas.items():
        writer.writerow([str(key), value])


def extract_centroid_coordinates(edo):
    ''' Extracs centroids coordinates of a polygon layer '''
    # Input: name of a state (edo)
    # Output: csv with centroids coordinates
    centroids_coordinates = {}
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + "_manzanas_urbanas_centroids.shp"
    centroids_layer = iface.addVectorLayer(path, "{}_centroids_layer".format(edo), "ogr")
    centroid_features = centroids_layer.getFeatures()
    for feature in centroid_features:
        manzana_ID = feature['CVEGEO']
        geometry = feature.geometry()
        x = geometry.asPoint().x()
        y = geometry.asPoint().y()
        if manzana_ID in centroids_coordinates.keys():
            centroids_coordinates[manzana_ID].append((x, y))
        else:
            centroids_coordinates[manzana_ID] = (x, y)
    writer = csv.writer(open('C:/Users/Andres/k1_centroid_coordinates_{}.csv'.format(edo), 'wb'))
    for key, value in centroids_coordinates.items():
        writer.writerow([str(key), value])


def extract_nodes_coordinates(edo):
    ''' Extracs nodes coordinates from a polygon layer '''
    # Input: name of a state (edo)
    # Output: csv with nodes coordinates
    nodes_coordinates = {}
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + "_manzanas_urbanas_nodes.shp"
    nodes_layer = iface.addVectorLayer(path, "{}_nodes_layer".format(edo), "ogr")
    nodes_features = nodes_layer.getFeatures()
    for feature in nodes_features:
        manzana_ID = feature['CVEGEO']
        geometry = feature.geometry()
        x = geometry.asPoint().x()
        y = geometry.asPoint().y()
        if manzana_ID in nodes_coordinates.keys():
            nodes_coordinates[manzana_ID].append((x, y))
        else:
            nodes_coordinates[manzana_ID] = [(x, y)]
        writer = csv.writer(open('C:/Users/Andres/k1_nodes_coordinates_{}.csv'.format(edo), 'wb'))
    for key, value in nodes_coordinates.items():
        writer.writerow([str(key), value])



# Processes --------------------------------------------------------------------
''' <Step 0> '''

for edo in edos:
    layer_path = path.format(edo, edo) + "_manzana.shp"
    basic_layer = iface.addVectorLayer(layer_path, edo + "_manzanas", "ogr") # This calls the layer from the folder
    
    # Filters the manzanas (blocks) layer keeping only the urban ones 
    filter_urban = basic_layer.getFeatures(QgsFeatureRequest().setFilterExpression ( u'"AMBITO" = \'U\'' ) )
    basic_layer.setSelectedFeatures([f.id() for f in filter_urban])
    
    # Writes the selection on a new layer
    QgsVectorFileWriter.writeAsVectorFormat(basic_layer,
    										path.format(edo, edo) + "_manzanas_urbanas.shp",
    										"ascii",
    										QgsCoordinateReferenceSystem(4326),
    										"ESRI Shapefile",
    										True) # Exporting only the selection (True)
    
    # Loading the new layer with manzanas urbanas
    manzanas_urbanas_layer = iface.addVectorLayer(path.format(edo, edo) + "_manzanas_urbanas.shp", edo + "_manzanas_urbanas", "ogr")
    
    # Extracts Centroids & Nodes
    input_layer = path.format(edo, edo) + "_manzanas_urbanas.shp"
    output_layer_centroids = path.format(edo, edo) + "_manzanas_urbanas_centroids.shp" 
    output_layer_nodes = path.format(edo, edo) + "_manzanas_urbanas_nodes.shp"
    # Extraction
    processing.runalg('qgis:convertgeometrytype', input_layer , 0, output_layer_centroids) # 0 for centroids
    processing.runalg('qgis:convertgeometrytype', input_layer , 1, output_layer_nodes) # 1 for nodes
    # Loading layers
    centroids_layer = iface.addVectorLayer(output_layer_centroids, "{}_centroids_layer".format(edo), "ogr")
    nodes_layer = iface.addVectorLayer(output_layer_nodes, "{}_nodes_layer".format(edo), "ogr")


# Given Step 0, through a loop on functions defined above we can perform steps 1 and 2 '''
for e in edos:
	''' <Step 1> '''
    extract_areas(e)
    ''' <Step 2> '''
    extract_centroid_coordinates(e)
    extract_nodes_coordinates(e)



''' Up to this point we should have CSV's with areas and coordinates for blocks, their centroids and nodes respectively '''

