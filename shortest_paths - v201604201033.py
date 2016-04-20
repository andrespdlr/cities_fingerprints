''' A city's physical structures, once established, may remain in place for more than 150 years '''

##############################
#   **Cities Fingerprint**   #
#        MODULE_PATH		 #                                                  
#                            #
#   		    codename: k1 #
#                            #
#   version = 201604201033   #
#    Issues known: none      #
#                            #
#   A script by AndresPDLR   #
##############################

# This script calculates the shortest path for a pair of points
# Documentation: http://bit.ly/1SfiBCQ

# Requirements -----------------------------------------------------------------
import processing
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
import itertools
import math as mt
import random as rd
import datetime
#from localidades import localidades_dict
# ------------------------------------------------------------------------------

edo_codes = {'01':  'ags', '02':  'bc', '03':  'bcs', '04':  'camp',
			 '05': 'coah', '06': 'col', '07': 'chis', '08':  'chih',
			 '09':   'df', '10': 'dgo', '11':  'gto', '12':   'gro',
			 '13':  'hgo', '14': 'jal', '15':  'mex', '16':  'mich',
			 '17':  'mor', '18': 'nay', '19':   'nl', '20':   'oax',
			 '21':  'pue', '22': 'qro', '23': 'qroo', '24':   'slp',
			 '25':  'sin', '26': 'son', '27':  'tab', '28': 'tamps',
			 '29': 'tlax', '30': 'ver', '31':  'yuc', '32':   'zac'}

edos = ["ags", "bc", "bcs",  "camp", "coah", "col", "chis", "chih",
        "df" , "dgo", "gro",   "gto",  "hgo", "jal",  "mex", "mich",
        "mor", "nay",  "nl",   "oax",  "pue", "qro", "qroo",  "slp",
        "sin", "son", "tab", "tamps", "tlax", "ver",  "yuc",  "zac"]

# Class ------------------------------------------------------------------------
class Crea_layer(object):
    ''' source: http://bit.ly/1VLwb3x '''
    def __init__(self,name,type):
        self.type = type
        self.name = name
        self.layer = QgsVectorLayer(self.type, self.name , "memory")
        self.pr = self.layer.dataProvider()
    def create_poly(self,points):
        self.seg = QgsFeature()  
        self.seg.setGeometry(QgsGeometry.fromPolygon([points]))
        self.pr.addFeatures( [self.seg] )
        self.layer.updateExtents()
    @property
    def disp_layer(self):
        QgsMapLayerRegistry.instance().addMapLayers([self.layer])
# ------------------------------------------------------------------------------

# Functions --------------------------------------------------------------------
def combinations(input_list, subset_size):
    '''
    Inputs:
    	- input_list
    	- subset_size such as 0 <= subset_size <= len(input_list)
    Output:
    	- list with all possible subsets of size subset_size given input_list
    '''
    rv = []
    for subset in itertools.combinations(input_list, subset_size):
    	rv.append(subset)
    return rv
# ------------------------------------------------------------------------------

# Instance =====================================================================
points = [(-99.1016,19.2719), (-99.0679,19.2856), (-99.2174,19.3627), (-99.0982,19.4171), (-99.1634,19.3197), (-99.1524,19.3098)]
points_combinations = combinations(points, 2)
# ==============================================================================

def localidad_layer(localidad, tipo):
    '''
    >> Creates a layer given a localidad
    
    Inputs:
        - localidad geocode (str)
        - tipo: type of layer
        
    Outputs:
        - a new layer
    '''
    edo = edo_codes[localidad[0:2]]
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'
    edo_layer = iface.addVectorLayer(path.format(edo, edo) + '_{}.shp'.format(tipo), '{}_{}'.format(edo, tipo), "ogr")
    filter_localidad = edo_layer.getFeatures(QgsFeatureRequest().setFilterExpression ( ' CVEGEO = ' + localidad ))
    edo_layer.setSelectedFeatures([f.id() for f in filter_localidad])
    QgsVectorFileWriter.writeAsVectorFormat(edo_layer,
                                            path.format(edo, localidad) + '_{}.shp'.format(tipo),
                                            "ascii",
                                            QgsCoordinateReferenceSystem(4326),
                                            "ESRI Shapefile",
                                            True)
    rv = iface.addVectorLayer(path.format(edo, localidad) + '_{}.shp'.format(tipo),
                              '{}_{}'.format(localidad, tipo),
                              "ogr")
    return rv
# ------------------------------------------------------------------------------

# Instances
localidad_layer('090070001', 'eje_vial')
localidad_layer('090070001', 'localidad_urbana_y_rural_amanzanada')

def create_grid(localidad, rows, columns):
    '''
    >> Creates a grid of size rows*columns for a given layer | based on: http://bit.ly/1VLwb3x
    
    Inputs:
        - layer
        - number of rows
        - number of columns
    
    Outputs:
        - a grid layer of size rows*columns
    '''
    layer = localidad_layer(localidad, 'localidad_urbana_y_rural_amanzanada')
    xmin = layer.extent().xMinimum()
    xmax = layer.extent().xMaximum()
    ymin = layer.extent().yMinimum()
    ymax = layer.extent().yMaximum()
    gridHeight = (ymax-ymin)/rows
    gridWidth = (xmax-xmin)/columns
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight
    pol = Crea_layer("{}_grid".format(localidad), "Polygon")
    pr = pol.layer.dataProvider()
    pr.addAttributes([QgsField("G_FID", QVariant.String)])
    pol.layer.updateFields()
    for i in range(int(columns)):
        ringYtop = ringYtopOrigin
        ringYbottom = ringYbottomOrigin
        for j in range(int(rows)):
            poly = [QgsPoint(ringXleftOrigin, ringYtop),
                        QgsPoint(ringXrightOrigin, ringYtop),
                        QgsPoint(ringXrightOrigin, ringYbottom),
                        QgsPoint(ringXleftOrigin, ringYbottom),
                        QgsPoint(ringXleftOrigin, ringYtop)] 
            pol.create_poly(poly)
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth
    iface.setActiveLayer(pol.layer)
    grid_features = pol.layer.getFeatures()
    pol.layer.setSelectedFeatures([f.id() for f in grid_features])
    edo = edo_codes[localidad[0:2]]
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'
    QgsVectorFileWriter.writeAsVectorFormat(pol.layer,
                                            path.format(edo, localidad) + '_grid_{}x{}.shp'.format(rows, columns),
                                            "ascii",
                                            QgsCoordinateReferenceSystem(4326),
                                            "ESRI Shapefile",
                                            True)
    rv = iface.addVectorLayer(path.format(edo, localidad) + '_grid_{}x{}.shp'.format(rows, columns),
                              '{}_grid_{}x{}'.format(localidad, rows, columns),
                              "ogr")
    return rv
# ------------------------------------------------------------------------------

# Instance
create_grid('090070001', 50, 50)
create_grid('090070001', 30, 30)
create_grid("160530001", 50, 50)
create_grid("160530001", 50, 50)
create_grid("230050001", 50, 50)

def random_direction(layer):
    ''''
    >> Populates a DIRECTION field in a given street layer:
        
        Inputs:
            - street layer
        
        Output:
            - Updated street layer wit
    '''
    # Create field
    layer.startEditing()
    provider = layer.dataProvider()
    provider.addAttributes([QgsField('DIRECTION', QVariant.String)])
    layer.updateFields()
    index = layer.fieldNameIndex('DIRECTION')
    # Populate field
    expression = QgsExpression(u'"SENTIDO" = \'Un Sentido\'')
    expression.prepare(layer.pendingFields())
    features = layer.getFeatures()
    for feature in features:
        value = expression.evaluate(feature)
        random_direction = rd.randint(0, 1)
        if random_direction == 0:
            direction = 'U'
        else:
            direction = 'D'
        if bool(value) == True:
            feature[index] = direction
            layer.updateFeature(feature)
        else:
            feature[index] = 'B'
            layer.updateFeature(feature)
    layer.commitChanges()
    print('Done at ', datetime.datetime.now().isoformat())
# ------------------------------------------------------------------------------

# Instance =====================================================================
edo = 'df'
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_eje_vial_sentido.shp'
layer = iface.addVectorLayer(path, 'df_eje_vial_sentido', "ogr")
random_direction(layer)
# ==============================================================================

# Process ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for e in edos:
    ''' For all state level layers, the expected time of running is 30 hrs'''
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(e, e) + '_eje_vial.shp'
    layer = iface.addVectorLayer(path, '{}_eje_vial'.format(e), "ogr")
    random_direction(layer)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def fid_creator(layer):
    '''
    >> Creates a feature id column, it assumes that G_FID already exists
        - G_FID comes from <create_grid(localidad, rows, columns)> method

    Inputs:
        - layer to update
    
    Output:
        - updated layer

    '''
    layer.startEditing()
    index = layer.fieldNameIndex("G_FID")
    localidad = layer.name()[0:9]
    # Populate new field
    features = layer.getFeatures()
    for feature in features:
        fid = feature.id()
        g_fid = str(localidad) + "_" + str(fid)
        feature[index] = g_fid
        layer.updateFeature(feature)
    layer.commitChanges()
    print('Done at ', datetime.datetime.now().isoformat())
# ------------------------------------------------------------------------------

# Instances ========================================
edo = 'df'
localidad = '090070001'
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, localidad) + '_grid_50x50.shp'
layer = iface.addVectorLayer(path, '{}_grid_50x50'.format(localidad), "ogr")
fid_creator(layer)

edo = 'df'
localidad = '090070001'
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, localidad) + '_grid_50x50_clipped.shp'
layer = iface.addVectorLayer(path, '{}_grid_50x50_clipped'.format(localidad), "ogr")
fid_creator(layer)
# ==============================================

def optimal_grid_size(layer, resolution):
    '''
    >> Calculates the optimal number of rows and columns to divide a layer in a grid
    
    Input:
        - layer to be optimized
        - resolution of the grid such as 0 > resolution <= 1 
    
    Output:
        - a dictionary with two values, optimal rows and optimal number of columns given a resolution
    '''
    rv = {}
    xmin = layer.extent().xMinimum()
    xmax = layer.extent().xMaximum()
    ymin = layer.extent().yMinimum()
    ymax = layer.extent().yMaximum()
    max_size = abs(xmax - xmin)*abs(ymax - ymin)
    desired_size = max_size * resolution
    # squared
    if abs(xmax - xmin) >= abs(ymax - ymin):
    	columns = mt.sqrt(desired_size)
    	rows = desired_size / columns
    	rv['C'] = columns
    	rv['R'] = rows
    # rectangular
    else:
    	rows = mt.sqrt(desired_size)
    	columns = desired_size / rows
    	rv['C'] = columns
    	rv['R'] = rows
    return rv
# ------------------------------------------------------------------------------

# Instance
edo = 'df'
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_localidad_urbana_y_rural_amanzanada_090070001.shp'
layer = iface.addVectorLayer(path, 'instance_polygon', "ogr")
optimal_grid_size(layer, .5)


def clipped_grid(localidad, resolution):
    ''' 
    >> Clips a layer given a clip_layer
    
    Inputs:
        - localidad
        - resolution -> tupple
    
    Outputs:
        - a grid layer of size rows*columns
    '''
    edo = edo_codes[localidad[0:2]]
    input_layer = create_grid(localidad, resolution[0], resolution[1])
    clip_layer = localidad_layer(localidad, 'localidad_urbana_y_rural_amanzanada')
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'
    output_path = path.format(edo, localidad) + '_grid_{}x{}_clipped.shp'.format(resolution[0], resolution[1])
    name = '{}_grid_{}x{}_clipped'.format(localidad, resolution[0], resolution[1])
    processing.runalg('qgis:clip', input_layer, clip_layer, output_path)
    rv = iface.addVectorLayer(output_path, name, "ogr")
    return rv
# ------------------------------------------------------------------------------

# Instance ===========================================
clipped_grid('090070001', (50, 50))



def clipped_grid_centroids(localidad, resolution):
    ''' 
    >> Clips a layer given a clip_layer
    
    Inputs:
        - localidad
        - resolution -> tupple
    
    Outputs:
        - a grid layer of size rows*columns
    '''
    edo = edo_codes[localidad[0:2]]
    input_layer = clipped_grid(localidad, resolution)
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'
    output_path = path.format(edo, localidad) + '_grid_{}x{}_clipped_centroids.shp'.format(resolution[0], resolution[1])
    name = '{}_grid_{}x{}_clipped_centroids'.format(localidad, resolution[0], resolution[1])
    processing.runalg('qgis:convertgeometrytype', input_layer, 0, output_path)
    rv = iface.addVectorLayer(output_path, name, "ogr")
    return rv
# ------------------------------------------------------------------------------

# Instance ===========================================
clipped_grid_centroids('090070001', (50, 50))


# Clip (NO FUNCTION)
localidad = '090070001'
edo = edo_codes[localidad[0:2]]
resolution = (50, 50)
input_layer = create_grid(localidad, resolution[0], resolution[1])
clip_layer = localidad_layer(localidad, 'localidad_urbana_y_rural_amanzanada')
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'
output_path = path.format(edo, edo) + '_grid_{}_{}x{}_clipped.shp'.format(localidad,resolution[0], resolution[1])
name = 'grid_{}_{}x{}_clipped'.format(localidad,resolution[0], resolution[1])
clip_layer(input_layer, clip_layer, output_path, name) # NOT WORKING
processing.runalg('qgis:clip', input_layer, clip_layer, output_path)
iface.addVectorLayer(output_path, name, "ogr")

# Extracts Centroids
output_path_centroids = path.format(edo, edo) + '_grid_{}_{}x{}_clipped_centroids.shp'.format(localidad,resolution[0], resolution[1])
processing.runalg('qgis:convertgeometrytype', output_path, 0, output_path_centroids)
iface.addVectorLayer(output_path_centroids, name + '_centroids', "ogr")

# Extracts Centroids Coordinates
def extract_coordinates(edo):
    ''' 
    >> Extracs coordinates of a polygon layer 
    '''
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

# Instance
edo = 'df'
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_localidad_urbana_y_rural_amanzanada_090070001.shp'
layer = iface.addVectorLayer(path, 'instance_polygon', "ogr")

clip_layer(input_layer, layer, output_path)

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


def shortest_path(layer, method, start_point, end_point):
    '''
    >> Finds the shortest path between two points moving along a street network
       Documentation: http://bit.ly/1SfiBCQ
    
    Inputs:
        - street_layer
        - method:
                - 0 two-way & one-way [U, D] 
                - 1 two_way
        - starting point (tupple)
        - ending point (tupple)
    
    Outputs:
    	- Visualize the path in the map
    	- List with nodes of shortest path
    '''
    pStart = QgsPoint(start_point[0], start_point[1])
    pStop = QgsPoint(end_point[0], end_point[1])
    
    # Building the network
    if method == 0:
        director = QgsLineVectorLayerDirector(layer, 10, 'U', 'D', 'B', 3)
        color = Qt.red
    else:
        director = QgsLineVectorLayerDirector(layer, -1, '', '', '', 3)
        color = Qt.blue
    properter = QgsDistanceArcProperter()
    director.addProperter(properter)
    crs = qgis.utils.iface.mapCanvas().mapRenderer().destinationCrs()
    builder = QgsGraphBuilder(crs)
    tiedPoints = director.makeGraph(builder, [pStart, pStop])
    graph = builder.graph()
    tStart = tiedPoints[0]
    tStop = tiedPoints[1]
    idStart = graph.findVertex(tStart)
    idStop = graph.findVertex(tStop)
    (tree, cost) = QgsGraphAnalyzer.dijkstra(graph, idStart, 0)
    if tree[idStop] == -1:
        print('Path not found')
    else:
        path_nodes = []
        curPos = idStop
        while curPos != idStart:
            path_nodes.append(graph.vertex(graph.arc(tree[curPos]).inVertex()).point())
            curPos = graph.arc(tree[curPos]).outVertex();
        path_nodes.append(tStart)
        rb = QgsRubberBand(qgis.utils.iface.mapCanvas())
        rb.setColor(color)
        rb.setWidth(2)
        for node in path_nodes:
            rb.addPoint(node)
    print('shortest path has been found at ', datetime.datetime.now().isoformat())
    return path_nodes
# ------------------------------------------------------------------------------

# Instances ====================================================================
edo = 'df'
path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_eje_vial_sentido.shp'
layer = iface.addVectorLayer(path, 'df_eje_vial_sentido', "ogr")
shortest_path(layer, 1, (-99.2174,19.3627), (-99.0982,19.4171))
shortest_path(layer, 0, (-99.2174,19.3627), (-99.0982,19.4171))
shortest_path(layer, 0, (-99.1016,19.2719), (-99.0679,19.2856))
shortest_path(layer, 1, (-99.1016,19.2719), (-99.0679,19.2856))
shortest_path(layer, 1, (-99.1634,19.3197), (-99.1524,19.3098))
shortest_path(layer, 0, (-99.1634,19.3197), (-99.1524,19.3098))
# ==============================================================================

# Instance Process +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for pair in points_combinations:
    ''' points_combinations comes from method <> '''
    shortest_path(layer, 0, pair[0], pair[1])
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

'''
Formula:

delta = average(all_possible_distances) / max{south_north_distance, west_east_distance}

Formula:
delta = average(all_(euclidean_distance/shorthest_path)) ******************************

''' 