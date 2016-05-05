''' A city's physical structures, once established, may remain in place for more than 150 years '''

##############################
#   **Cities Fingerprint**   #
#        MODULE_PATH		 #                                                  
#                            #
#  codename: k1 			 #
#                            #
#   version = 201604241009   #
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
import math as mt
import csv

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

def random_direction(layer):
    ''''
    >> Populates a DIRECTION field in a given street layer:
        
        Inputs:
            - street layer
        
        Output:
            - updated street layer with DIRECTION field
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
    fid_creator(input_layer)
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'
    output_path = path.format(edo, localidad) + '_grid_{}x{}_clipped_centroids.shp'.format(resolution[0], resolution[1])
    name = '{}_grid_{}x{}_clipped_centroids'.format(localidad, resolution[0], resolution[1])
    processing.runalg('qgis:convertgeometrytype', input_layer, 0, output_path)
    rv = iface.addVectorLayer(output_path, name, "ogr")
    return rv

def extract_centroid_coordinates(layer):
    '''
    >> Extracs coordinates of a centroid layer
    
    Inputs:
        - centroid layer
    Outputs:
        - dictionary of centroid coordinates with g_fid as keys
    '''
    centroids_coordinates = {}
    features = layer.getFeatures()
    for feature in features:
        g_fid = str(feature['G_FID'])
        geometry = feature.geometry()
        x = geometry.asPoint().x()
        y = geometry.asPoint().y()
        centroids_coordinates[g_fid] = (x, y)
    return centroids_coordinates

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

def reverse_dictionary(dictionary):
    ''' 
    >> transform values into keys and kayes into values
    
    Input:
        - dictionary
    
    Output:
        - reversed dictionary
    '''
    rv = {}
    for k in dictionary.keys():
        rv[dictionary[k]] = k
    return rv

def combinations_calculator(n, r):
    n_f = mt.factorial(n)
    r_f = mt.factorial(r)
    nr_f = mt.factorial(n - r)
    rv = n_f / (r_f * nr_f)
    return rv

def euclidean_distance(p1, p2):
    ''' Calculates euclidiean distance between two points '''
    x_1 = p1[0]
    y_1 = p1[1]
    x_2 = p2[0]
    y_2 = p2[1]
    distance = mt.sqrt(pow((x_2-x_1), 2) + pow((y_2-y_1), 2))
    return distance

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
        pass
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
        distances = []
        for i, n in enumerate(path_nodes):
            if i != len(path_nodes) - 1:
                d = euclidean_distance(path_nodes[i], path_nodes[i + 1])
                distances.append(d)
        dijkstra_d = sum(distances)
        euclidean_d = euclidean_distance(path_nodes[0], path_nodes[-1])
        delta = euclidean_d / dijkstra_d
        return (euclidean_d, dijkstra_d)
# ------------------------------------------------------------------------------

def city_simulation(localidad, resolution):
    '''
    >> 
    Input:
        - localidad
        - resolution
    Output:
        - csv file with e_distance and d_distance
    '''
    print(localidad + " started at ", datetime.datetime.now().isoformat())
    rv = {}
    edo = edo_codes[localidad[0:2]]
    input_layer = clipped_grid_centroids(localidad, resolution)
    input_dictionary = extract_centroid_coordinates(input_layer)
    reversed_input_dictionary = reverse_dictionary(input_dictionary)
    input_list = input_dictionary.values()
    output_list = combinations(input_list, 2)
    localidad_layer(localidad, 'eje_vial')
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, localidad) + '_eje_vial.shp'
    layer = iface.addVectorLayer(path, '{}_eje_vial'.format(localidad), "ogr")
    for pair in output_list:
        g_fid_1 = reversed_input_dictionary[pair[0]]
        g_fid_2 = reversed_input_dictionary[pair[1]]
        distance_tuple = shortest_path(layer, 0, pair[0], pair[1])
        rv[(g_fid_1, g_fid_2)] = distance_tuple
    writer = csv.writer(open('C:/Users/Andres/{}_simulation_results.csv'.format(localidad), 'wb'))
    for key, value in rv.items():
        writer.writerow([str(key), value])
    print(localidad + " completed at ", datetime.datetime.now().isoformat())