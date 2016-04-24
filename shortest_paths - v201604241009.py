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

localidades = {"090070001": 1815786, "150330001": 1655015, "140390001": 1495182,
			   "211140001": 1434062, "080370001": 1321004, "020040001": 1300983,
			   "110200001": 1238962, "090050001": 1185772, "141200001": 1142483,
			   "190390001": 1135512, "150580001": 1104585, "080190001": 809232,
			   "150570001": 792211, "310500001": 777615, "090100001": 726664,
			   "240280001": 722772, "010010001": 722250, "260300001": 715061,
			   "050300001": 709671, "020020001": 689775, "250060001": 675773,
			   "190260001": 673616, "120010001": 673479, "151040001": 653410,
			   "230050001": 628306, "220140001": 626495, "090030001": 620416,
			   "150310001": 612383, "050350001": 608836, "160530001": 597511,
			   "280320001": 589466, "140980001": 575942, "090120001": 574577,
			   "071010001": 537102, "090150001": 531831, "100050001": 518709,
			   "150130001": 489160, "151210001": 484573, "190060001": 467157,
			   "280220001": 449815, "190460001": 443273, "090170001": 430978,
			   "301930001": 428323, "300870001": 424755, "090020001": 414711,
			   "141010001": 408759, "090130001": 407885, "090140001": 385439,
			   "090060001": 384326, "250120001": 381583, "280270001": 373725,
			   "090160001": 372889, "151220001": 356352, "270040001": 353577,
			   "190210001": 352444, "110070001": 340387, "170070001": 338650,
			   "180170001": 332863, "280410001": 305155,
			   "090110001": 305076, "260180001": 298625, "280380001": 297284,
			   "150600001": 281799, "020010001": 279765, "150200001": 277959,
			   "190480001": 268347, "161020001": 264439, "100070001": 257352,
			   "130480001": 256584, "200670001": 255029, "240350001": 255015,
			   "211560001": 248716, "150810019": 242272, "090080001": 238431,
			   "300390001": 235983, "040020001": 220389, "050180001": 215271,
			   "030030001": 215178, "260430001": 212533, "151090003": 206081,
			   "140670001": 203342, "070890001": 202672, "280090001": 197216,
			   "120290001": 187251, "301310001": 185242, "150290001": 172919,
			   "040030001": 169466, "150250001": 168720, "170110001": 162427,
			   "090040001": 160491, "110270001": 160169, "260550001": 158089,
			   "070780001": 158027, "151090025": 156191, "170060001": 154358,
			   "190310001": 151893, "050250001": 150178, "230080001": 149923,
			   "300440001": 140896, "220160001": 138878, "060020001": 137383,
			   "050020001": 134233, "060070001": 130035, "320560001": 129011,
			   "300280037": 126507, "240130001": 124644, "320170001": 124623,
			   "190190001": 122627, "150370071": 121470, "320100001": 120944,
			   "301180001": 120844, "280030122": 118614, "120350001": 118468,
			   "080210001": 118071, "080170001": 114007, "260420001": 113836,
			   "260290001": 113082, "301080001": 112046, "150240001": 108449,
			   "150990001": 105165, "080320001": 104836, "150020015": 102667,
			   "201840001": 101810}

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
    
city_simulation("150390001", (10, 10))

# Process
for localidad in localidades.keys():
    city_simulation(localidad, (10, 10))

