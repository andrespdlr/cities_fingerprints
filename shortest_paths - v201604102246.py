''' A city's physical structures, once established, may remain in place for more than 150 years '''

##############################
#   **Cities Fingerprint**   #
#        MODULE_PATH		 #                                                  
#                            #
#   		    codename: k1 #
#                            #
#   version = 201604102246   #
#    Issues known: none      #
#                            #
#   A script by AndresPDLR   #
##############################

# This script calculates the shortest path for a pair of points
# Documentation: http://bit.ly/1RAJI94


# Requirements -----------------------------------------------------------------
import processing
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *
# ------------------------------------------------------------------------------

# Class ------------------------------------------------------------------------
class Crea_layer(object):
    ''' source: http://bit.ly/1VLwb3x '''
    def __init__(self,name,type):
        self.type=type
        self.name = name
        self.layer =  QgsVectorLayer(self.type, self.name , "memory")
        self.pr =self.layer.dataProvider() 
    def create_poly(self,points):
        self.seg = QgsFeature()  
        self.seg.setGeometry(QgsGeometry.fromPolygon([points]))
        self.pr.addFeatures( [self.seg] )
        self.layer.updateExtents()
    @property
    def disp_layer(self):
        QgsMapLayerRegistry.instance().addMapLayers([self.layer])
# ------------------------------------------------------------------------------

edo = 'df' # Instance

path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_eje_vial.shp' # Defining the path
vl = iface.addVectorLayer(path, '{}_streets'.format(edo), "ogr") # Uploads the street layer

# Shortest path ----------------------------------------------------------------
''' It works at localidad urbana level only '''
# Points to be evaluated
pStart = QgsPoint(-99.1016,19.2719)
pStop = QgsPoint(-99.0679,19.2856)

director = QgsLineVectorLayerDirector(vl, -1, '', '', '', 3)
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
    rb.setColor(Qt.red)
    rb.setWidth(2)
    
    for node in path_nodes:
        rb.addPoint(node)

# ================================================================================== #

# ================================== points ================================== #
# GRID
min_x = vl.extent().xMinimum()
max_x = vl.extent().xMaximum()
min_y = vl.extent().yMinimum()
max_y = vl.extent().yMaximum()


centerx = (vl.extent().xMinimum() + vl.extent().xMaximum()) / 2
centery = (vl.extent().yMinimum() + vl.extent().yMaximum()) / 2
width = (vl.extent().xMaximum() - vl.extent().xMinimum())
height = (vl.extent().yMaximum() - vl.extent().yMinimum())

grid='C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_grid.shp'
processing.runalg("qgis:creategrid", 1, width, height, 0, 0, centerx, centery, crs, grid)

random_points = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_random_points.shp'
processing.runalg('qgis:randompointsinextent', min_x, max_x, min_y, max_y, 50, 0.0, random_points)

processing.runalg('qgis:randompointsinextent', 0,1,0,1,1,'C:/Users/Andres/jdvjdvnodv.shp')
# ======================================================================== #

# ================================== random_points ================================== #
import random,sys
# This parameter sets the amount of random points to be generated
pointId = 10
# Test if active layer is vector layer and of type polygon
layer = iface.activeLayer()
sys.stdout.write('RPG: ')
if not (layer and layer.type() == 0 and layer.geometryType() == 2):
    print("No polygon layer selected.")
else:
    # Prepare new temporary editable memory layer
    pointLayer = iface.addVectorLayer("Point?crs="+layer.crs().toWkt(), "random_points", "memory")
    xmin=xmax=ymin=ymax = 0.0
    # Create global bounding box from polygons/features
    for polygon in layer.getFeatures():
        bounds = polygon.geometry().boundingBox()
        xmin = bounds.xMinimum() if bounds.xMinimum() < xmin else xmin
        xmax = bounds.xMaximum() if bounds.xMaximum() > xmax else xmax
        ymin = bounds.yMinimum() if bounds.yMinimum() < ymin else ymin
        ymax = bounds.yMaximum() if bounds.yMaximum() > ymax else ymax                            
    # Iterate until N random points found
    while pointId > 0:
        # Create random point
        xRandom = xmin + (random.random() * (xmax-xmin))
        yRandom = ymin + (random.random() * (ymax-ymin))
        randomPoint = QgsPoint(xRandom,yRandom)
        randomPointGeometry = QgsGeometry.fromPoint(randomPoint)
        # if random_point is inside polygon feature, create new point feature in temporary layer
        for polygon in layer.getFeatures():
            if polygon.geometry().contains(randomPointGeometry):
                pointFeature = QgsFeature()
                pointFeature.setGeometry(randomPointGeometry)
                pointLayer.dataProvider().addFeatures([pointFeature])
                pointId -= 1
                sys.stdout.write('.')
                break
    print(" Ok.")


# ================================== grid_points ================================== #



from math import ceil
canvas= qgis.utils.iface.mapCanvas()
# first layer
layer = canvas.layer(0)
xmin = layer.extent().xMinimum()
xmax = layer.extent().xMaximum()
ymin = layer.extent().yMinimum()
ymax = layer.extent().yMaximum()
gridWidth = .01
gridHeight = .01
rows = ceil((ymax-ymin)/gridHeight)
cols = ceil((xmax-xmin)/gridWidth)
ringXleftOrigin = xmin
ringXrightOrigin = xmin + gridWidth
ringYtopOrigin = ymax
ringYbottomOrigin = ymax-gridHeight
pol = Crea_layer("grid", "Polygon")
for i in range(int(cols)):
    # reset envelope for rows
    ringYtop = ringYtopOrigin
    ringYbottom =ringYbottomOrigin
    for j in range(int(rows)):
        poly = [QgsPoint(ringXleftOrigin, ringYtop),QgsPoint(ringXrightOrigin, ringYtop),QgsPoint(ringXrightOrigin, ringYbottom),QgsPoint(ringXleftOrigin, ringYbottom),QgsPoint(ringXleftOrigin, ringYtop)] 
        pol.create_poly(poly) 
        ringYtop = ringYtop - gridHeight
        ringYbottom = ringYbottom - gridHeight
    ringXleftOrigin = ringXleftOrigin + gridWidth
    ringXrightOrigin = ringXrightOrigin + gridWidth

pol.disp_layer





'''
Formula:

delta = average(all_possible_distances) / max{south_north_distance, west_east_distance}

Formula:
delta = average(all_(euclidean_distance/shorthest_path))

'''