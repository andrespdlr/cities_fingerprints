# This script calculates the shortest path for a pair of points 
# Documentation: http://bit.ly/1RAJI94

import processing
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
from qgis.networkanalysis import *

edo = 'df' # Instance

path = 'C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_eje_vial.shp' # Defining the path
vl = iface.addVectorLayer(path, '{}_streets'.format(edo), "ogr") # Uploads the street layer

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
    p = []
    curPos = idStop
    while curPos != idStart:
        p.append(graph.vertex(graph.arc(tree[curPos]).inVertex()).point())
        curPos = graph.arc(tree[curPos]).outVertex();
        
    p.append(tStart)
    rb = QgsRubberBand(qgis.utils.iface.mapCanvas())
    rb.setColor(Qt.red)
    rb.setWidth(2)
    
    for pnt in p:
        rb.addPoint(pnt)

# GRID
cellsize = 50
input = processing.getObject(vl.name())
centerx = (input.extent().xMinimum() + input.extent().xMaximum()) / 2
centery = (input.extent().yMinimum() + input.extent().yMaximum()) / 2
width = (input.extent().xMaximum() - input.extent().xMinimum())
height = (input.extent().yMaximum() - input.extent().yMinimum())
grid='C:/Users/Andres/geo_mex_2015/{}/{}'.format(edo, edo) + '_grid.shp'
processing.runalg("qgis:creategrid", 1, width, height, 0, 0, centerx, centery, crs, grid)
print(height)

processing.alglist()

Random points in extent------------------------------>qgis:randompointsinextent
Random points in layer bounds------------------------>qgis:randompointsinlayerbounds
Random points inside polygons (fixed)---------------->qgis:randompointsinsidepolygonsfixed
Random points inside polygons (variable)------------->qgis:randompointsinsidepolygonsvariable
