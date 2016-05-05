##############################
#   **Cities Fingerprint**   #
#        MODULE_FIELD        #                                                  
#                            #
#  codename: k1 			 #
#                            #
#   version = 201604291414   #
#    Issues known: none      #
#                            #
#   A script by AndresPDLR   #
##############################

# Calculates a new float field and deletes previous  from join

import processing
import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
from qgis.gui import *
import datetime

edos = ["ags", "bc", "bcs",  "camp", "coah", "col", "chis", "chih",
        "df" , "dgo", "gro",   "gto",  "hgo", "jal",  "mex", "mich",
        "mor", "nay",  "nl",   "oax",  "pue", "qro", "qroo",  "slp",
        "sin", "son", "tab", "tamps", "tlax", "ver",  "yuc",  "zac"]


def phi_field(layer):
    ''''
    >> Creates and populates a PHI field in a given manzanas layer:
        
        Inputs:
            - manzanas layer
        
        Output:
            - updated manzanas layer with PHI field
    '''
    # Create field
    layer.startEditing()
    provider = layer.dataProvider()
    provider.addAttributes([QgsField('PHI', QVariant.Double)])
    layer.updateFields()
    index = layer.fieldNameIndex('PHI')
    # Populate field
    expression = QgsExpression("field_2")
    expression.prepare(layer.pendingFields())
    features = layer.getFeatures()
    for feature in features:
        value = expression.evaluate(feature)
        feature[index] = float(value)
        layer.updateFeature(feature)
    layer.commitChanges()
    print('Done at ', datetime.datetime.now().isoformat())

# Process 1
for edo in edos:
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}_manzanas_urbanas_joined.shp'.format(edo, edo)
    layer = iface.addVectorLayer(path, '{}_manzanas'.format(edo), 'ogr')
    phi_field(layer)
    index_field_1 = layer.fieldNameIndex('field_1')
    index_field_2 = layer.fieldNameIndex('field_2')
    layer.dataProvider().deleteAttributes([index_field_1, index_field_2])
    layer.updateFields()
    print(edo + " completed at ", datetime.datetime.now().isoformat())

def localidad_field(layer):
    ''''
    >> Creates and populates a LOCALIDAD field in a given manzanas layer:
        
        Inputs:
            - manzanas layer
        
        Output:
            - updated manzanas layer with PHI field
    '''
    # Create field
    layer.startEditing()
    provider = layer.dataProvider()
    provider.addAttributes([QgsField('CVEGEO_LOC', QVariant.String)])
    layer.updateFields()
    index = layer.fieldNameIndex('CVEGEO_LOC')
    # Populate field
    expression = QgsExpression("CVEGEO")
    expression.prepare(layer.pendingFields())
    features = layer.getFeatures()
    for feature in features:
        value = expression.evaluate(feature)
        feature[index] = value[0:9]
        layer.updateFeature(feature)
    layer.commitChanges()
    print('Done at ', datetime.datetime.now().isoformat())

# Process 2
for edo in edos:
    path = 'C:/Users/Andres/geo_mex_2015/{}/{}_manzanas_urbanas_joined.shp'.format(edo, edo)
    layer = iface.addVectorLayer(path, '{}_manzanas'.format(edo), 'ogr')
    localidad_field(layer)
    print(edo + " completed at ", datetime.datetime.now().isoformat())