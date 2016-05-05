##############################
#   **Cities Fingerprint**   #
#        MODULE_JOIN         #                                                  
#                            #
#  codename: k1              #
#                            #
#   version = 201604291218   #
#    Issues known: none      #
#                            #
#   A script by AndresPDLR   #
##############################

'''
Documentation: 
            http://bit.ly/1O0Cgm6
            http://bit.ly/1rnPjt3
'''

import processing
import datetime

edos = ["ags",  "bc", "bcs",  "camp", "coah", "col", "chis", "chih",
        "df" , "dgo", "gro",   "gto",  "hgo", "jal",  "mex", "mich",
        "mor", "nay",  "nl",   "oax",  "pue", "qro", "qroo",  "slp",
        "sin", "son", "tab", "tamps", "tlax", "ver",  "yuc",  "zac"]

for edo in edos:
	path_layer = 'C:/Users/Andres/geo_mex_2015/{}/{}_manzanas_urbanas.shp'.format(edo, edo)
	path_csv = 'C:/Users/Andres/k1_Output_{}.csv'.format(edo)
	path_output = 'C:/Users/Andres/geo_mex_2015/{}/{}_manzanas_urbanas_joined.shp'.format(edo, edo)
	processing.runalg("qgis:joinattributestable",path_layer,path_csv,"CVEGEO","field_1",path_output)
	print(edo + " completed at ", datetime.datetime.now().isoformat())
