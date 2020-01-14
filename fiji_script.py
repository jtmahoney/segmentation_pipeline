from ij import IJ, ImagePlus
from ij.process import FloatProcessor, ImageProcessor
from ij.io import DirectoryChooser
from ij.gui import Roi, PolygonRoi
from ij.plugin.frame import RoiManager
import csv
import os
import sys
from array import array
import re
import zipfile

import ast

import platform

gene_list = ['dapi', 'syfp', 'sst', 'gad1', 'pvalb']




dir = DirectoryChooser('Select Directory')


fs = os.listdir(dir.getDirectory())
csv_path = os.path.join(dir.getDirectory(), 'cell_coordinates.csv')
#img_path = os.path.join(dir.getDirectory(), 'input', 'dapi_max-z.png')

print(csv_path)
#print(img_path)

#imp = IJ.openImage(img_path)

#imp.show()


roi_manager = RoiManager()

for gene in gene_list:
	roi_manager.reset()
	with open(csv_path) as csvfile:
		reader = csv.DictReader(csvfile)
		for n, row in enumerate(reader):
	#		print(row['cell_n'])
			poly_name = row['gene_name']
#			poly_name = ast.literal_eval(poly_name)
			if gene == poly_name:
				print(gene, poly_name)
				rr = row['row_pixels']
				cc = row['col_pixels']
				rs = ast.literal_eval(rr)
				cs = ast.literal_eval(cc)
				proi = PolygonRoi(cs, rs, len(rs), Roi.POLYGON)
				roi_manager.addRoi(proi)
		roi_manager.runCommand("Deselect")
		roi_save_path = os.path.join(dir.getDirectory(), gene+"_RoiSet.zip")
		print(roi_save_path)
		
		if not os.path.exists(roi_save_path):
			with zipfile.ZipFile(roi_save_path, "w") as file:
				pass
			file.close()
		
		roi_manager.runCommand("Save", roi_save_path)
	

