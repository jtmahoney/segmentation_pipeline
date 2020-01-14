# this script works with python 3.7

import os
# import cv2
import copy
import tifffile as tif
import numpy as np
import json
import time
from skimage.feature import register_translation

datapath = './cal_images/'
jsonpath = 'shiftset.json'





def readimgs(imagepath):
	return tif.imread(imagepath)


def get_pixel_shift(img1,img2,upscalefactor = 10**3):
    shift, error, phasediff = register_translation(img1,img2,upscalefactor)
    return np.asarray(shift)

def read_json(jsonfilepath):
	ifile = open(jsonfilepath,'r')
	thestring = ifile.read()
	shiftdict = json.JSONDecoder().decode(thestring)
	ifile.close()
	for key,value in shiftdict.items():
		shiftdict[key] = np.asarray(value)
	return shiftdict

def write_json(outpath,thedict):
	for key,value in thedict.items():
		thedict[key] = value.tolist()
	thestring = json.JSONEncoder().encode(thedict)
	# filename = 'theshifts.json' #JT
	filename = outpath
	with open(filename,'w') as f:
		f.write(thestring)
		f.close()
	return


shifts = read_json(jsonpath)


Imagepath = [f for f in os.listdir(datapath) if f.endswith('tif')]
print(Imagepath)
QBimg = [f for f in os.listdir(datapath) if 'dmqb' in f.lower()]
print(os.path.join(datapath,QBimg[0]))
QBimg = readimgs(os.path.join(datapath,QBimg[0]))


for key in shifts.keys():
	print(key.split('_')[-1])

	for item in Imagepath:
		# print item
		if key.split('_')[-1].lower() in item.lower():
			img = readimgs(os.path.join(datapath,item))
			shift = get_pixel_shift(QBimg,img)
			shifts[key] = shift
			print(type(shift))

write_json(jsonpath,shifts)





