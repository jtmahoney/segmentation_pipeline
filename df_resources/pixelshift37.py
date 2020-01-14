# I did not write most of this.  I took code from Rob and modified it so it would run on Python3.7
# This class is called in order to allign multiple channels

from skimage import data, color, io
from skimage import filters
from skimage.feature import canny, register_translation
from skimage.measure import regionprops
import os
import copy
# import tifffile
import numpy as np
from scipy import signal
from scipy import ndimage as nd
import json
import time
import subprocess
import sys
import collections
# from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
# import PyQt5.QtGui as QtGui
# from PyQt5.QtCore import *

cwd =  os.getcwd()
# jsonfilepath = './data_files/shiftset.json'

class ImageConverter():
	def __init__(self,root = None,imagepathID = None,filetype = 'None'):

		self.root = root
		self.imagepathID = imagepathID
		self.filetype = filetype
		# print self.root

		if self.root is not None:
			self.imagepath = os.path.join(self.root,self.imagepathID)
			self.imgpathlist = sorted([f for f in os.listdir(os.path.join(root,self.imagepath)) if f.endswith(str(self.filetype))])
			# print self.imgpathlist, self.imagepath
	@classmethod
	def readimgs(cls,imagepath,filetype = 'tif'):
		return(io.imread(imagepath))


	def convertimgs(self):
		pass

	@classmethod
	def createPixelShifter(cls,imgpathlist):
		imglist = ImageConverter.readimgs(imgpathlist)
		return(PixelShifter(images=imglist, imagepathlist = imgpathlist))


	@staticmethod
	def createimageset(imagefilenames,shift,aligned_to_dichroic = 'DMQB'):
		dataset = {
					'Image': image,
					'Shift': shift,
					'Aligned to': aligned_to_dichroic, 
					'Filenames': filenames}
		return(dataset)

	@staticmethod
	def writeimageset(path,foldername,filename,image,endset = True):
		if endset:
			savepath = os.path.join(path,foldername)
		else:
			savepath = path
		if not os.path.exists(savepath):
			os.makedirs(savepath)

		io.imsave(os.path.join(savepath,filename),image)

		return




class PixelShifter():
	def __init__(self,images = None,imagepathlist = None,jsonfilepath = None):
		self.images = images
		# print type(self.images)
		self.imagepathlist = imagepathlist
		self.jsonfilepath = jsonfilepath
		if jsonfilepath is not None:
			# jsondir = os.path.join(os.getcwd(),jsondir)
			# jsonfilepath = [f for f in os.listdir(jsondir) if f.endswith('json')]
			# print jsonfilepath
			ifile = open(jsonfilepath,'r')
			thestring = ifile.read()
			self.shiftdict = json.JSONDecoder().decode(thestring)
			ifile.close()
			# print self.shiftdict.values()
			for key,value in self.shiftdict.items():
				# print value
				if value is not None:
					# print value[1]
					self.shiftdict[key] = np.asarray(value)
			# print self.shiftdict
		if images is not None:
			self.images = self.shape_images(self.images)


	def get_centroid(self,img):
	    #get centroid of input image, returns xy location as tuple
	    threshold = self.thresholdimg(img)
	    props = regionprops(threshold,img)
	    
	    return(props[0].centroid)

	def readmetadata(self,path):
		metadata = {}
		i = 1
		# print path
		with open(path,'r') as f:
		    for line in f:
		        content = line.split('\t')
		        if 'excitation' in content[0].lower():
		        	# print content[1]
		        	content[1] = content[1].rstrip()
		        	# print content[1]
		        	metadata['Color' + str(i)] = content[1].translate(None,'"')
		        	i += 1
		f.close()
		# print(collections.OrderedDict(sorted(metadata.items()))) #JT
		return(collections.OrderedDict(sorted(metadata.items())))

	def medfilt(self,img,sigma):
	    return(nd.filters.median_filter(img,size = (sigma,sigma)))

	def thresholdimg(self,img):
	    threshold = filters.threshold_otsu(img) #uses otsu method to calculate foreground value
	    thresholdimg = (img>threshold).astype(int) #returns thresholded img as binary array
	    return(thresholdimg)


	def pad_images(self,img1,padsizex = 50,padsizey = 50):

	   	img1 = np.pad(img1,((padsizex,padsizex),(padsizey,padsizey)),'constant',constant_values = (0,0))
	   	return(img1)

	def shape_images(self,images = None):
		if images is None:
			return(None)
		else:
			if len(images.shape) > 2:
				# print len(images)
				newimgs = []
				for i in range(len(images)):
					# print 'image shape:',images[i].shape, type(images[i,:,:])
					img = self.pad_images(images[i,:,:])
					newimgs.append(img)
				return(np.asarray(newimgs))
			else:
				return(self.pad_images(images))


	def align_images(self,img1,img2, weighted = False, threshold = []):    
	    shift = self.get_pixel_shift(img1,img2)
	    # print 'shift size:',shift
	    


	    img2freq = np.fft.fftn(img2)
	    reg = nd.fourier.fourier_shift(img2freq,shift)
	    reg = np.abs(np.fft.ifftn(reg))
	    reg = reg.astype(np.uint16)

	    return(reg, shift)



	def get_pixel_shift(self,img1,img2,upscalefactor = 10**3):
	    shift, error, phasediff = register_translation(img1,img2,upscalefactor)
	    # print error
	    return(shift)

	def align_with_shift_matrix(self,image,shiftval = None):
		
		if shiftval is None:
			return(image)
		# print shiftval
		# print image.shape

		imgfreq = np.fft.fftn(image)
		reg = nd.fourier.fourier_shift(imgfreq,shiftval)
		reg = np.abs(np.fft.ifftn(reg))
		image = reg.astype(np.uint16)
		return(image)

	def align_tile_set(self,tileset,keys):

		for i in range(len(tileset)):
			# print keys[i]
			tileset[i] = self.align_with_shift_matrix(tileset[i],keys[i])
			
		return



# class Executor(QWidget):
# 	def __init__(self,datapath = None,jsonpath = None):
# 		super(Executor,self).__init__()
# 		self.jsonpath = jsonpath
# 		self.shifter = PixelShifter(jsonfilepath = self.jsonpath)
# 		self.datapath = datapath
# 		self.outputfolder = 'Aligned'
# 		if self.datapath == None:
# 			self.datapath = self.getdatapath()
		
		
# 		self.runregistration()
# 		self.quit()
		


# 	def getdatapath(self):
# 		datapath = QFileDialog.getExistingDirectory(caption='choose directory',options=QFileDialog.ShowDirsOnly)
# 		datapath = str(datapath)
# 		print(os.listdir(datapath))
# 		return(datapath)


# 	def quit(self):
# 		print('quitting')
# 		QCoreApplication.quit()

# 	def runregistration(self):
# 		if self.datapath == None:
# 			self.quit()
# 		if self.jsonpath == None:
# 			self.quit()
# 		Imagepathlist = []
# 		path = os.path.join(cwd,self.datapath)
# 		for item in os.listdir(path):
# 			if item.endswith('tif'):
# 				Imagepathlist.append(os.path.join(path,item))
# 			elif item.endswith('txt'):
# 				metadatapath = os.path.join(path,item)
# 				metadata = self.shifter.readmetadata(metadatapath)
# 		Imagepathlist = sorted(Imagepathlist)
# 		# print 'image paths',Imagepathlist
# 		# for item in Imagepathlist:
# 		# 	print item.split('/')[-1]





# 		print(self.shifter.shiftdict.keys())
# 		time0 = time.time()
# 		for i in range(len(Imagepathlist)):

# 			# print i
# 			# print Imagepathlist[i]

# 			filename = Imagepathlist[i].split(os.sep)[-1]
# 			for key in self.shifter.shiftdict.keys():
# 				# print key.split('_')[-1].lower()
# 				# print filename.lower()
# 				if key.split('_')[-1].lower() in filename.lower():
# 					Image = ImageConverter.readimgs(Imagepathlist[i])
# 					# print Image.shape #JT
# 					Image = self.shifter.shape_images(Image)
# 					shiftval = self.shifter.shiftdict[key]
# 					# print shiftval
# 					break
# 			Image = self.shifter.align_with_shift_matrix(Image,shiftval)
# 			ImageConverter.writeimageset(self.datapath,self.outputfolder,filename,Image)
# 		print('total time:', time.time()-time0)
# 		return



# if __name__ == '__main__':
# 	app = QApplication(sys.argv)
# 	print('hello')
# 	ex = Executor(jsonpath = jsonfilepath)
# 	sys.exit()
