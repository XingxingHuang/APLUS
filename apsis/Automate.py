'''
This is the code to Automate Wei's steps for the CLASH pipeline.

CREATED: 05/10/11 (Viana)
UPDATED: 05/19/11 (Viana)
'''

import glob
import os
import pyfits
import shutil
import subprocess



#import prep2
#from prep2 import *


def gather_files(fileset):
	'''
	Finds the correct directories and copies the flt files.
	
	CREATED: 05/19/11 (Viana)
	UPDATED: 05/19/11 (Viana)
	'''
	
	dataset = 'macs1149'
	
	# Gather all the files.
	path = '/data01/cipcal/' + dataset + '/flt/' + dataset + '/'
	# file_list = glob.glob(path + '*_flt.fits')
	file_list = glob.glob(path + '*_fl*.fits')

	# Check/create a directory to work in.
	query = os.access('../Data/',os.F_OK)
	if query == False:
		os.mkdir('../Data/')
	
	# Copy files.
	for filename in file_list:
		root = string.split(filename,'/')[-1]
		query = os.access('../Data/' + root,os.F_OK)
		if query == True:
			print '../Data/' + root + ' already exists.'
		elif query == False:
			print 'Copying: ' + filename + ' --> clashdms1:/Users/viana/Data/'
			shutil.copy(filename,'../Data/')
	
	# Run prep2.py
	run_prep2(dataset,'/Users/viana/Data/')


def acex(filter = 'f555w'):
	'''
	CREATED: 07/??/11 (Viana)
	UPDATED: 08/09/11 (Viana)
	'''
	# Check and create a ref.fits file.
	ingest_path = os.environ['INGEST']
	ref = ingest_path + '/' + filter + '/ref.fits'
	query = os.access(ref,os.F_OK)	
	if query == False:
		os.symlink('/user/viana/CLASH_TEST/ingest/ref.fits',ref)

	# Clean out the output path.
	output_path = os.environ['DATASETS']
	output_path += '/' + filter
	query = os.access(output_path,os.F_OK)
	if query == True:
		shutil.rmtree(output_path)
	
	# Select and execute the process
	if filter in ['f555w']:
		process = '/Users/zheng/acex/bin/acex'
	else:
		print 'Filter not recognized: ' + filter
	subprocess.call([process,filter])


def sextractor(filter = 'f555w'):
	'''
	CREATED: 08/10/11 (Viana)
	UPDATED: 08/11/11 (Viana)
	'''
	# Define the output path.
	output_path = os.environ['DATASETS']
	output_path += '/' + filter + '/'
	
	# Define the input path.
	input_path = os.environ['INGEST']
	input_path += '/' + filter + '/'
	
	# Create files
	query = output_path + 'Images/*_fl*.fits'
	print query
	flt_list = glob.glob(query)
	
	# Check/Create default.param
	query = os.access(os.getcwd(),os.F_OK)
	if query == False:
		shutil.copy2('/user/viana/CLASH/defaul.param',os.getcwd())
	
	f = open('default.sex')
	data = f.readlines()
	f.close()
	

	# Run SExtractor
	ref = input_path + 'ref.fits'
	for flt in flt_list:
		# pdb.set_trace()
		root = string.split(flt,'/')[-1]
		root = string.split(root,'.')[0]
		catalog_name = output_path + root + '.sex'
		catalog_string = '-CATALOG_NAME ' + catalog_name
		print 'sex',flt,catalog_string
		subprocess.call(['sex',flt,catalog_string])


def read_sextractor(filter = 'f555w'):
	'''
	CREATED: 08/11/11 (Viana)
	UPDATED: 08/15/11 (Viana)
	'''
	
	# Define the output path.
	output_path = os.environ['DATASETS']
	output_path += '/' + filter + '/'
	
	# Gather the files.
	query = output_path + '*.sex'
	catalog_list_1 = glob.glob(query)

	# Zero the counters.	
	x_image_max = '0'
	x_image_min = '500000'
	y_image_max = '0'
	y_image_min = '500000'
	
	for catalog in catalog_list:
		f = open(catalog,'r')
		data = f.readlines()
		f.close()
		
		for line in data:
			line = string.split(line)
			if line[0] == '#':
				if line[2] == 'X_IMAGE':
					x_image_column = int(line[1]) - 1
				if line[2] == 'Y_IMAGE':
					y_image_column = int(line[1]) - 1
			else:
				#print line[0],line[x_image_column],line[y_image_column]
				if line[x_image_column] > x_image_max:
					x_image_max = line[x_image_column]
				elif line[x_image_column] < x_image_min:
					x_image_min = line[x_image_column]
				
				if line[y_image_column] > y_image_max:
					y_image_max = line[y_image_column] 
				elif line[y_image_column] < y_image_min:
					y_image_min	= line[y_image_column]
			
		print x_image_max,x_image_min,y_image_max,y_image_min
	
	
def clip_file(filter = 'f555w'):
	'''
	CREATED: 08/11/11 (Viana)
	UPDATED: 08/11/11 (Viana)	
	'''
	### NOTE THIS WILL HAVE TO UPDATED FOR ACEX/WFEX/UVEX RUNS ###
	
	output_name = XXX
	
	array = pyfits.open(file1)
	data = array1[0].data
	array.close()
	data_out = data[x_min:x_max,y_min:y_max]
	
	# Write the fits file.
	hdu = pyfits.PrimaryHDU(data_out)
	hdu.writeto(output_name)

	
def find_shifts('f555w'):
	'''
	CREATED: 08/11/11 (Viana)
	UPDATED: 08/11/11 (Viana)	
	'''
	
	### CREATED inp.lis
	f = open(output_path + 'inp.list','w')
	f.write(input_path + 'ref.fits' + '\n')
	file_list = (output_path +'f*.fits')
	for filename in file_list:	
		f.write(filename + '\n')
	f.close()
	
	subprocess.call(['cl','<','inp.lis'])
	
	
	
if __name__ == '__main__':	
	#acex()
	#sextractor()
	read_sextractor()
