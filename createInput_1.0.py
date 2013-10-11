#!/usr/bin/env python
#
# Create input files to calculate de-biased object discovery probabilities
# If used in research, please cite: 
#   Adams, E. R., A. A. S. Gulbis, J. L. Elliot, S. D. Benecchi, M. W. Buie, D. E. Trilling, and L.H. Wasserman.
#   De-biased Populations of Kuiper Belt Objects from the Deep Ecliptic Survey, submitted to AJ Oct 2013
#
# Created 2013 October by era
#   Create input files to be run by findProb_X.py
#   Note: input can be real objects or a grid for calculation
#
# Note to users: check variables between rows of sharps, e.g.
###################################
# iProbablyWantToBeChanged = toTheDesiredValue
###################################

# System packages (more or less standard depending on your python installation)
import sys
import os
import string
import math
import random
import numpy as np
from time import clock, time
# quad requires Fortran! Does your mac come with Fortran? Go install it if not ###
from scipy.integrate import quad  # integrator
from scipy.optimize import brentq # root finder

###################################
# User packages are stored here on your machine
debiasingModuleDirectory = "/Users/era/Research/Projects/Python/modules/"
###################################
sys.path.append(debiasingModuleDirectory) 
import p4Functions as p4


# Test case values:
testA = 46.027
testE = 0.1081
testI = 11.0025
testH = 7.0
testDiscMag = 20.2
testDiscDist = 42.464
testR = 40.
testObj = "19521"

## extra string for testing
extra=""
## diagnostic test
outputTestObj = False
if outputTestObj:
	print "Running test case:"
	print "a=",testA,"e=",testE,"i=",testI,"H=",testH,"discMag=",testDiscMag,"discR=",testDiscDist,"testR=",testR



#################### Definitions ###############################################

###################################
# Set the highest level directory for the code on your local machine
# All subfolders are automatically created using the code below
mainProjectDir = "/Remote/From_ASTRON/Local/Projects/Library/Bodies/KBOs_and_Centaurs/Deep Ecliptic Survey/Analyses/Post Paper II (2005)/Class Populations/"
scriptDir = mainProjectDir + "Create Probabilities/"
###################################

###################################
# Different probability runs have been assigned different batch numbers (arbitrary string okay below)
ourBatch = "Batch32" 
###################################

mainProbDir = mainProjectDir+"Probabilities/"

realProbDir = mainProbDir + ourBatch + "/"
if os.path.exists(realProbDir) == False:
	print "Creating batch directory"
	os.mkdir(realProbDir)

###################################
# Grids are used to calculate detected fractions
# Note: the _ is a convention so that the grids appear to the top of the folder 
gridDir3to2 = realProbDir+"_Grid_3to2_Big/"
gridDir7to4 = realProbDir+"_Grid_7to4/"
gridDir5to2 = realProbDir+"_Grid_5to2/"
gridDir5to3 = realProbDir+"_Grid_5to3/"
gridDirScattered = realProbDir+"_Grid_Scattered_Big/"
gridDirClassical = realProbDir+"_Grid_Classical_Big/"
gridTest = realProbDir+"_Grid_Test/"

## Classical big grid parameters
#possA = np.arange(37.0, 51.0, 1.0)
#possE = np.arange(0.0, 0.22, 0.02)
#possI = np.arange(2.0, 12.0, 2.0)
#possH = np.arange(4.0, 10.2, 0.2)

## Test grid
possA = np.arange(42.,43.,1)
possE = np.arange(0,0.1,0.1)
possI = np.arange(0.,10,5)
possH = np.arange(6.,7,1)

###################################
# Which directory do you want to run now? It will be made if it doesn't exist
# This can be either a set of real objects...
probabilityDir = realProbDir

# ...or a grid
makeGrid = False
#probabilityDir = gridTest
###################################



print "\nUsing directory: ",probabilityDir,"\n"
if makeGrid == True:
	print "Making a grid with the following parameters:"
	print "a", possA
	print "e", possE
	print "i", possI
	print "h", possH
if os.path.exists(probabilityDir) == False:
	os.mkdir(probabilityDir)
# The input files go in a subdirectory, which gets created if it doesn't exist
inputDir = probabilityDir + "_Input/"
if os.path.exists(inputDir) == False:
	os.mkdir(inputDir)

### How are input files named?
def inputFile(obj,extra):
	return inputDir + extra + obj +".tsv"




## This function is written with default values so you can run a test without real objects
def printInputFile(filename="test.tsv", obj="test", survey="none", objClass="NA", quality="NA", discMag=testDiscMag, hhh=testH, aaa=testA, aErr=0, eee=testE, eErr=0, iii=testI, iErr=0, discDist=testDiscDist, minLong=0, maxLong=360):
	ff = open (inputDir+filename,"w")
	print >>ff, "object\t"+obj
	print >>ff, "survey\t"+survey
	print >>ff, "class\t"+objClass
	print >>ff, "quality\t"+quality
	print >>ff, "discMag\t"+str(discMag)
	print >>ff, "Hmag\t"+str(hhh)
	print >>ff, "a_now\t"+str(aaa)
	print >>ff, "aE_now\t"+str(aErr)
	print >>ff, "e_now\t"+str(eee)
	print >>ff, "eE_now\t"+str(eErr)
	print >>ff, "i_now\t"+str(iii)
	print >>ff, "iE_now\t"+str(iErr)
	print >>ff, "discR\t"+str(discDist)
	print >>ff, "minLong\t"+str(minLong)
	print >>ff, "maxLong\t"+str(maxLong)
	ff.close()


### Functions for estimating median distance and approximate magnitude
def minR(a, e):
	return a * ( 1. - e)
def maxR(a, e):
	return a * ( 1. + e)

## Probability of finding an object at a given distance R
def probObjAtDistR(R,a,e):
	if minR(a, e) < R < maxR(a, e): 
		return ( 3.*R**2 / (2*e* a**3 * (3 + e**2)) )
	else:
		return 0.

## Find the median distance, ie, where 50% of time is spent closer to the sun and 50% further
## renormalized so range is -0.5 to 0.5 for root finding
def percentObjLessThanR(R, aa, ee):
	res, err = quad(probObjAtDistR, minR(aa, ee), R ,args=(aa, ee))
	return res - 0.5

def findMedianDist(aa, ee):
	## for circular orbits
	if (ee <= 10**-6):
		medianDist = aa
	else:
		medianDist = brentq(percentObjLessThanR, minR(aa,ee), maxR(aa,ee), args=(aa,ee))
	return medianDist

aa = testA
ee = testE
#print findMedianDist(aa, ee)


## given distance and H you can estimate the magnitude
def magAtDistRfromHmag(H, R): 
 	return H + 5* math.log(R / 1, 10) + 5 * math.log((R - 1)/ 1, 10)

if makeGrid:
	print "Creating grid inputs..."
	## Run through all a/e/i/H values
	for aaa in possA:
		for eee in possE:
			for iii in possI:
				for hhh in possH:
					filename = "A"+str(aaa)+"E"+str(eee)+"I"+str(iii)+"H"+str(hhh)+".tsv"
					gridObj = "A"+str(aaa)+"E"+str(eee)+"I"+str(iii)+"H"+str(hhh)
					### grid objects need an estimated discovery distance (we use the time-median)
					medianDist = round(findMedianDist(aaa, eee), 2)
					### and an estimated discovery magnitude
					gridDiscMag = round(magAtDistRfromHmag(float(hhh), float(medianDist)), 2)
					### We leave the longitude range as 0-360 even though you might want to limit it to something for resonant classes (not a huge factor)
					printInputFile(filename=filename, obj=gridObj, survey="none", objClass="NA", quality="NA", discMag=gridDiscMag, hhh=hhh, aaa=aaa, aErr=0, eee=eee, eErr=0, iii=iii, iErr=0,  discDist=medianDist, minLong=0, maxLong=360)


if outputTestObj:
	print "Creating test input..."
	printInputFile(filename="test.tsv", obj="test", survey="none", objClass="NA", quality="NA", discMag=testDiscMag, hhh=testH, aaa=testA, aErr=0, eee=testE, eErr=0, iii=testI, iErr=0, discDist=testDiscDist, minLong=0, maxLong=360)
	

### Create real object input files
inputfile = "objectData.tsv"
ff = open(scriptDir+inputfile,"r")
lines = ff.readlines()
for line in lines[1:]:
	elems = line.split("\t")
	obj = elems[0]
	survey = elems[1]
	objClass = elems[2]
	quality = elems[3]
	discMag = elems[4]
	hhh = elems[5]
	aaa = elems[6]
	aErr = elems[7]
	eee = elems[8]
	eErr = elems[9]
	iii = elems[10]
	iErr = elems[11]
	discDist = elems[12]
	minLong = elems[13]
	maxLong = elems[14]
	printInputFile(filename=obj+".tsv", obj=obj, survey=survey, objClass=objClass, quality=quality, discMag=discMag, hhh=hhh, aaa=aaa, aErr=aErr, eee=eee, eErr=eErr, iii=iii, iErr=iErr, discDist=discDist, minLong=minLong, maxLong=maxLong)
print "Created "+str(len(lines)-1)+" input files"

ff.close()
