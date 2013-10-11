#!/usr/bin/env python
#
# Calculate de-biased object discovery probabilities
# If used in research, please cite: 
#   Adams, E. R., A. A. S. Gulbis, J. L. Elliot, S. D. Benecchi, M. W. Buie, D. E. Trilling, and L.H. Wasserman.
#   De-biased Populations of Kuiper Belt Objects from the Deep Ecliptic Survey, submitted to AJ Oct 2013
#
# Created 2012 June by era
#   Calculate the detection probability of objects by the DES
#	Takes the output of Mathematica notebooks for real objects
#   v2.0, 2012-10-01
#		First stab at accounting for longitude bias, using libration/Kozai amplitudes for a normalized "allowed" region
#		(It's just a step function.)
#	v2.1, 2012-10-05
#		Output intermediate probabilities (which omit one or more of the bias factors)
#   v3.0, 2013-10-07
#		Split into two scripts: createInput.py handles making the files, findProb.py only calculates probabilities

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


#### Test case values:
useA = 46.027
useE = 0.1081
useI = 11.0025
useH = 7.0
useDiscM = 20.2
useDiscR = 42.464
testR = 40.
testObj = "19521"

extra=""
diag = False
if diag:
	print "Running test case:"
	print "a=",useA,"e=",useE,"i=",useI,"H=",useH,"discMag=",useDiscM,"discR=",useDiscR,"testR=",testR



#################### Definitions ###############################################

###################################
### If this is False, then only objects without existing probability files will be run (saves time!)
### If this is True, then all objects are rerun 
rerunObjects=False
###################################

###################################
# The highest level directory for the code; subfolders are created using the code below
mainProjectDir = "/Remote/From_ASTRON/Local/Projects/Library/Bodies/KBOs_and_Centaurs/Deep Ecliptic Survey/Analyses/Post Paper II (2005)/Class Populations/"
scriptDir = mainProjectDir + "Create Probabilities/"
###################################

###################################
# Different probability runs have been assigned different batch numbers (arbitrary string okay below)
ourBatch = "Batch32" 
###################################

### You likely don't need to modify this (though it should match what's in createInput.py)
mainProbDir = mainProjectDir+"Probabilities/"
realProbDir = mainProbDir + ourBatch + "/"

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


###################################
# Which directory do you want to run now? It will be made if it doesn't exist
# This can be either a set of real objects...
probabilityDir = realProbDir

# ...or a grid
#probabilityDir = gridTest
###################################


print "\nUsing directory: ",probabilityDir,"\n"
if os.path.exists(probabilityDir) == False:
	exit(["Probability directory doesn't exist:  "+probabilityDir])

# The input files are in a subdirectory; the script exits if it doesn't exist
inputDir = probabilityDir + "_Input/"
if os.path.exists(inputDir) == False:
	exit(["Probability input directory doesn't exist:  "+inputDir])


#### Read in all of the extant input files
allPotentialInput = os.listdir(inputDir)

##### Settings from Elliot 2005 aka Paper II --- these values might change at some point...
effMax = 1.0 #(* fixed value in Paper II *)
magRangeFromPaperII = 0.58 #(* value for all VR KBOs from Paper II *)

## To adapt this code to your own survey, you need to know the following for each field:
##  field = identifying string for the chip (DES format: FIELD_DATE_N, where N=1-8 since each field is composed of 8 CCDs)
##  lat and long =central latitude and longitude (ecliptic)
##  halfEff = half-efficiency magnitude
##  tilt = angle frame is tilted with respect to the Kuiper Belt plane
##  longFrac = longitude fraction (span of field in degrees, divided by 360, accounting for frame tilt)
##  usableFrac = amount of field that is searchable, i.e., overlaps two frames, not covered by bright stars/bad pixels/etc.
###########################################
#fieldFile = mainProjectDir+"Summary Files/"+"allFields_Batch31.tsv"  ## full listing of DES fields
fieldFile = scriptDir+"sampleDESfields.tsv"  ## Contains a sample of 5 DES fields x 8 chips for testing
###########################################


## import field info
def importFieldInfo():
	fieldData= open(fieldFile,"r")
	lines = fieldData.readlines()
	fieldNames=[]; fieldLats=[]; fieldLongs=[]; fieldHalfEffMags=[]; fieldTiltAngles=[];  fieldLongFracs=[];  fieldUsableFracs=[]
	for line in lines[1:]:
		elems=line.rstrip("\n").split("\t")
		fieldNames.append(elems[0])
		fieldLats.append(eval(elems[1]))
		fieldLongs.append(eval(elems[2]))
		fieldHalfEffMags.append(eval(elems[3]))
		fieldTiltAngles.append(eval(elems[4]))
		fieldLongFracs.append(eval(elems[5]))
		fieldUsableFracs.append(eval(elems[6]))
	fieldData.close()
	return fieldNames, fieldLats, fieldLongs, fieldHalfEffMags, fieldTiltAngles, fieldLongFracs, fieldUsableFracs

fieldNames, fieldLats, fieldLongs, fieldHalfEffMags, fieldTiltAngles, fieldLongFracs, fieldUsableFracs = importFieldInfo()

############# Equations and helpers ##################
## Location of output probability file
def probFileName(obj,extra):
	return	probabilityDir+ extra+obj+ ".tsv";


## Detection efficiency (Equation 1 from Adams2013)
def detectionEff(fieldN, mag): 
	return (0.5*effMax * (1 + math.tanh( (fieldHalfEffMags[fieldN] - mag) / magRangeFromPaperII)))
if diag:
	print "Detection efficiency of first frame for 22 mag obj: "+str(round(detectionEff(0,22.),4))


## Location probability (Equations 4-5 from Adams2013)
def minR(a, e):
	return a * ( 1. - e)
def maxR(a, e):
	return a * ( 1. + e)

## (Equation 6 from Adams2013)
def probObjAtDistR(R,a,e):
	if minR(a, e) < R < maxR(a, e): 
		return ( 3.*R**2 / (2*e* a**3 * (3 + e**2)) )
	else:
		return 0.

if diag:
	before=time()
	prob=probObjAtDistR(useDiscR,useA, useE)
	after = time()
	print "\nProbability object with a/e is at discR:",prob," which took ", round(after-before,8)," seconds\n"	


def magAtDistRfromHmag(H, R): 
	return (H + 5.*math.log(R/1, 10) + 5. * math.log((R - 1)/1,10))

if diag:
	print "Magnitude at testR given H: "+str(round( magAtDistRfromHmag(useH, testR), 3))


## Maginutde of object at given distance (Equation 3 from Adams2013)
def magObjAtDistR(discMag, discR, R):
	return (discMag + 5. * math.log(R/discR, 10) + 5. * math.log((R-1)/(discR-1), 10))

if diag:
	print "Magnitude at testR: "+str(round( magObjAtDistR(useDiscM, useDiscR, testR), 3))

def integrateMe(R, a, e, discMag, discR, fieldN):
	return probObjAtDistR(R, a,e)*detectionEff(fieldN, magObjAtDistR(discMag, discR, R))

#if diag:
#	print "For instance, at some random values integrate me returns:",integrateMe(42, 43., 0.05, 22., 44., 0)

### Likelihood object is in field given its magnitude, which depends on a/e/H/R/M, also field limiting mag
## (Equation 7 from Adams2013)
def magLikelihoodObjInField(a, e, discMag, discR, useNfields=len(fieldLats)): 
	likelyList=[]
#	print "Finding mag likelihood", a, e, discMag, discR
	## if ecc is very low:
	if e < 10**-6:
		for nn in range(useNfields):
			likelyList.append(detectionEff(nn, magObjAtDistR(discMag, discR, a)))
	else:
		for nn in range(useNfields):
			### This numerical intergration method was plucked from scipy with little thought to alternatives
			res, err = quad(integrateMe, minR(a,e), maxR(a,e),args=(a, e, discMag, discR, nn))
			likelyList.append(res)
	return likelyList


### Likelihood factor based on object inclination, field latitude, tilt of field
## (Equation 9 from Adams2013)
def incLikelihoodObjInField(i, useNfields=range(len(fieldNames))): 
	likelyList=[]
	for nn in useNfields:
		likelyList.append(p4.p4RelativeFieldLikelihood(fieldTiltAngles[nn],fieldLats[nn],i))
	return likelyList	



### NEW: restrict the longitude likelihood based on the allowed longitudes of resonant objects
## (Equation 10 from Adams2013)
def longLikelihoodObjInField(minLong, maxLong, useNfields=len(fieldLats)):
	likelyList=[]
	for nn in range(useNfields):
		thisFieldLong = fieldLongs[nn]
		if (minLong <= thisFieldLong <= maxLong):
			likelyList.append(360/abs(maxLong-minLong))
		else:
			likelyList.append(0)
	return likelyList	



### How are object input files named?
def inputFile(obj,extra):
	return inputDir + extra + obj +".tsv"


##### Functions for reading in object elements
def importElems(obj,extra):
	elemFile=inputFile(obj,extra)
#	print obj, elemFile
	values={}; names=[]
	data = open(elemFile,"r")
	lines = data.readlines()
	for line in lines:
		elems = line.rstrip("\n").split("\t")
		names.append(elems[0])
		values[elems[0]]=elems[1]
	data.close()
	return values,names


## Multiply all items in a list -- really a grabBag function
def multiplyAllListElems(listA):
	return reduce(lambda x,y:x*y, listA)


####### Tests
def runMagTimeTest():
	before=time()
	magLikely=magLikelihoodObjInField(useA, useE, useDiscM, useDiscR,useNfields=100)
	after = time()
	print "\nMagnitude likelihood (a/e/discR/discMag but NOT i),"," which took ", round(after-before,8)," seconds\n"	
	print magLikely

def runIncTimeTest():
	before=time()
	incLikely=incLikelihoodObjInField(useI,useNfields=100)
	after = time()
	print "\nInclination likelihood for field,"," which took ", round(after-before,8)," seconds\n"	
	print incLikely
	
#print "\nImporting test object",testObj
#values,names=importElems(testObj,extra)
#print values
#print names


### Some lazy preprocessing: assume inclination calculation for chip 2 applies to all chips 1-8
### cuts time by 1/8 with little change in results
chipNumbers=[]; chip2Fields=[]; onlyChip2s=[]
proxyChip={}
for ii,ff in enumerate(fieldNames):
	chip = ff[-1]
	field = ff[:-2] ### ditch the chip
	chipNumbers.append(chip)
	if chip == "2":
		onlyChip2s.append(ii)
		chip2Fields.append(field)
for ii in range(len(onlyChip2s)):
	field = chip2Fields[ii]
	proxyChip[field] = ii
#	print field, proxyChip[field]


########## Main Functions #############

### Wrapper function to create a probability file
def createProbabilityFile(obj,extra="",lazy=False):
	probabilityFile = probFileName(obj, extra);
	
	values,names = importElems(obj,extra)
#	print eval(values["a_now"]), eval(values["e_now"]), eval(values["discMag"]), eval(values["discR"])
	magFactor = magLikelihoodObjInField(eval(values["a_now"]), eval(values["e_now"]), eval(values["discMag"]), eval(values["discR"]))
#	print magFactor
	if lazy == True:
		incFactorLazy = incLikelihoodObjInField(eval(values["i_now"]),useNfields=onlyChip2s)
		incFactor=[]
		for ff in fieldNames:
			field = ff[:-2] ### ditch the chip
			incFactor.append( incFactorLazy[ proxyChip[field] ] )
#		print len(incFactorLazy)
#		print len(incFactor)
	else:
		incFactor = incLikelihoodObjInField(eval(values["i_now"]))
	longFactor = longLikelihoodObjInField(eval(values["minLong"]), eval(values["maxLong"]))
		
#	print len(magFactor), len(incFactor),len(fieldUsableFracs),len(fieldLongFracs)
	
	probsDidntFind = []; interProbLong = []; interProbObscured = []; interProbInc = []; interProbMag = []; fracs = 0
	for ii in range(len(magFactor)):
		fracs = fracs +fieldLongFracs[ii]
		probsDidntFind.append(1 - (magFactor[ii] * incFactor[ii] * fieldUsableFracs[ii] * longFactor[ii]) * fieldLongFracs[ii] )
		interProbLong.append(1 - (magFactor[ii] * incFactor[ii] * fieldUsableFracs[ii]) * fieldLongFracs[ii])
		interProbObscured.append(1 - (magFactor[ii] * incFactor[ii] * longFactor[ii])  * fieldLongFracs[ii])
		interProbInc.append(1 - (magFactor[ii] * fieldUsableFracs[ii] * longFactor[ii]) * fieldLongFracs[ii])
		interProbMag.append(1 - (incFactor[ii] * fieldUsableFracs[ii] * longFactor[ii]) * fieldLongFracs[ii])
	totalProbability = 1 - multiplyAllListElems(probsDidntFind)
	# Also output a list of probabilities that omit ONE factor, to see which is the most important in final probability
	totalProbSansLong = 1 - multiplyAllListElems(interProbLong)
	totalProbSansObscured = 1 - multiplyAllListElems(interProbObscured)
	totalProbSansInc = 1 - multiplyAllListElems(interProbInc)
	totalProbSansMag = 1 - multiplyAllListElems(interProbMag)
#	if diag:
#		print "\nTEST", magFactor[0], incFactor[0], fieldLongFracs[0],fieldUsableFracs[0]
	
#	print "Total longitude fraction", fracs
	
	g = open(probabilityFile,"w")
	for nn in names:
		print >>g, nn+"\t"+values[nn]
	print >>g, "totalProbability"+"\t"+str(round(totalProbability,10))
	print >>g, "totalProbSansLong"+"\t"+str(round(totalProbSansLong,10))
	print >>g, "totalProbSansObscured"+"\t"+str(round(totalProbSansObscured,10))
	print >>g, "totalProbSansInc"+"\t"+str(round(totalProbSansInc,10))
	print >>g, "totalProbSansMag"+"\t"+str(round(totalProbSansMag,10))
	g.close()
	return totalProbability

#before=time()
#prob = createProbabilityFile(testObj)
#after = time()
#print "\nProbabilty for ",testObj," is ",prob," and took ", round(after-before,8)," seconds\n"	
######## 9 seconds!!!!! not 1389 or 1600 or 2000!!!!

extra = ""

################## Run 'em all and let me sort them out
# If useTheseObjects is NOT specified, will run ALL input files in inputDir without probability files in probabilityDir
# areWeLazy speeds up code by factor of 1/8 (not applicable to non-DES fields!)
def runAllInputObjects(useTheseObjects=[], areWeLazy=False):
	if useTheseObjects == []:
		allPotentialInput = os.listdir(inputDir+extra)
	else:
		allPotentialInput = useTheseObjects
	
	print "\nFound potential files:",len(allPotentialInput),"at",inputDir+extra
	
	### Which objects don't have probability files yet?
	useObj=[]
	count=0
	for aa in allPotentialInput:
		if ((aa != "_allFields.tsv") & (aa != ".DS_Store")):
			count = count+1
			obj = string.replace(aa,".tsv","")
			if ((rerunObjects == False) & (os.path.isfile(probFileName(obj, extra)))):
				if diag:
					print "File already exists for ",obj
			else:
				useObj.append(obj)
#	useObj = ["fake_5","fake_6","fake_7","fake_8","fake_9"]
	
	print "\nObjects to run:",len(useObj)," (out of", count,")"
	
	before=time()
	for obj in useObj:
		prob = createProbabilityFile(obj,extra=extra,lazy=areWeLazy)
		print "\nProbability for ",obj," is ",round(prob,10)
	after = time()
	
	print "\n\n Running the whole thing took ", round(after-before,8)," seconds\n\n"	


############## Predefined sets
allDESclassicals = ["119067", "126719", "129772", "134860", "135182", "138537", "148780",  "149348", "160091", "160256", "182933", "183963", "184314", "19521",  "1998WG24", "1998WV24", "1998WW31", "1998WX24", "1998WX31",  "1998WY24", "1998WY31", "1999HH12", "1999HJ12", "1999HS11",  "1999HV11", "2000CE105", "2000CF105", "2000CJ105", "2000CL104",  "2000CL105", "2000CN105", "2000CN114", "2000CP104", "2000CQ114",  "2000OH67", "2000ON67", "2000OU69", "2001DD106", "2001FK185",  "2001FK193", "2001FL185", "2001FO185", "2001KA77", "2001KF76",  "2001KH76", "2001KT76", "2001OQ108", "2001QB298", "2001QJ298",  "2001QO297", "2001QP297", "2001QQ297", "2001QQ322", "2001QR297",  "2001QS322", "2001QX297", "2001QZ297", "2001RW143", "2001RZ143",  "2001UN18", "2002CB225", "2002CD251", "2002CS154", "2002CT154",  "2002CU154", "2002CY154", "2002CY248", "2002CZ154", "2002CZ224",  "2002FX36", "2002PD155", "2002PO149", "2002PQ145", "2002VB131",  "2002VD131", "2002VE130", "2002VF131", "2002VT130", "2002XH91",  "2003FD128", "2003FK127", "2003FM127", "2003GF55", "2003LB7",  "2003LD9", "2003QA91", "2003QA92", "2003QB112", "2003QD91",  "2003QE112", "2003QE91", "2003QF113", "2003QF91", "2003QG91",  "2003QJ91", "2003QL91", "2003QN91", "2003QO91", "2003QQ91",  "2003QR91", "2003QS91", "2003QT91", "2003QU90", "2003QV90",  "2003QX90", "2003QY111", "2003QY90", "2003QZ111", "2003UN284",  "2003UN292", "2003UT291", "2003UY291", "2003UZ291", "2004DH64",  "2004DM64", "2004DM71", "2004DN64", "2004EP95", "2004ES95",  "2004EU95", "2004OL12", "2004PX107", "2004PY107", "2004UD10",  "2004VU75", "2005EC318", "2005EF298", "2005EM303", "2005EN302",  "2005EO304", "2005EX297", "2005GC187", "2005GX186", "2005JH177",  "2005JP179", "2005JR179", "2005JZ174", "275809", "307616", "53311",  "60454", "69987", "80806", "82157", "88267", "88268", "88611"]
allDES3to2 = ["119069", "119473", "133067", "139775", "169071", "1998UR43",  "1998WS31", "1998WV31", "1998WZ31", "2000CK105", "2001KB77",  "2001KD77", "2001KY76", "2001QF298", "2001QH298", "2001RU143",  "2001RX143", "2002CE251", "2002CW224", "2002GE32", "2002GF32",  "2002GL32", "2002GV32", "2002GW31", "2002GY32", "2002VD138",  "2002VX130", "2003FF128", "2003FL127", "2003QB91", "2003QH91",  "2003QX111", "2003UT292", "2003UV292", "2003WA191", "2004EH96",  "2004EJ96", "2004EV95", "2004VT75", "2004VZ75", "2005EZ300",  "2005GA187", "2005GB187", "2005GE187", "2005GF187", "2005GV210",  "28978", "306792", "307463", "69986", "69990", "91205"]
allDESscattered = ["118379", "118702", "134210", "138628", "148209", "181855", "181874",  "182934", "184212", "2000CG105", "2000CO105", "2000CQ105",  "2001FN185", "2001FT185", "2001FU185", "2001FV185", "2001KE77",  "2001KG77", "2001KO77", "2001KW76", "2001QA298", "2002CX224",  "2002GB32", "2002GH32", "2002VF130", "2003FH127", "2003FJ127",  "2003QA112", "2003QK91", "2003UB292", "2004DG77", "2004OJ14",  "2004OL14", "2004PA112", "2004PB108", "2004PS107", "2004PT107",  "2004PZ107", "2004TF282", "2005EF304", "2005JA175", "38083", "60458",  "82158"]
allDESresonant = ["118378", "119066", "119068", "119069", "119070", "119473",  "119878", "119956", "126619", "127871", "133067", "135024", "136120",  "139775", "149349", "160147", "160148", "169071", "182294", "182397",  "183964", "1998UR43", "1998UU43", "1998WS31", "1998WV31", "1998WZ31",  "1999HG12", "2000CK105", "2000OP67", "2000QL251", "2000QN251",  "2001FQ185", "2001KB77", "2001KD77", "2001KG76", "2001KL76",  "2001KY76", "2001QE298", "2001QF298", "2001QH298", "2001QR322",  "2001RU143", "2001RX143", "2001UP18", "2002CE251", "2002CW224",  "2002CZ248", "2002GD32", "2002GE32", "2002GF32", "2002GL32",  "2002GP32", "2002GS32", "2002GV32", "2002GW31", "2002GW32",  "2002GY32", "2002VD138", "2002VV130", "2002VX130", "2003FE128",  "2003FF128", "2003FL127", "2003LA7", "2003QB91", "2003QB92",  "2003QH91", "2003QW111", "2003QX111", "2003UT292", "2003UV292",  "2003WA191", "2004EG96", "2004EH96", "2004EJ96", "2004EV95",  "2004OK14", "2004PW107", "2004TT357", "2004TV357", "2004TX357",  "2004VK78", "2004VS75", "2004VT75", "2004VZ75", "2005EO297",  "2005ER318", "2005EZ300", "2005GA187", "2005GB187", "2005GE187",  "2005GF187", "2005GV210", "28978", "306792", "307463", "38084",  "42301", "69986", "69988", "69990", "91205", "95625"]

#runAllInputObjects(["2000CO105", "2002VT130", "2003QE91", "2003YU179", "2004UD10"])
#runAllInputObjects(allDESresonant)

### lazy  means I want to only calculate inc likelihoods for chip 2
# Grids are run lazy=True; real objects aren't (tiny difference: for 19521, P=0.07183 lazy and 0.07181 not)
runAllInputObjects(areWeLazy=False)
