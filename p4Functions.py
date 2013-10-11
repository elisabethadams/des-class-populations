#!/usr/bin/env python
#
# Python version of the relevant p4 Mathematica functions developed for the jleGroup suite of packages
# Mathematica version by A. A. S. Gulbis (last modified 2010) first published in
#   Gulbis, A. A. S. et al. 2010. Unbiased Inclination Distributions for Objects in the Kuiper Belt. AJ, 140, 350
# Python translation/tests by E. R. Adams (June 2012)


import sys
import os
import string
import math
import numpy as np
### quad requires Fortran! Does your mac come with Fortran? Go install it if not ###
from scipy.integrate import quad

### Single DES chip height and width, in arcsec
p4ChipHeightAndWidth = [0.29582, 0.147911]

degree = 2.0 * math.pi / 360.

### Turn on/off diagnostic tests at end
diag=False


### Subfunctions called by the main dealy-job
def pBi(inc, lat):
	if ((inc == 0.) & (lat == 0.)):
		return 1
	elif math.sin(inc) <= abs(math.sin(lat)):
		return 0
	else:
 		return math.cos(lat)/(math.pi*math.sqrt(math.sin(inc)**2 - math.sin(lat)**2))

def delLambMaxFcn(thetaRad):
	if (math.tan(thetaRad) >= (p4ChipHeightAndWidth[0]/p4ChipHeightAndWidth[1])  ):
		return p4ChipHeightAndWidth[0]*degree/math.sin(thetaRad)
	else:
		return p4ChipHeightAndWidth[1]*degree/math.cos(thetaRad)


def deltaLambda(betaRad, betaminRad, betamaxRad, thetaRad): 
	if (math.tan(thetaRad) >= (p4ChipHeightAndWidth[0]/ p4ChipHeightAndWidth[1])):
		beta1 = p4ChipHeightAndWidth[0] * degree * math.cos(thetaRad)
	else:
		beta1 = p4ChipHeightAndWidth[1] * degree * math.sin(thetaRad)
	deltaLambdaMax = delLambMaxFcn(thetaRad) 
	if -(math.pi/2) <= betaRad < betaminRad:
		return 0.
	elif betamaxRad < betaRad <= math.pi/2:
	 	return 0.
	elif betaminRad <= betaRad < beta1 + betaminRad:
		return deltaLambdaMax * (betaRad - betaminRad)/beta1
	elif beta1 + betaminRad <= betaRad <= betamaxRad - beta1: 
	 	return deltaLambdaMax
	elif betamaxRad - beta1 <= betaRad <= betamaxRad:
		return deltaLambdaMax * (betamaxRad - betaRad)/beta1

def p4DeltaLambdaMax(thetaRad):
	if (math.tan(thetaRad) >= p4ChipHeightAndWidth[0]/p4ChipHeightAndWidth[1]):
		return p4ChipHeightAndWidth[0]*degree/math.sin(thetaRad)
	else:
		return p4ChipHeightAndWidth[1]*degree/math.cos(thetaRad)

## Use this form for numerical integration
def integrateDeltaLambdaTimesPBi(betaRad, betaminRad, betamaxRad, thetaRad, inc):
	return deltaLambda(betaRad, betaminRad, betamaxRad, thetaRad) * pBi(inc, betaRad)



########## Finding the normalization function #############		
## This number was derived in Mathematica using (NIntegrate form preferred):
##(*normFactor=Re[Integrate[deltaLambda[lat,-90 *degree,90 [Degree],0.] pBi[1. *degree,lat],{lat,-90 *degree,90 [Degree]}]];*)
## normFactor=Re[NIntegrate[deltaLambda[lat,-90 *degree,90 *degree,0.] pBi[1. *degree,lat],{lat,-90 \*degree,90 *degree},MaxRecursion->60]]; 
normFactor = 0.0025815339499173327

## Note that we can get the same result in python! (good...)
##  however it requires splitting the integral at zero to avoid a singularity
#	res1, err1 = quad(integrateDeltaLambdaTimesPBi, -90. * degree, 0. * degree, args=(-90. *degree, 90.* degree, 0., 1.* degree))
#	res2, err2 = quad(integrateDeltaLambdaTimesPBi, 0. * degree, 90. * degree, args=(-90. *degree, 90.* degree, 0., 1.* degree))
#	res = res1 + res2
#	print res
##########################################################


### The main dealy-job, the inclination likelihood factor
def p4RelativeFieldLikelihood(thetaDeg, midLatDeg, incDeg): 
	b0 = abs(midLatDeg) 
	thetaRad = round(thetaDeg * degree, 8) 
	inc = round(incDeg*degree, 8) 
	fh = p4ChipHeightAndWidth[0] * degree 
	fw = p4ChipHeightAndWidth[1] * degree 
	h = fh * math.cos(thetaRad) + fw * math.sin(thetaRad) 
	if (math.tan(thetaRad) >= (fh/fw)):
		deltaB1 = fh * math.cos(thetaRad)
	else:
		deltaB1 = fw * math.sin(thetaRad) 
	deltaB2 = 0.5*(h-2*deltaB1)
	halfHeight = deltaB1 + deltaB2
	latMax = round(b0 *degree + halfHeight, 8)
	latMin = round(b0 *degree - halfHeight, 8)
	edge1 = round(b0 *degree - deltaB2, 8)
	edge2 = round(b0 *degree + deltaB2, 8)
	if (latMin < 0 <= edge1):
		step1 = 0; step2 = edge1; step3 = edge2;
	elif (edge1 < 0 <= edge2):
		step1 = edge1; step2 = 0; step3 = edge2;
	elif (edge2 < 0 <= latMax):
		step1 = edge1; step2 = edge2; step3 = 0;
	if ((inc == 0.) & (latMin < inc < latMax)):
		return 1./math.cos(thetaRad)
	elif ((thetaDeg == 0.) | (thetaDeg == 90.)): 
	 	return p4RelativeFieldLikelihoodBruteForce(thetaDeg, b0, incDeg), 
	elif latMin < 0 < latMax:
		res1, err1 = quad(integrateDeltaLambdaTimesPBi, latMin, step1, args=(latMin, latMax, thetaRad, inc))
		res2, err2 = quad(integrateDeltaLambdaTimesPBi, step1, step2, args=(latMin, latMax, thetaRad, inc))
		res3, err3 = quad(integrateDeltaLambdaTimesPBi, step2, step3, args=(latMin, latMax, thetaRad, inc))
		res4, err4 = quad(integrateDeltaLambdaTimesPBi, step3, latMax, args=(latMin, latMax, thetaRad, inc))
		return 1/normFactor *  (res1 + res2 + res3 + res4)
	else: 
		res1, err1 = quad(integrateDeltaLambdaTimesPBi, latMin, edge1, args=(latMin, latMax, thetaRad, inc))
		res2, err2 = quad(integrateDeltaLambdaTimesPBi, edge1, edge2, args=(latMin, latMax, thetaRad, inc))
		res3, err3 = quad(integrateDeltaLambdaTimesPBi, edge2, latMax, args=(latMin, latMax, thetaRad, inc))
		return 1/normFactor *  (res1 + res2 + res3)

	
### and the brute force version (noy as accurate, but faster):
## First some helpers
def p4EfficiencyForPeakingOrbits(tiltAngle, eclipticLatitude, inclination):
	b0 = eclipticLatitude
	i = inclination
	theta = tiltAngle
	fh = p4ChipHeightAndWidth[0] * degree
	fw = p4ChipHeightAndWidth[1] * degree
	h = fh*math.cos(theta) + fw*math.sin(theta)
	
	if ((b0 < h/2) & (i > b0 + h/2) & (i > abs(b0 - h/2))): 
		lowEnd = abs(b0 - h/2)
	elif b0 > h/2:
		lowEnd = abs(b0 - h/2)
	else:
		lowEnd = 0
		
	if b0 <= h/2:
		return 0.5
	else:
		return 0.5 - 1 / math.pi * math.arctan(1/math.sin(i) * math.sin(lowEnd) / math.sqrt(1 - (1/math.sin(i) * math.sin(lowEnd))**2))


def p4EfficiencyForCrossingOrbits(tiltAngle, eclipticLatitude, inclination):
	b0 = abs(eclipticLatitude) # (* note that the ecliptic lat. read into this fcn. is the absolute value*)
	i = inclination
	theta = tiltAngle
	fh = p4ChipHeightAndWidth[0] * degree
	fw = p4ChipHeightAndWidth[1] * degree
	h = fh*math.cos(theta) + fw*math.sin(theta)
	if ((b0 < h/2) & (i <= max([b0 + h/2, abs(b0 - h/2)])) ):
		lowEnd = 0
		highEnd = min([abs(b0 + h/2), abs(b0 - h/2)])
	else:
		lowEnd = b0 - h/2 
		highEnd = b0 + h/2
	return (fh/h * 1/math.pi*(math.asin(1/math.sin(i) * math.sin(highEnd)) - math.asin(1/math.sin(i) * math.sin(lowEnd))))


## this one is big!
def p4TiltedProportion(tiltAngle, eclipticLatitude, inclination, region):
	b0 = abs(eclipticLatitude)
	i = abs(inclination)
	theta = tiltAngle
	fh = p4ChipHeightAndWidth[0] * degree
	fw = p4ChipHeightAndWidth[1] * degree
	h = fh*math.cos(theta) + fw*math.sin(theta)
	if (math.tan(theta) >= (fh/fw)):
		deltaB1 = fh * math.cos(theta)
	else:
		deltaB1 = fw * math.sin(theta) 
	deltaB2 = 1/2 (h - 2*deltaB1)
	if i == 0:
		area = fw
	else:
		area = fw*i
	allSky =  1/((i - (b0 - h/2))*fw) 
	
	## Define nonEclipCross
	if tiltAngle == 0.:  # (*frame is not tilted so prob. is one - without this the trig functions explode,fw*(i-(b0-h/2))* allSky *)
		nonEclipCross = 1.
	elif tiltAngle == 90.: # (*frame is tilted fully so prob. is fraction difference between frame height and width - without this the trig functions explode fh*(i-(b0-h/2))*allSky
		nonEclipCross = fh/fw
	elif ((b0 - h/2) <= i & (i <= b0 - deltaB2)): # (*orbit peaks below parallelogram*)
		nonEclipCross = ((i - (b0 - h/2))**2 / (2*math.sin(theta)*math.cos(theta)))*allSky,
	elif ((b0 - deltaB2) < i & (i <= (b0 + deltaB2))): #(*orbit peaks in parallelogram*)
		nonEclipCross = ((-deltaB2 + h/2)**2 / (2*math.cos(theta)*math.sin(theta)) + (i - (b0 - deltaB2))*p4DeltaLambdaMax(theta))*allSky
	elif ((b0 + deltaB2) < i & (i <= b0 + h/2)): #(*orbit peaks above parallelogram*)
		nonEclipCross = (fw*fh - (b0 + h/2 - i)**2 / (2*math.cos(theta)*math.sin(theta)))*allSky
	
	## Define negSide
	if tiltAngle == 0.:
		negSide = 1 # (*geometric terms blow up at 0 and Pi *)
	elif tiltAngle == 90.:
		negSide = fh/fw
	elif ((i == 0) & (b0 <= deltaB2)):
		negSide = p4DeltaLambdaMax(theta)/fw # (* if we don't account for i==0, then it's always 0 *)
	elif ((i == 0) & (b0 > deltaB2)):
		negSide = (2*abs(b0 - h/2) * math.csc(2*theta))/fw
	elif ((b0 > deltaB2) & (i <= abs(b0 - h/2))):
		negSide = (2*i*abs(b0 - h/2) - i**2) / (2*math.cos(theta) * math.sin(theta)) / area	
	elif ((b0 > deltaB2) & (abs(b0 - h/2) < i)):
		negSide = (abs(b0 - h/2) / (2 * math.sin(theta) * math.cos(theta)))/fw
	elif ((b0 <= deltaB2) & (i <= abs(b0 - deltaB))):
		negSide = p4DeltaLambdaMax(theta)*i/area
	elif ((b0 <= deltaB2) & (abs(b0 - deltaB2) < i <= abs(b0 - h/2))):
		negSide = (p4DeltaLambdaMax(theta)*abs(b0-deltaB2) + ((abs(b0-h/2)-abs(b0-deltaB2))**2 - (abs(b0-h/2)-i)**2)/(2*math.sin(theta) * math.cos(theta)))/area
	elif ((b0 <= deltaB2) & (i > abs(b0 - h/2))):
		negSide = (p4DeltaLambdaMax(theta)*abs(b0-deltaB2) + ((abs(b0-h/2)-abs(b0-deltaB2))**2) / (2*math.sin(theta)*math.cos(theta))) / (fw*abs(b0-h/2))
		
	## Define posSide
	if tiltAngle == 0.:
		posSide = 1
	elif tiltAngle == 90.:
		posSide = fh/fw
	elif ((b0 <= deltaB2) & (i <= abs(b0 + deltaB2))): 
		posSide = i*p4DeltaLambdaMax(theta)/area
	elif ((b0 <= deltaB2) & abs(b0 + deltaB2) < i <= abs(b0 + h/2)):
		posSide = (p4DeltaLambdaMax(theta)*(b0 + deltaB2) + (((b0 + h/2) - (b0 + deltaB2))**2 - ((b0 + h/2) - i)**2)/(2*math.sin(theta) * math.cos(theta)))/area
	elif ((b0 <= deltaB2) & (i > abs(b0 + h/2))):
		posSide = (p4DeltaLambdaMax(theta)*(b0 + deltaB2) + (((b0 + h/2) - (b0 + deltaB2))**2)/(2*math.sin(theta) *math.cos(theta))) / (fw*(b0 + h/2))
	elif ((b0 > deltaB2) & (i <= abs(b0 - deltaB2))):
		posSide = ((abs(b0 - h/2) + i)**2 - (abs(b0 - h/2))**2)/(2*math.cos(theta)*math.sin(theta)*area)
	elif ((b0 > deltaB2) & (abs(b0 - deltaB2) < i <= abs(b0 + deltaB2))):
		posSide = (p4DeltaLambdaMax(theta)*(i-(b0-deltaB2)) + (((b0 - deltaB2) + abs(b0 - h/2))**2 - (abs(b0 - h/2))**2) / (2*math.sin(theta)*math.cos(theta))) / area
	elif ((b0 > deltaB2) & (abs(b0 + deltaB2) < i <= abs(b0 + h/2))):
		posSide = (fh*fw - ((abs(b0 - h/2))**2 + ((b0 + h/2) - i)**2) / (2*math.sin(theta)*math.cos(theta))) / area
	elif ((b0 > deltaB2) & (abs(b0 + h/2) < i)):
		posSide = (fh*fw - (abs(b0-h/2))**2 / (2*math.sin(theta) * math.cos(theta))) / (fw*(b0 + h/2))
	### Finally done!
	if b0 > h/2:
		return nonEclipCross
	elif region == "negative":
		return negSide
	elif region == "positive":
		return posSide
	elif region == "both":
		return negSide + posSide


#### The brute force dude
def p4RelativeFieldLikelihoodBruteForce(tiltAngle, eclipticLatitude, inclination):
	b0 = abs(eclipticLatitude) * degree
	rawb=eclipticLatitude * degree
	i = abs(inclination) * degree
	fh = p4ChipHeightAndWidth[0] * degree
	fw=p4ChipHeightAndWidth[1] * degree
	theta = tiltAngle * degree
	h = fh*math.cos(theta) + fw*math.sin(theta)
	if math.tan(theta)>=fh/fw:
		deltaB1 = fh*math.cos(theta)
	else:
		deltaB1 = fw*math.sin(theta)
	deltaB2 = (h-2*deltaB1)/2
	
	if i <= b0-h/2: #	(* inclination never within range of latitudes *)
		return 0.0, 
	elif ((b0<=deltaB2) & (i<=abs(b0-deltaB2)) & (i<=b0+deltaB2)): # (*inclination always within range of latitudes*)
		return 1./math.cos(theta)
	elif ((b0 <= h/2) & (i <= min([b0+h/2,abs(b0-h/2)]))): # (*orbit peaks in tilted, ecliptic-crossing frame *)
		return p4EfficiencyForPeakingOrbits(theta,b0,i) * p4TiltedProportion(theta,b0,i,"both")
	elif ((b0 <= h/2) & (b0+h/2>=i>=abs(b0-h/2))): # (*orbit peaks in positive and crosses  in negative ecliptic-crossing frame *)
		return (p4EfficiencyForPeakingOrbits(theta,b0,i) * p4TiltedProportion(theta,b0,i,"positive") + p4EfficiencyForPeakingOrbits(theta,b0,i) * p4TiltedProportion(theta,b0,i,"negative"))
	elif ((b0 <= h/2) &  (abs(b0-h/2)>=i>=b0+h/2)): # (*orbit peaks negative and crosses  in positive ecliptic-crossing frame *)
		return (p4EfficiencyForPeakingOrbits(theta,b0,i) *p4TiltedProportion(theta,b0,i,"negative") + p4EfficiencyForCrossingOrbits(theta,b0,i) * p4TiltedProportion(theta,b0,"positive"))
	elif ((b0 <= h/2) & (i > max([abs(rawb + h/2), abs(rawb-h/2)]))): #, (*orbit crosses completely through ecliptic-crossing frames*)
		return (p4EfficiencyForCrossingOrbits(theta,b0,i) * p4TiltedProportion(theta,b0,i,"both") / 2)
	elif ((b0 > h/2) & ((b0-h/2)<i<(b0+h/2))): # (* frame does not cross ecliptic, and orbit has peak latitude within range *)
		return (p4EfficiencyForPeakingOrbits(theta,b0,i) * p4TiltedProportion(theta,b0,i,"both"))
	elif ((b0 > h/2) & ((b0+h/2)<=i)): #	(* frame does not cross ecliptic, orbit crosses range of latitudes *)
		return p4EfficiencyForCrossingOrbits(theta,b0,i)



########## DIAGNOSTICS ##########
if diag:
	print "\nDiagnostics\n"
	from time import clock, time
	before=time()
	test=p4RelativeFieldLikelihood(12.,3.4,23.4)
	after = time()
	print "Inc likelihood:",test," which took ", round(after-before,8)," seconds\n"	
	before=time()
	test=p4RelativeFieldLikelihoodBruteForce(12.,3.4,23.4)
	after = time()
	print "...brute force:",test," which took ", round(after-before,8)," seconds"
	print "Inclination angle test:"
	testTilt = 8.2
	testLat = 0.2
	print "Test tilt",testTilt,"test lat",testLat
	for ii in [0.0,0.01,0.1,1.,2.,5.,8.,9.]:
		test=p4RelativeFieldLikelihood(testTilt, testLat, ii)
		print ii, test
	for ii in np.arange(10.,180.,10):
		test=p4RelativeFieldLikelihood(testTilt, testLat, ii)
		print ii, test
	for ii in [171,175,179,179.9,179.99,180.0]:
		test=p4RelativeFieldLikelihood(testTilt, testLat, ii)
		print ii, test