#!/usr/bin/env python

import os
from dateutil import relativedelta
from datetime import datetime
from dateutil.parser import *
import scipy.io as sio
import numpy as np
import pandas as pd
import string

def cleanString(incomingString):
    # this is bad :( ... use regular expressions you fool
    newstring = incomingString
    newstring = newstring.replace("!","")
    newstring = newstring.replace("@","")
    newstring = newstring.replace("#","")
    newstring = newstring.replace("$","")
    newstring = newstring.replace("%","")
    newstring = newstring.replace("^","")
    newstring = newstring.replace("&","and")
    newstring = newstring.replace("*","")
    newstring = newstring.replace("(","")
    newstring = newstring.replace(")","")
    newstring = newstring.replace("+","")
    newstring = newstring.replace("=","")
    newstring = newstring.replace("?","")
    newstring = newstring.replace("\'","")
    newstring = newstring.replace("\"","")
    newstring = newstring.replace("{","")
    newstring = newstring.replace("}","")
    newstring = newstring.replace("[","")
    newstring = newstring.replace("]","")
    newstring = newstring.replace("<","")
    newstring = newstring.replace(">","")
    newstring = newstring.replace("~","")
    newstring = newstring.replace("`","")
    newstring = newstring.replace(":","")
    newstring = newstring.replace(";","")
    newstring = newstring.replace("|","")
    newstring = newstring.replace("\\","")
    newstring = newstring.replace("/","")        
    newstring = newstring.replace(" ","")     
    newstring = newstring.strip()
    return newstring

def getNumMonths(startDate, endDate):
    r = relativedelta.relativedelta(endDate, startDate)
    if r.years==0:
        numMonths = r.months
    if r.years>=1:
        numMonths = 12*r.years+r.months
    return numMonths

def parseGLTdatetimeStr(dtStr, maxDate):

    fmtNum = 0
    for fmt in ('%Y-%m-%d', '%m/%d/%y'):
        try:
            dtObj = datetime.strptime(dtStr, fmt)
            if dtObj > maxDate:
                dtObj -= relativedelta.relativedelta(years = 100)
            return (dtObj,fmtNum)
        except ValueError:
            fmtNum = fmtNum + 1
    raise ValueError('no valid date format found')

def normalizeGlobalLandTemperatures():
    globalLandTemperaturesFile = './raw_files/GlobalLandTemperaturesByCountry.csv'
    uniqueCountriesFile = './raw_files/UniqueCountriesGlobalLandTemperatures.txt'
    outputFile = './normalized_files/landTemperatures.mat'

    # only process if the output file doesn't exist already
    if(os.path.isfile(outputFile)):
        return

    maxDateFile = parse("2013-09-01")
    (endDate,fmtNum) = parseGLTdatetimeStr("2012-12-01",maxDateFile)
    
    f = open(uniqueCountriesFile)
    countries = f.readlines()
    f.close()
    numCountries = len(countries)
    printable = set(string.printable)

    gltfDict = {}
    countryBeingProcessed = ''
    with open(globalLandTemperaturesFile) as infile:
        for line in infile:
            lineSplit = line.split(',')
            if(lineSplit[0]=='dt'):
                # this is the header, we can ignore
                pass
            else:
                currDate = lineSplit[0]
                (currDateObj,fmtNum) = parseGLTdatetimeStr(currDate, maxDateFile)
                avgTemp = -999.0
                avgTempStr = lineSplit[1]
                if(avgTempStr!=''):
                    avgTemp = float(avgTempStr)
                currentCountry = lineSplit[3]
                # if we have any unfriendly characters in 'currentCountry', remove them
                currentCountry = filter(lambda x: x in printable, currentCountry)
                currentCountry = cleanString(currentCountry)

                if(currentCountry==countryBeingProcessed):
                    # this is the same country .., so we just insert the datapoint into the vector
                    # for this country
                    if(fmtNum==1 and currDateObj.year!=nextExpectedDate.year):
                        yearsDiff = relativedelta.relativedelta(currDateObj,nextExpectedDate)
                        currDateObj = currDateObj-relativedelta.relativedelta(years=yearsDiff.years)
                    nextExpectedDate = currDateObj+relativedelta.relativedelta(months=1)
                    if(currDateObj<=endDate):
                        idxToInsert = getNumMonths(startDate, currDateObj)
                        gltfDict[countryBeingProcessed][idxToInsert] = avgTemp
                else:
                    countryBeingProcessed = currentCountry 
                    
                    # make corrections ... b/c the datetime is ambiguously coded in this file :/
                    if(fmtNum==1 and currDateObj.year>2000):
                        currDateObj = currDateObj-relativedelta.relativedelta(years=100)

                    startDate = currDateObj
                    numMonths = getNumMonths(startDate, endDate)
                    gltfDict[countryBeingProcessed] = np.zeros(numMonths+1)
                    gltfDict[countryBeingProcessed][0] = avgTemp
                    print 'Processing ' + currentCountry + ' StartDate=' + str(startDate) + \
                          ' StopDate=' + str(endDate) + ' numMonths=' + str(numMonths)
                    nextExpectedDate = startDate+relativedelta.relativedelta(months=1)

    sio.savemat(outputFile, gltfDict)

def parseElninoDate(dateStr):
    dateObj = datetime.strptime(dateStr, '%y%m%d')
    # if the year is not in the 1900's ... 
    # then we mis-parsed so we correct for this
    if(dateObj.year>2000):
        dateObj = dateObj - relativedelta.relativedelta(years=100)

    return dateObj

def parseElninoNum(floatStr):
    try:
        return float(floatStr)
    except:
        return -999.0

def determineIfSameBuoy(prevDateObj, currDateObj):
    # this is the heuristic we use based on studying the data
    z = currDateObj-prevDateObj
    if(z.total_seconds()>0):
        return True
    else:
        return False

def normalizeElNino():
    elninoFile = './raw_files/elnino.csv'
    outputFile = './normalized_files/elnino.mat'

    # only process if the output file doesn't exist already
    if(os.path.isfile(outputFile)):
        return

    # define the subkeys, as they are static
    subkey_zonalWinds = 'zonalWinds'
    subkey_meridionalWinds = 'meridionalWinds'
    subkey_humidity = 'humidity'
    subkey_airTemp = 'airTemp'
    subkey_seaSurfaceTemp = 'seaSurfaceTemp'
    subkey_datecode = 'datecode'

    initNpArrSz = 1000

    elninoDict = {}
    buoyNumber = 0 
    with open(elninoFile) as infile:
        for line in infile:
            lineSplit = line.split(',')
            if(lineSplit[0]=='Observation'):
                # ignore the header
                pass
            else:
                zonalWindsVal = parseElninoNum(lineSplit[7])
                meridionalWindsVal = parseElninoNum(lineSplit[8])
                humidityVal = parseElninoNum(lineSplit[9])
                airTempVal = parseElninoNum(lineSplit[10])
                seaSurfaceTempVal = parseElninoNum(lineSplit[11])
                currDate = lineSplit[4]
                currDateObj = parseElninoDate(currDate)

                # if this is our first time
                if(buoyNumber==0):
                    prevDateObj = currDateObj
                sameBuoy = determineIfSameBuoy(prevDateObj,currDateObj)

                if(sameBuoy and buoyNumber>0):
                    # resize the array's if needed
                    if(numRecordsStored>=arrSize):
                        elninoDict[key][subkey_zonalWinds] = np.resize(elninoDict[key][subkey_zonalWinds], arrSize+initNpArrSz)
                        elninoDict[key][subkey_meridionalWinds] = np.resize(elninoDict[key][subkey_meridionalWinds], arrSize+initNpArrSz)
                        elninoDict[key][subkey_humidity] = np.resize(elninoDict[key][subkey_humidity], arrSize+initNpArrSz)
                        elninoDict[key][subkey_airTemp] = np.resize(elninoDict[key][subkey_airTemp], arrSize+initNpArrSz)
                        elninoDict[key][subkey_seaSurfaceTemp] = np.resize(elninoDict[key][subkey_seaSurfaceTemp], arrSize+initNpArrSz)
                        elninoDict[key][subkey_datecode] = np.resize(elninoDict[key][subkey_datecode], arrSize+initNpArrSz)
                        arrSize = arrSize + initNpArrSz

                    elninoDict[key][subkey_zonalWinds][numRecordsStored] = zonalWindsVal
                    elninoDict[key][subkey_meridionalWinds][numRecordsStored] = meridionalWindsVal
                    elninoDict[key][subkey_humidity][numRecordsStored] = humidityVal
                    elninoDict[key][subkey_airTemp][numRecordsStored] = airTempVal
                    elninoDict[key][subkey_seaSurfaceTemp][numRecordsStored] = seaSurfaceTempVal
                    elninoDict[key][subkey_datecode][numRecordsStored] = (currDateObj-datetime(1980,1,1)).days
                    numRecordsStored = numRecordsStored + 1
                else:
                    if(buoyNumber!=0):
                        # lets truncate the numpy array's to the amount of records we have
                        # to make post-processing easier :)
                        elninoDict[key][subkey_zonalWinds] = np.resize(elninoDict[key][subkey_zonalWinds], numRecordsStored)
                        elninoDict[key][subkey_meridionalWinds] = np.resize(elninoDict[key][subkey_meridionalWinds], numRecordsStored)
                        elninoDict[key][subkey_humidity] = np.resize(elninoDict[key][subkey_humidity], numRecordsStored)
                        elninoDict[key][subkey_airTemp] = np.resize(elninoDict[key][subkey_airTemp], numRecordsStored)
                        elninoDict[key][subkey_seaSurfaceTemp] = np.resize(elninoDict[key][subkey_seaSurfaceTemp], numRecordsStored)
                        elninoDict[key][subkey_datecode] = np.resize(elninoDict[key][subkey_datecode], numRecordsStored)

                    buoyNumber = buoyNumber + 1
                    # we are recording a new buoy's measurements, setup the dictionary
                    key = 'buoy%d'%buoyNumber
                    zonalWindsArr = np.zeros(initNpArrSz)
                    meridionalWindsArr = np.zeros(initNpArrSz)
                    humidityArr = np.zeros(initNpArrSz)
                    airTempArr = np.zeros(initNpArrSz)
                    seaSurfaceTempArr = np.zeros(initNpArrSz)
                    datecodeArr = np.zeros(initNpArrSz)

                    print 'Processing ' + key

                    elninoDict[key] = {subkey_zonalWinds:zonalWindsArr, \
                                       subkey_meridionalWinds:meridionalWindsArr, \
                                       subkey_humidity:humidityArr, \
                                       subkey_airTemp:airTempArr, \
                                       subkey_seaSurfaceTemp:seaSurfaceTempArr, \
                                       subkey_datecode:datecodeArr}

                    elninoDict[key][subkey_zonalWinds][0] = zonalWindsVal
                    elninoDict[key][subkey_meridionalWinds][0] = meridionalWindsVal
                    elninoDict[key][subkey_humidity][0] = humidityVal
                    elninoDict[key][subkey_airTemp][0] = airTempVal
                    elninoDict[key][subkey_seaSurfaceTemp][0] = seaSurfaceTempVal
                    elninoDict[key][subkey_datecode][0] = (currDateObj-datetime(1980,1,1)).days

                    numRecordsStored = 1
                    arrSize = initNpArrSz
                
                prevDateObj = currDateObj

    sio.savemat(outputFile, elninoDict)


def normalizePollution():
    pollutionFile = './raw_files/pollution_us_2000_2016.csv'
    outputFile = './normalized_files/pollution.csv'

    # only process if the output file doesn't exist already
    if(os.path.isfile(outputFile)):
        return

    data = pd.read_csv(pollutionFile)
    data = data.drop(["Unnamed: 0", "State Code", "County Code", \
                      "State", "City", "County", "Address", "NO2 Units", \
                      "O3 Units", "SO2 Units", "CO Units"],axis=1)
    aqiData = data[['Site Num', 'Date Local', 'NO2 AQI', 'O3 AQI', 'SO2 AQI', 'CO AQI']]
    aqiData = aqiData.dropna(axis="rows")
    aqiData['Date Local'] = (pd.to_datetime(aqiData['Date Local'],format='%Y-%m-%d')-datetime(1980,1,1)).dt.days
    aqiData = aqiData.groupby(['Site Num','Date Local']).mean()
    
    aqiData.to_csv(outputFile)


if __name__=='__main__':
    normalizeGlobalLandTemperatures()
    normalizeElNino()
    normalizePollution()