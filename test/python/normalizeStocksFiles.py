#!/usr/bin/env python

import glob
from datetime import datetime
from datetime import timedelta
from dateutil.parser import *
import os

def isValidDay(dateObj):
    isWeekday = dateObj.isoweekday() in range(1, 6)
    
    isChristmas = (dateObj.month==12) and (dateObj.day==25)
    isNewYears = (dateObj.month==1) and (dateObj.day==1)
    
    isSpecialDay = isChristmas or isNewYears
    
    if(isWeekday):
        if(isSpecialDay):
            return False
        else:
            return True
    else:
        return False
    

def getNextValidDay(dateObj):
    while(True):
        dateObj = dateObj - timedelta(days=1)
        if(isValidDay(dateObj)):
            break

    return dateObj

def normalizeFile(fileToProcess):
    expectedDate = '2017-01-06'
    expectedDateObj = parse(expectedDate)

    f = open(fileToProcess)
    print 'Processing file=' + str(fileToProcess)
    
    # create filename etc... for the normalized file
    fromFileDir = os.path.dirname(fileToProcess)
    fname = os.path.basename(fileToProcess)
    toFileDir = os.path.join(fromFileDir, '../normalized_files')
    toFile = os.path.join(toFileDir,fname)
    fOut = open(toFile, 'w')
    
    z = f.readlines()
    fOut.write(z[0])        # keep the header
    
    missingDates = []
    prevData = []
    ii = 1
    while ii<len(z):
        line = z[ii]
        data = line.split(',')
        date = data[0]
        dateObj = parse(date)

        delta = expectedDateObj-dateObj
        deltaDays = delta.days
        if(deltaDays==0):
            # this is a valid record, we will print it
            fOut.write(line)
            
            # get the next valid day
            expectedDateObj = getNextValidDay(expectedDateObj)
        else:
            if(deltaDays>0):
                missingDates.append(expectedDateObj)
                #print 'Missing Day=' + str(expectedDateObj)
                # NOTE: we are doing a sample&hold of the data ... is this ideal?
                tmpImputedData = prevData
                tmpImputedData[0] = datetime.strftime(expectedDateObj,'%Y-%m-%d')    # this format might not match?
                fOut.write(','.join(tmpImputedData))
                
                while(True):
                    expectedDateObj = expectedDateObj - timedelta(days=1)
                    delta = expectedDateObj-dateObj
                    deltaDays = delta.days
                    
                    if(deltaDays==0):
                        if(isValidDay(expectedDateObj)):
                            # this is a valid record, add it in to the data that is "normalized"
                            fOut.write(line)
                        expectedDateObj = getNextValidDay(expectedDateObj)
                        break
                    else:
                        # if it is a valid date, mark it as a missing date
                        if(isValidDay(expectedDateObj)):
                            missingDates.append(expectedDateObj)
                            #print 'Missing Day=' + str(expectedDateObj)
                            # NOTE: we are doing a sample&hold of the data ... is this ideal?
                            tmpImputedData = prevData
                            tmpImputedData[0] = datetime.strftime(expectedDateObj,'%Y-%m-%d')    # this format might not match?
                            fOut.write(','.join(tmpImputedData))
                    
            else:
                # means a date exists in the file that is not considered a "valid" date
                # by our isValidDay function.  We impute this data to align with the rest
                # of the files
                
                #print 'Fast Forwarding >>>>'
                #print 'Expected Date = ' + str(expectedDateObj) + ' actualDate = ' + str(dateObj)
                                
                # fast-foward the file to catch-up to the next expected date
                while True:
                    ii = ii + 1
                    if(ii>len(z)):
                        break
                    line = z[ii]
                    data = line.split(',')
                    date = data[0]
                    dateObj = parse(date)
                    #print 'dateObj = ' + str(dateObj)
                    
                    if(dateObj<expectedDateObj):
                        # fill in the missing values between dateObj and expectedDateObj
                        while(expectedDateObj>=(dateObj+timedelta(days=1))):
                            # NOTE: we are doing a sample&hold of the data ... is this ideal?
                            if(isValidDay(expectedDateObj)):
                                #print 'Imputing day=' + str(expectedDateObj)
                                tmpImputedData = prevData
                                tmpImputedData[0] = datetime.strftime(expectedDateObj,'%Y-%m-%d')    # this format might not match?
                                fOut.write(','.join(tmpImputedData))
                            expectedDateObj = expectedDateObj - timedelta(days=1)
                        
                        # find the next valid date that is lower than or equal to this date
                        expectedDateObj = dateObj
                        if(isValidDay(expectedDateObj)):
                            break
                        else:
                            expectedDateObj = getNextValidDay(expectedDateObj)
                    elif(dateObj==expectedDateObj):
                        break
                    
                ii = ii - 1     # we do this, b/c the +1 below, inefficient :(
            
        ii = ii + 1
        prevData = data

    fOut.close()
    
def testNormalization(folder):
    # makes sure that all the files have the same dates, which means that
    # they are all normalized
    filelist = glob.glob(os.path.join(folder, '*.csv'))
    numFiles = len(filelist)
        
    # this is super inefficient, but its easy :/
    for ii in range(0,numFiles):
        for jj in range(ii+1,numFiles):
            print 'Comparing ' + str(filelist[ii]) + ' to ' + str(filelist[jj])
            f1 = open(filelist[ii])
            f2 = open(filelist[jj])
            
            f1Data = f1.readlines()
            f2Data = f2.readlines()
            
            numLinesToCheck = min(len(f1Data),len(f2Data))
            for kk in range(1,numLinesToCheck):
                f1Line = f1Data[kk]
                f2Line = f2Data[kk]
                f1Date = f1Line.split(',')[0]
                f2Date = f2Line.split(',')[0]
                f1DateObj = parse(f1Date)
                f2DateObj = parse(f2Date)
                if(f1DateObj!=f2DateObj):
                    print '*** ERROR ***'
                    print f1DateObj
                    print f2DateObj
                    print '*************'
                    raise Exception('ERROR')
                    
            
            print '*** Passed ***'
    

if __name__=='__main__':
    # list all the files to be checked
    
    filelist = glob.glob('./raw_files/*.csv')
    for fileToProcess in filelist:
        normalizeFile(fileToProcess)
    
        
    testNormalization('./normalized_files')