#!/usr/bin/env python

import os
import sys
import glob
import csv

# Converts a GCT to a CSV file for easier analysis.  The specification
# for the GCT file is given here:
# http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT

def convertFile(gctFile, csvFile):
    print 'Converting ' + gctFile + ' to ' + csvFile
    with open(gctFile) as f:
        lines = f.readlines()
    # according to the file format, the data starts on Line 4, after col 2

    writer = csv.writer(open(csvFile, 'w'),lineterminator='\n')
    for ii in range(3,len(lines)):
        line = lines[ii].strip()
        lineSplit = line.split('\t')
        writer.writerow(lineSplit[2:])

if __name__=='__main__':
    dataDir = sys.argv[1]
    gctFiles = glob.glob(os.path.join(dataDir, 'gct_files', '*.gct'))
    for gctFile in gctFiles:
        fname = os.path.basename(gctFile)
        fnameWithoutExt = os.path.splitext(fname)[0]
        csvFile = os.path.join(dataDir, 'csv_files', fnameWithoutExt) + '.csv'
        convertFile(gctFile, csvFile)