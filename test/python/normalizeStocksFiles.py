#!/usr/bin/env python

import glob
import pandas as pd
import os

def normalizeFile(fileToProcess):

    fname = os.path.basename(fileToProcess)
    toFile = os.path.join('./normalized_files',fname)

    df = pd.read_csv(fileToProcess)

    if(fname=='ftse.csv'):
        df['Date'] = pd.to_datetime(df['Date'],format='%Y%m%d')
    else:
        df['Date'] = pd.to_datetime(df['Date'])
    minDate = df.iloc[len(df)-1]['Date']
    maxDate = df.iloc[0]['Date']
    idx = pd.date_range(minDate, maxDate)

    df.set_index(['Date'],inplace=True)
    df = df.reindex(idx, method='nearest')
    df.to_csv(toFile)

if __name__=='__main__':
    # list all the files to be checked
    
    
    filelist = glob.glob('./raw_files/*.csv')
    for fileToProcess in filelist:
        print('Processing -- ' + fileToProcess)
        normalizeFile(fileToProcess)
    