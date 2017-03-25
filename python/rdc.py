#!/usr/bin/env python

import warnings

import numpy as np
import scipy.stats

def pobs(x):
    if(type(x)!=np.ndarray):
        raise Exception('Input x must be a numpy ndarray datatype!')

    shape_x = x.shape
    if(len(shape_x)<2):
        x = x.reshape((shape_x[0],1))
    n = x.shape[0]
    d = x.shape[1]

    u = np.zeros(shape=x.shape)
    for ii in range(d):
        u[:,ii] = scipy.stats.rankdata(x)/(n+1.)

    return u

def cca(X, Y):
    """
    Canonical Correlation Analysis
    Currently only returns the canonical correlations.
    """
    n,p1 = X.shape
    n,p2 = Y.shape

    # center X and Y
    meanX = X.mean(axis=0)
    meanY = Y.mean(axis=0)
    X = X-meanX[np.newaxis,:]
    Y = Y-meanY[np.newaxis,:]

    Qx,Rx = np.linalg.qr(X)
    Qy,Ry = np.linalg.qr(Y)

    rankX = np.linalg.matrix_rank(Rx)
    if rankX == 0:
        raise Exception('Rank(X) = 0! Bad Data!')
    elif rankX < p1:
        #warnings.warn("X not full rank!")
        Qx = Qx[:,0:rankX]  
        Rx = Rx[0:rankX,0:rankX]

    rankY = np.linalg.matrix_rank(Ry)
    if rankY == 0:
        raise Exception('Rank(X) = 0! Bad Data!')
    elif rankY < p2:
        #warnings.warn("Y not full rank!")
        Qy = Qy[:,0:rankY]
        Ry = Ry[0:rankY,0:rankY]

    d = min(rankX,rankY)
    svdInput = np.dot(Qx.T,Qy)

    U,r,V = np.linalg.svd(svdInput)
    r = np.clip(r,0,1)
    #A = np.linalg.lstsq(Rx, U[:,0:d]) * np.sqrt(n-1)
    #B = np.linalg.lstsq(Ry, V[:,0:d]) * np.sqrt(n-1)
    
    # TODO: resize A to match inputs

    #return (A,B,r)
    return r

def rdc(x,y,k=20,s=1./6.):
    """
    Computes the RDC between two sets of (possibly multivariate)
    random variables.  x and y are matrices, both with n samples, and
    d1 and d2 (possibly w/ d1=d2) dimensions.  They must be numpy 
    datatypes!

    See: https://papers.nips.cc/paper/5138-the-randomized-dependence-coefficient.pdf
    """
    if(type(x)!=np.ndarray or type(y)!=np.ndarray):
        raise Exception('Inputs x and y must be numpy ndarray datatypes!')
    shape_x = x.shape
    shape_y = y.shape
    if( len(shape_x)>2 or len(shape_y)>2 ):
        raise Exception('Tensor inputs not allowed!')
    if( len(shape_x)==1 ):
        shape_x = (shape_x[0], 1)
        x = x.reshape((shape_x[0],1))
    if( len(shape_y)==1 ):
        shape_y = (shape_y[0], 1)
        y = y.reshape((shape_y[0],1))
    # ensure we have the same # of observations for x and y
    if(shape_x[0]!=shape_y[0]):
        raise Exception('The number of observations for RVs x and y must be the same!')

    # compute pseudo-observations & insert into CCA computation vector
    xx = np.ones(shape=(shape_x[0],shape_x[1]+1))
    yy = np.ones(shape=(shape_y[0],shape_y[1]+1))
    u = pobs(x)
    v = pobs(y)
    # copy the data into xx & yy
    for ii in range(shape_x[1]):
        xx[:,ii] = u[:,ii]
    for ii in range(shape_y[1]):
        yy[:,ii] = v[:,ii]

    xn = np.random.normal(size=(xx.shape[1],k))
    yn = np.random.normal(size=(yy.shape[1],k))

    xs = np.sin(s/xx.shape[1]*np.dot(xx,xn))    # dot is matrix multiplication for ndarray datatypes
    ys = np.sin(s/yy.shape[1]*np.dot(yy,yn))    # dot is matrix multiplication for ndarray datatypes

    xu = np.ones(shape=(xs.shape[0],xs.shape[1]+1))
    yu = np.ones(shape=(ys.shape[0],ys.shape[1]+1))

    # copy the data into xu and yu
    for ii in range(xs.shape[1]):
        xu[:,ii] = xs[:,ii]
    for ii in range(ys.shape[1]):
        yu[:,ii] = ys[:,ii]

    r = cca(xu, yu)
    res = r[0]
    return res

if __name__=='__main__':
    import scipy.io as sio

    np.set_printoptions(precision=4,
                       threshold=10000,
                       linewidth=150,
                       suppress=True)

    M = 500

    rho = 0.5
    mu = [0,0]
    R = [[1,rho],[rho,1]]
    x = np.random.multivariate_normal(mu,R,M)
    rdcVal = rdc(x[:,0],x[:,1])
    print(rdcVal)

    # test perfect non-monotonic functional dependence
    x = np.random.uniform(low=-1.0,high=1.0,size=(M,1))
    y = np.power(x,2)
    rdcVal = rdc(x,y)
    print(rdcVal)