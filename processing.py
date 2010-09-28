#!/usr/bin/env python
# encoding: utf-8
"""
fakeFid.py

Created by Daniel O'Donovan on 2010-09-27.
Copyright (c) 2010 Daniel O'Donovan. All rights reserved.
"""

import numpy

# A first stab at some useful NMR processing functions to be kept in one place

def real2complex( d ):
    return numpy.array( [complex( d[i], d[i+1] ) for i in range( 0, len(d), 2 ) ] )

def nmr_fft( d, phase=None, ph00=False ):
    
    if ph00:
        d[0] = (d[-1] + d[0]) / 2

    d = numpy.fft.fft( d )

    if phase is not None:
        for i in xrange( len(d) ):
            d[i] *= numpy.exp( (1.j) * (phase[0] + float(i) * phase[1] / float(len(d)) ) )

    d = numpy.concatenate( (numpy.split( d, 2 )[1], numpy.split( d, 2 )[0]) )

    d = d[::-1]

    return d

def zerofill( d, prop ):
    z = numpy.zeros( len( d ), dtype=complex )
    nz = int( prop * len( d ) )  
    for i in xrange( len(d) - nz ): z[i] = d[i]
    return z
    
def cumulativeArray(array):

  ndim = len(array)
  cumul = ndim * [0]
  n = 1
  for i in xrange(len(array)):
    cumul[i] = n
    n = n * array[i]

  return (n, cumul)

def arrayOfIndex(index, cumul):

  ndim = len(cumul)
  array = ndim * [0]
  for i in xrange(ndim-1, -1, -1):
    c = cumul[i]
    array[i], index = divmod(index, c)

  return tuple(array)

def indexOfArray(array, cumul):

  index = 0
  for i in xrange(len(array)):
    index += array[i] * cumul[i]

  return index


def nmr_fft_Nd( data, dataShape, phase=None, ph00=False ):
  """ dataShape is actual number of points (which are hypercomplex in > 1 D) """

  def getTransposeDims( dim, ndim ):
    """ [0,1,2] -> [0,1,2] : [0,1,2] -> [1,0,2] : etc. """
    trans       = range( ndim ) 
    trans[0]    = dim
    trans[dim]  = 0
    return trans

  def complexToRealFID( fid ):
    realFid = numpy.array( [ [f.real, f.imag] for f in fid] ).ravel()
    return realFid

  def realToComplexFID( realFid ):
    fid = []
    for i in xrange( 0, len( realFid ), 2 ):
      fid.append( complex( realFid[i], realFid[i+1] ) )
    fid = numpy.array( fid )
    return fid

  data = numpy.array( data )
  data = data.reshape( dataShape, order='F' )

  dataShape = numpy.array( dataShape )

  # do NMR fft in each dim
  for dim in xrange( len( dataShape ) ):

    trans = getTransposeDims( dim, len( dataShape ) )
    data = data.transpose( trans )

    newShape = numpy.array( data.shape )
    dimRem = newShape.prod() / newShape[0]

    # reshape 
    data = data.reshape( [newShape[0], dimRem], order='F' )

    for i in xrange( dimRem ):
      
      fid = realToComplexFID( data[:,i] )
      fid = nmr_fft( fid, phase=phase, ph00=ph00 )
      data[:,i] = complexToRealFID( fid )

    # return to post-transpose shape
    data = data.reshape( newShape, order='F' )

    # untranspose the data
    data = data.transpose( trans )

  return data

# NOTE THIS REDUCE ISN'T EXPECTED TO WORK YET!!!
def reduceData( data, dimSize, dims=None):
  """ Remove imaginary parts of dims=dimList from data """

  ndims = len( dimSize )

  if not dims:
    dims = range( ndims )

  for dim in dims:

    trans = getTransposeDims( dim, len( dataShape ) )
    data = data.transpose( trans )

    newShape = numpy.array( data.shape )
    dimRem = newShape.prod() / newShape[0]

    # reshape 
    data = data.reshape( [newShape[0], dimRem], order='F' )

    for i in xrange( dimRem ):

      data[:,i] = complexToRealFID( fid )

    # return to post-transpose shape
    data = data.reshape( newShape, order='F' )

    # untranspose the data
    data = data.transpose( trans )

  return data






