#!/usr/bin/env python
# encoding: utf-8
"""
fakeFid.py

Created by Daniel O'Donovan on 2010-09-27.
Copyright (c) 2010 Daniel O'Donovan. All rights reserved.
"""

import sys, os

import numpy as np
import matplotlib.pyplot as mp

import processing

class FID:
  """ Fake Free Induction Decay generation class """

  def __init__(self, size, time, dtype=np.complex):

    self.pi2  = np.pi * 2.
    self.dt   = time / float(size)
    self.size = size
    self.time = time
    self.td   = { 'time':np.arange( 0., time, self.dt ), \
                  'td':np.empty( size, dtype=np.complex ) }
    self.fd   = { 'freq':-1.*np.arange( -1./(time*2.), 1./(time*2.), (2./(time*2.))/size ), \
                  'fd':np.empty( size, dtype=np.complex ) }

  def addSignal( self, height, frequency, T2 ):
    """ Add an FID signal at given height, frequency and T2, phase 0 0 """

    if frequency < ( -1./(self.time*2.) ) or frequency > (1./(self.time*2.)):
      print '&&& Warning signal at %f is out of spectral range, signal will be aliased.' % frequency
      print '&&& Frequency must lie in %f < f < %f' % ( -1./(self.time*2.), 1./(self.time*2.))

    ( h, f ) = ( height, ((self.pi2 * self.time) * frequency) / self.dt )

    t               = self.td['time']
    self.td['td']   += h * np.exp( + 1.j * t * f ) * np.exp( - t / T2 )

  def fft( self ):

    self.fd['fd']    = processing.nmr_fft( self.td['td'] )

if __name__ == '__main__':

  fid = FID( size=2048, time=2.2 )
  
  fid.addSignal( 10.,  0.00, .5 )
  
  fid.addSignal( 10.,  0.10, .5 )
  fid.addSignal( 10., -0.10, .5 )
  
  fid.addSignal( 20.,  0.18, .5 )
  fid.addSignal( 30., -0.18, .5 )
  
  fid.addSignal( 20., -0.20, .5 )

  # fid = FID( size=2048, time=2048. )
  # 
  # fid.addSignal( 10.,  0.0, 1.0E3 )
  # 
  # fid.addSignal( 10.,  0.0001, 1.0E3 )
  # fid.addSignal( 10., -0.0001, 1.0E3 )
  # 
  # fid.addSignal( 20.,  0.00018, 1.0E3 )
  # fid.addSignal( 30., -0.00018, 1.0E3 )
  # 
  # fid.addSignal( 20., -0.0002, 1.0E3 )

  fid.fft()

  mp.subplot( 211 )
  mp.plot( fid.td['time'], fid.td['td'].real )
  mp.plot( fid.td['time'], fid.td['td'].imag )
  
  mp.subplot( 212 )
  mp.plot( fid.fd['freq'], fid.fd['fd'].real )
  mp.plot( fid.fd['freq'], fid.fd['fd'].imag )







