#!/usr/bin/env python

import sys
import numpy
import math
import random

def bootstrapping(data, nsample=200):
	ndata = data.size
	m_re = numpy.zeros(nsample, float)
	for i in xrange(0, nsample):
		for n in xrange(0, ndata):
			rnd_ndx=random.randint(0,ndata-1)
			m_re[i] += data[rnd_ndx]
	m_re /= ndata
	return m_re.var()

def statisticalInefficiencyMultiple(A_kn, fast=False, return_correlation_function=False):
  """
  Estimate the statistical inefficiency from multiple stationary timeseries (of potentially differing lengths).

  REQUIRED ARGUMENTS  
    A_kn (Python list of numpy arrays) - A_kn[k] is the kth timeseries, and A_kn[k][n] is nth value of timeseries k.  Length is deduced from arrays.

  OPTIONAL ARGUMENTS  
    fast can be set to True to give a less accurate but very quick estimate (default False)
    return_correlation_function - if True, will also return estimates of normalized fluctuation correlation function that were computed (default: False)

  RETURNS
    g is the statistical inefficiency (equal to 1 + 2 tau, where tau is the integrated autocorrelation time).
    Ct (list of tuples) - Ct[n] = (t, C) with time t and normalized correlation function estimate C is returned as well if return_correlation_function is set to True
    
  NOTES 
    The autocorrelation of the timeseries is used to compute the statistical inefficiency.
    The normalized fluctuation autocorrelation function is computed by averaging the unnormalized raw correlation functions.
    The fast method described in Ref [1] is used to compute g.

  REFERENCES  
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
    histogram analysis method for the analysis of simulated and parallel tempering simulations.
    JCTC 3(1):26-41, 2007.

  EXAMPLES

  Estimate statistical efficiency from multiple timeseries of different lengths.

  >>> import testsystems
  >>> N_k = [1000, 2000, 3000, 4000, 5000]
  >>> tau = 5.0 # exponential relaxation time
  >>> A_kn = [ testsystems.generateCorrelatedTimeseries(N=N, tau=tau) for N in N_k ]
  >>> g = statisticalInefficiencyMultiple(A_kn)

  Also return the values of the normalized fluctuation autocorrelation function that were computed.

  >>> [g, Ct] = statisticalInefficiencyMultiple(A_kn, return_correlation_function=True)

  """

  # Convert A_kn into a list of arrays if it is not in this form already.
  if (type(A_kn) == numpy.ndarray):
    A_kn_list = list()
    if A_kn.ndim == 1:      
      A_kn_list.append(A_kn.copy())
    else:
      [K,N] = A_kn.shape
      for k in range(K):
        A_kn_list.append(A_kn[k,:].copy())
    A_kn = A_kn_list
      
  # Determine number of timeseries.
  K = len(A_kn)

  # Get the length of each timeseries.
  N_k = numpy.zeros([K], numpy.int32)
  for k in range(K):
    N_k[k] = A_kn[k].size

  # Compute average timeseries length.
  Navg = numpy.array(N_k, numpy.float64).mean()

  # Determine total number of samples.
  N = sum(N_k)

  # Initialize statistical inefficiency estimate with uncorrelated value.
  g = 1.0
    
  # Compute sample mean.
  mu = 0.0
  for k in range(K):
    mu += A_kn[k].sum()
  mu /= float(N)
  
  # Construct and store fluctuation timeseries.
  dA_kn = list()
  for k in range(K):
    dA_n = A_kn[k] - mu
    dA_kn.append(dA_n.copy())

  # Compute sample variance from mean of squared fluctuations, to ensure that C(0) = 1.
  sigma2 = 0.0
  for k in range(K):
    sigma2 += (dA_kn[k]**2).sum()
  sigma2 /= float(N)
  
  # Initialize statistical inefficiency estimate with uncorrelated value.
  g = 1.0

  # Initialize storage for correlation function.
  Ct = list() # Ct[n] is a tuple (t, C) of the time lag t and estimate of normalized fluctuation correlation function C
    
  # Accumulate the integrated correlation time by computing the normalized correlation time at
  # increasing values of t.  Stop accumulating if the correlation function goes negative, since
  # this is unlikely to occur unless the correlation function has decayed to the point where it
  # is dominated by noise and indistinguishable from zero.
  t = 1
  increment = 1  
  while (t < N_k.max()-1):
    # compute unnormalized correlation function
    numerator = 0.0
    denominator = 0.0
    for k in range(K):
      if (t >= N_k[k]): continue # skip trajectory if lag time t is greater than its length
      dA_n = dA_kn[k] # retrieve trajectory
      x = dA_n[0:(N_k[k]-t)] * dA_n[t:N_k[k]]
      numerator += x.sum() # accumulate contribution from trajectory k
      denominator += float(x.size) # count how many overlapping time segments we've included

    C = numerator / denominator

    # compute normalized fluctuation correlation function at time t
    C = C / sigma2
    #print "C[%5d] = %16f (%16f / %16f)" % (t, C, numerator, denominator)    

    # Store estimate of correlation function.
    Ct.append( (t,C) )

    # Terminate if the correlation function has crossed zero.
    # Note that we've added a hack (t > 10) condition to avoid terminating too early in correlation functions that have a strong negative peak at 
    if (C <= 0.0) and (t > 10):
      break
  
    # Accumulate contribution to the statistical inefficiency.
    g += 2.0 * C * (1.0 - float(t)/Navg) * float(increment)

    # Increment t and the amount by which we increment t.
    t += increment

    # Increase the interval if "fast mode" is on.
    if fast: increment += 1

  # g must be at least unity
  if (g < 1.0): g = 1.0

  # Return statistical inefficency and correlation function estimate, if requested.
  if return_correlation_function:
    return (g, Ct)

  # Return the computed statistical inefficiency.
  return g

sys.stdout.softspace=False
if __name__ == "__main__" :
	lines = open(sys.argv[1])
	nskip = 1
	if (len(sys.argv)>2):
		nskip = int(sys.argv[2])
	nbootstrap = 200
	if (len(sys.argv)>3):
		nbootstrap = int(sys.argv[3])
	data = []
	
	nline = 0
	for line in lines:
		numbers = line.split()
		sample = []
		nline = nline+1
		if (nline % nskip <> 0): continue
		for number in numbers:
			sample.append( float(number) )
		data.append(sample)

	nrow = len(data)
	ncol = len(data[0])

	for ndx in range(ncol):
		series = numpy.zeros(nrow, float)
		sum2 = 0
		delta = 0
		for i in range(nrow):
			series[i] = data[i][ndx]
			if i > 0:
				delta = series[i] - series[i-1]
			sum2 += delta*delta
		if sum2 > 0.0001:	
			g = statisticalInefficiencyMultiple(series)
			m = series.mean()
			v = bootstrapping(series, nbootstrap)
			print "col %2d skip %3d mean %8.6f var %8.6f g %8.6f" % (ndx, nskip, m, v*g, g)

