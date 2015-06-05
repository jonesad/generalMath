# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 18:20:33 2015

general matrix manipulations

@author: jonesad
"""
#List columns if their norm is less than tol
def getZeroCols(a, tol=10**-5):
   rmlist=[]
   import numpy
   for nIdx in range(a.shape[1]) :
     if numpy.linalg.norm(a[:,nIdx])<tol:
       rmlist.append(nIdx)
   return rmlist
#List the set of row labels of a 2d array that are repeats of previous row values to a given tolerance
def getRepeatSlices(a,axis=0 ,tol=10**-5):
  import numpy as np
  rmList=[]
  slc=[slice(None)]*len(a.shape)
  for nRIdx in range(a.shape[axis]):
    temp1=list(slc)
    temp1[axis]=nRIdx
    for nTest in range(nRIdx + 1,a.shape[axis]):
      temp2=list(slc)
      temp2[axis]=nTest
#      print nRIdx, nTest
      if np.allclose(a[temp1],a[temp2],atol=tol):
#        print temp2==temp1
        rmList.append(nTest)
  rmList=sorted(list(set(rmList)))
  return rmList  

#remove slices of a matrix a that are listed in rm list along axis 
def rmSlice(rmlist, a, nAxis):     
   import numpy
   inlist=[aa for aa in range(a.shape[nAxis]) if aa not in rmlist] 
#   print 'inlist', len(inlist)   
   npaCondition=numpy.zeros(a.shape[nAxis], dtype=bool)
   nInIdx=0
   for nIdx in range(npaCondition.size):
     if nInIdx>=len(inlist):
       break
     elif inlist[nInIdx]==nIdx:
       npaCondition[nIdx]=1
       nInIdx+=1
   a=numpy.compress(npaCondition, a,nAxis)
   return  a
#     end of remove the zeros
 
 #remove the elements of a list indexed by the elements in another list
def rmElem(rmList,datList):
  inlist=[aa for aa in range(len(datList)) if aa not in rmList]
  b=[]
  for nIdx in inlist:
    if b!=[]:
      b.append(b,datList[nIdx])
    else:
      b=[datList[nIdx]]
  datList=b
  return datList
  
# construct the wighted least square solution for correlated data
#return  the solution vector and the solution vector errors.
     
def weightedlsq(npaCoef, npaWeights, npaData, npaCorr=[]):
  import numpy as np
  xtw=np.dot(np.transpose(npaCoef),npaWeights)
  print np.dot(xtw,npaCoef)
  xtwx_inv=np.linalg.inv(np.dot(xtw,npaCoef))
  xtwx_invxtw=np.dot(xtwx_inv,xtw)
  npaSolution=np.dot(xtwx_invxtw,npaData)
  if npaCorr==[]:
    npaErrors=np.diag(xtwx_inv)
  else:
    temp1=np.dot(xtwx_invxtw,npaCorr)
    temp1=np.dot(temp1,np.transpose(npaWeights))
    temp1=np.dot(temp1,npaCoef)
    temp2=np.dot(np.transpose(npaCoef),np.transpose(npaWeights))
    temp2=np.dot(temp2,npaCoef)
    npaErrors=np.dot(temp1, np.linalg.inv(temp2))
  return npaSolution, npaErrors
