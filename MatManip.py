# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 18:20:33 2015

general matrix manipulations

@author: jonesad
"""


# List columns if their norm is less than tol
def getZeroCols(a, tol=10**-5):
    rmlist = []
    import numpy
    for nIdx in range(a.shape[1]):
        if numpy.linalg.norm(a[:, nIdx]) < tol:
            rmlist.append(nIdx)
            #       print nIdx
    return rmlist


# List the set of row labels of a 2d array that are repeats of previous row
# values to a given tolerance
def getRepeatSlices(a, axis=0, tol=10**-5):
    import numpy as np
    rmList = []
    slc = [slice(None)] * len(a.shape)
    for nRIdx in range(a.shape[axis]):
        temp1 = list(slc)
        temp1[axis] = nRIdx
        for nTest in range(nRIdx + 1, a.shape[axis]):
            temp2 = list(slc)
            temp2[axis] = nTest
            if np.allclose(a[temp1], a[temp2], atol=tol):
                rmList.append(nTest)
    rmList = sorted(list(set(rmList)))
    return rmList


# remove slices of a matrix a that are listed in rm list along axis
def rmSlice(rmlist, a, nAxis):
    import numpy
    inlist = [aa for aa in range(a.shape[nAxis]) if aa not in rmlist]
    npaCondition = numpy.zeros(a.shape[nAxis], dtype=bool)
    nInIdx = 0
    for nIdx in range(npaCondition.size):
        if nInIdx >= len(inlist):
            break
        elif inlist[nInIdx] == nIdx:
            npaCondition[nIdx] = 1
            nInIdx += 1
    a = numpy.compress(npaCondition, a, nAxis)
    return a


# remove the elements of a list indexed by the elements in another list
def rmElem(rmList, datList):
    inlist = [aa for aa in range(len(datList)) if aa not in rmList]
    b = []
    for nIdx in inlist:
        if b != []:
            b.append(b, datList[nIdx])
        else:
            b = [datList[nIdx]]
    datList = b
    return datList


# construct the wighted least square solution for correlated data
# return  the solution vector and the solution vector errors.
def weightedlsq(npaCoef, npaWeights, npaData, npaCorr=[]):
    import numpy as np
    xtw = np.dot(np.transpose(npaCoef), npaWeights)
    xtwx_inv = np.linalg.inv(np.dot(xtw, npaCoef))
    xtwx_invxtw = np.dot(xtwx_inv, xtw)
    npaSolution = np.dot(xtwx_invxtw, npaData)
    if npaCorr == []:
        npaErrors = np.diag(xtwx_inv)
    else:
        temp1 = np.dot(xtwx_invxtw, npaCorr)
        temp1 = np.dot(temp1, np.transpose(npaWeights))
        temp1 = np.dot(temp1, npaCoef)
        temp2 = np.dot(np.transpose(npaCoef), np.transpose(npaWeights))
        temp2 = np.dot(temp2, npaCoef)
        npaErrors = np.dot(temp1, np.linalg.inv(temp2))
    return npaSolution, npaErrors


# take 2 lists/npas of labels and corresponding 2d npas of values and make a
# consitent list/npa of labels and a 2d npa that has all the elements of the
# original set with zeros in the places where the labels are inserted in the
# list.
def combinedLabeledColumns(npaLab1, npaVal1, npaLab2, npaVal2):
    import numpy as np
    if type(npaLab1) != type(npaLab2):
        try:
            npaLab1 = np.array(npaLab1)
            npaLab2 = np.array(npaLab2)
        except:
            print 'Error: Labels must be able to be cast to same type.'
            print 'type(Label1): ', type(npaLab1)
            print 'type(Label2): ', type(npaLab2)
            return None
    bComp1 = isinstance(npaLab1, type(np.array([])))
    bComp2 = isinstance(npaLab2, type(np.array([])))
    if bComp1 and bComp2:
        npaLab1 = list(npaLab1)
        npaLab2 = list(npaLab2)
    elif not(isinstance(npaLab1, list) and isinstance(npaLab2, list)):
        print ('Error: Labels must be of type: ', type(np.array([])),
               ' or type:', list)
        return None
    import MatManip
    lLabUnion = MatManip.lnpaUnion(npaLab1, npaLab2)
    MatManip.sortnpaList(lLabUnion)
    npaValNew = np.zeros([npaVal1.shape[0] + npaVal2.shape[0],
                          len(lLabUnion)])
    for nIdx, lab in enumerate(lLabUnion):
        for nIdx1, lab1 in enumerate(npaLab1):
            if np.all(lab == lab1):
                npaValNew[:npaVal1.shape[0], nIdx] = npaVal1[:, nIdx1]
        for nIdx2, lab2 in enumerate(npaLab2):
            if np.all(lab == lab2):
                npaValNew[(npaVal1.shape[0]):, nIdx] = npaVal2[:, nIdx2]
    return lLabUnion, npaValNew


# Test if the intersection of 2 lists occurs in the same order in each list.
def intersectInOrder(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = list(set1 & set2)
    lOrder1 = []
    lOrder2 = []
    for elem1 in intersection:
        for nIdx, elem2 in enumerate(list1):
            if elem1 == elem2:
                lOrder1.append(nIdx)
        for nIdx, elem2 in enumerate(list2):
            if elem1 == elem2:
                lOrder2.append(nIdx)
    import numpy as np
    return np.all(lOrder1 == lOrder2)


# take a list of numpy arrays and sorts them in oprder of their flattened
# elements. Returns 
def sortnpaList(List1):
    import numpy as np
    if len(List1) == 0 or len(List1) == 1:
        return List1
    for nIdx, elem in enumerate(List1):
        try:
            List1[nIdx] = np.array(elem)
        except:
             ''               
        if not (isinstance(elem, type(np.array([]))) or isinstance(elem, list)):
            print ('Error: element ', nIdx, ' of input is of type ',
                   type(elem), 'only numpy arrays and lists are valid.')
            return None
    for npaElem in List1:
        if npaElem.size != List1[0].size:
            print 'Error: All arrays must be of the same size.'            
            return None
    for nIdx in range(List1[0].size):
        for nIdx1 in range(len(List1)):
            temp = range(len(List1)-nIdx1)
            lnIdxs = []
            for num in temp:
                lnIdxs.append(num + nIdx1)
            for nIdx2 in lnIdxs:
                bComp1 = List1[nIdx1][nIdx] > List1[nIdx2][nIdx]
                bComp2 = np.all(List1[nIdx1][0:nIdx] == List1[nIdx2][0:nIdx])
                if bComp1 and bComp2:
                    temp = List1[nIdx2]
                    List1[nIdx2] = List1[nIdx1]
                    List1[nIdx1] = temp
    return None


# take a list of npas and return the list minus the repitions
def lnpaSet(lnpa):
    lnpaOut = []
    import numpy as np
    for elem1 in lnpa:
        bIsIt = False
        for elem2 in lnpaOut:
            bIsIt = np.all(elem1 == elem2)
            if bIsIt:
                break
        if not bIsIt:
            lnpaOut.append(elem1)
    return lnpaOut


def lnpaIntersect(lnpa1, lnpa2):
    lnpaOut = []
    import numpy as np
    for elem1 in lnpa1:
        for elem2 in lnpa2:
            if np.all(elem1 == elem2):
                bIsIt = False
                for elem3 in lnpaOut:
                    bIsIt = np.all(elem1 == elem3)
                    if bIsIt:
                        break
                if not bIsIt:
                    lnpaOut.append(elem1)
    return lnpaOut

# takes union of lnpa1 and lnpa2 only pass sets of distinguishable elements
def lnpaUnion(lnpa1, lnpa2):
    import numpy as np
    lnpaOut = [elem for elem in lnpa1]
    for elem2 in lnpa2:
        bIsIt = False
        for elem1 in lnpa1:
            bIsIt = np.all(elem1 == elem2)
            if bIsIt:
                break
        if not bIsIt:
            lnpaOut.append(elem2)
    return lnpaOut
