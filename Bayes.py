# -*- coding: utf-8 -*-
"""
Created on Sun Sep 06 19:33:00 2015

bayesian code

@author: Adam
"""


def selectModel(npafIndepVar, npafData, funErrorPDF, lfunPriorParamPDF,
                lllfPriorParamBounds, lfunModels, lfPriorOdds=[]):
    '''
        Select between two models using a naive quadrature method for
        integration of probability distributions. Inputs are the data to model,
        the error pdf for the data, the parameter probability density, a list
        of 2 element lists of the [lower bound, upper bound] of the arguments
        of parameter distribution, the list of models to select from and a list
        of lists of 2 element lists of the model parameter bounds. For discrete
        parameters, select between functions defined with the different values
        of the deiscreet parameters. The outputs are the index of the best
        model and a list of the odds ratios for each model relative to the best
        model.
    '''
    from scipy.integrate import nquad
    import numpy as np
    npaOR = []
    if lfPriorOdds == []:
        lfPriorOdds = np.ones(len(lfunModels))
    nIdx = 0
    for funParamPDF, llfParamBounds, funModel in zip(lfunPriorParamPDF,
                                                     lllfPriorParamBounds,
                                                     lfunModels):
        funIntegrand = lambda *x: (funParamPDF(x) *
                                   funErrorPDF(npafData -
                                   funModel(npafIndepVar, x)))
        npaOR.append(nquad(funIntegrand, llfParamBounds)*lfPriorOdds[nIdx])
        nIdx += 1
    npaOR = np.array(npaOR)
    npaOR = npaOR/np.max(npaOR)
    nBestModel = np.argmin(npaOR)
    return nBestModel, npaOR


def normGauss(fX, fMu, fSigma):
    import numpy as np
    num = np.exp(-(fX-fMu)**2(2.*fSigma**2))
    denom = fSigma*np.sqrt(2.*np.pi)
    return num/denom


def uniformDist(lfX, fLow, fHigh):
    if type(lfX) == float or type(lfX) == int:
        if lfX > fLow and lfX < fHigh:
            return 1. / (fHigh - fLow)
        else:
            return 0.
    else:
        import numpy as np
        lfX = np.array(lfX)
        if np.all(lfX > fLow) and np.all(lfX < fHigh):
            return 1./(fHigh-fLow)**lfX.size
        else:
            return 0

funCon = lambda (x, a): a
funLin = lambda (x, a, b): a*x + b
funQuad = lambda (x, a, b, c): a*x**2 + b*x + c
lfunModels = [funCon, funLin, funQuad, normGauss]

import numpy as np
npafIndepVar = np.linspace(-10, 10, 100)
npafData = np.random.normal(0., 1., npafIndepVar.shape)
funErrorPDF = lambda x: normGauss(x, 0., 1.)
lfunPriorParamPDF = [lambda x: uniformDist(x, -10, 10)]*len(lfunModels)
lllfPriorParamBounds = [[[-10, 10]], [[-10, 10], [-10, 10]], [[-10, 10],
                         [-10, 10], [-10, 10]], [[-10, 10], [-10, 10]]]
selectModel(npafIndepVar, npafData, funErrorPDF, lfunPriorParamPDF,
            lllfPriorParamBounds, lfunModels, lfPriorOdds=[])