# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.30

model.py describes the models to infer parameters and to select optimal mechanistic model. For actual mechanistic models please refer to modelBase.ph
'''

import warnings
import time
import numpy as np
from sys import stdout
from enum import Enum
from scipy.optimize import minimize, differential_evolution
from functools import reduce
from biologiclib import ioUtils, modelBase, plotUtils
from biologiclib.modelBase import ModelType, ModelSpec, ModelSet
from sympy import symbols, lambdify, diff, pretty
import multiprocessing as mp

# ModelSolver defines the possible procedures to estimate the parameters of a single mechanistic model
ModelSolver = Enum("ModelSolver", ("Nelder_Mead", "N_M", "BFGS", "COBYLA", "SLSQP"))

# ModelCriterion defines the standard to evaluate a given model{expression and the parameters}, basically Infromation Criteria
ModelCriterion = Enum("ModelCriterion", ("AIC"))

class Solution:
    '''
    To Store the metadata of a modol solution
    '''

    def __init__(self, modelType, modelSpecs, expression, thetaKey, thetaVal, residue, IC):
        self.modelType = modelType
        self.modelSpecs = modelSpecs
        self.expression = expression
        self.thetaKey = thetaKey
        self.thetaVal = thetaVal
        self.residue = residue
        self.IC = IC

# thetaList is necessary for buidling sse func. However as thetaList should be considered as a part of the model, no additional info. is needed.
# inducer E R^(k*N), reporter E R^N, reporterStd E R^N
def genSSE(inducer, reporter, reporterStd, model, thetaList):
    def SSE(X):
        # Check # parameter
        if len(X) != len(thetaList):
            raise Exception('Incorrect number of parameters. Expected %d but %d parameters are given'%(len(thetaList), len(X)))

        # Check dimension of inducer & reporter
        if len(inducer) != len(reporter):
            raise Exception('Dimension of inducer and reporter does not match. Inducer(%d) but reporter(%d). Abort.'%(len(inducer), len(reporter)))

        # Composing theta **Note that the order of X is important**
        theta = {key: value for key, value in zip(thetaList, X)}

        # mu is the "real" reporter value, assuming data is disturbed
        mu = model(np.array(inducer), theta)
        sse = sum([(m - x) ** 2 for m, x in zip(mu, reporter)])

        return sse

    return SSE

def genJac(inducer, reporter, reporterStd, expression, thetaList):
    '''
    Generate the Jacobian of the *SSE*, mainly for optimization's purpose.
    Consider the function: P_st = expression(inducer; theta), where theta E R^k are the parameters. genJac renders a vector of partial derivitive to each parameter.
    '''

    # Check parameter dimensions
    try:
        if len(inducer) != len(reporter) or len(inducer) != len(reporterStd):
            raise Exception("Usage: genJac(inducer, reporter, reporterStd, expression, thetaList), where inducer, reporter, and reporterStd should be of the same dimension. However, inducer(%d), reporter(%d), and reporterStd(%d)"%(len(inducer), len(reporter), len(reporterStd)))
    except TypeError:
        raise Exception("Usage: genJac(inducer, reporter, reporterStd, expression, thetaList), where inducer, reporter, and reporterStd should be iterable")

    # Claim symbols
    thetaSymbols = symbols(thetaList)
    # Original function lambdify
    PstFunc = lambdify((symbols('A'), thetaSymbols), expression, "numpy")
    # Symbolic derivitives
    jacExp = [diff(expression, x) for x in thetaSymbols]
    # Derivitives lambdify
    jacFunc = [lambdify((symbols('A'), thetaSymbols), exp, "numpy") for exp in jacExp]

    # Construct the function
    def jacobian(theta):
        jac = []
        for i in range(len(thetaList)):
            # derivitive of sse is a linear combinition
            try:
                len(inducer[0])    # empty imput is not allowed
                dSSE = 2 * sum(
                    [(PstFunc(A[0], theta) - x) * jacFunc[i](A[0], theta)\
                        for A, x in zip(inducer, reporter)]
                )
            except TypeError:
                dSSE = 2 * sum(
                    [(PstFunc(A, theta) - x) * jacFunc[i](A, theta)\
                        for A, x in zip(inducer, reporter)]
                )
            jac.append(dSSE)
        return np.array(jac)    # Some scipy optimization funcs, e.g. BFGS, insist on using np.array jacobian

    # Return
    return jacobian

def estimatePara(sse, X0, jacobian = None, constraints = None, bounds = None, method = ModelSolver.Nelder_Mead):
    # Check input
    if type(method) != ModelSolver:
        raise Exception('Usage: paraEstimator(sse, X0, method), where method should be ModelSolver type Enum.')

    # Tentatively add bounds to ensure non-negative solutions
    bounds = [(0, x * 10) for x in X0]

    # Always first use Nelder_Mead to estimate initial theta
    NMres = minimize(sse, X0, method = "Nelder-Mead", options = {
        "maxiter": 20,
        "disp": False
    })
    newX0 = NMres.x

    # Alternative solvers
    if method == ModelSolver.Nelder_Mead or method == ModelSolver.N_M:
        res = NMres
    elif method == ModelSolver.COBYLA:
        res = minimize(sse, newX0, constraints = constraints, method = "COBYLA", options = {
            "tol": 1E-4,
            "maxiter": 50,
            "disp": False
        })
    elif method == ModelSolver.BFGS:
        res = minimize(sse, newX0, jac = jacobian, method = "BFGS", options = {
            "maxiter": 50,
            "gtol": 1E-4,
            "disp": False
        })
    elif method == ModelSolver.SLSQP:
        try:
            res = minimize(sse, newX0, constraints = constraints, jac = jacobian, method = "SLSQP", options = {
                "maxiter": 50,
                "ftol": 1E-2,    # to avoid "ComplexInfinity" error
                "disp": False
            })
        except ValueError:    # hack to continue
            warnings.warn("Algorithm overflow. This model renders no feasible solutions")
            return X0, 1E10

    # Return
    return res.x, res.fun

def infoCriteria(sse, theta, reporter, method = ModelCriterion.AIC):
    if method == ModelCriterion.AIC:
        k, n = len(theta), len(reporter)
        AIC = 2 * k + n * np.log(sse / n)
        # AIC correction when small sample size
        #AICc = AIC + 2 * k * (k + 1) / (n - k - 1)
        return AIC
    else:
        return sse

gInducer, gReporter, gReporterStd = [], [], []
gModelSolver, gModelCrterion = None, None

def __fitModel(modelKeys):
    '''
    Paralleled version of model fitting
    '''

    modelType, modelSpecs = modelKeys

    # Generate the model
    exp, model, constraints, thetaList = modelBase.genModel(modelType, modelSpecs, plain_print=False)

    # Get SSE
    sse = genSSE(gInducer, gReporter, gReporterStd, model, thetaList)
    # Also get jacobian of SSE
    jacobian = genJac(gInducer, gReporter, gReporterStd, exp, thetaList)

    # Is this repression or activation model?
    repression = False
    if ModelSpec.Repression in modelSpecs and ModelSpec.Inducer_Repression not in modelSpecs:
        repression = True
    elif ModelSpec.Repression not in modelSpecs and ModelSpec.Inducer_Repression in modelSpecs:
        repression = True

    # Parameterization
    startTime = time.time()
    inferredTheta, residue = estimatePara(
            sse, modelBase.defaultPara(thetaList, gInducer, gReporter, repression = repression),
            jacobian, constraints, method = gModelSolver
    )
    duration = time.time() - startTime

    # Calculate Infomation Criteria
    IC = infoCriteria(residue, inferredTheta, gReporter, gModelCriterion)

    # Compose model meta
    meta = Solution(
            modelType, modelSpecs,
            pretty(exp, use_unicode=False),
            thetaList, inferredTheta,
            residue, IC)

    return meta, duration

def selectModel(inducer, reporter, reporterStd,
        inducerTags, replicateTags, reporterTag,
        modelSolver = ModelSolver.N_M, modelSet = ModelSet.Simple_Inducible_Promoter, modelCriterion = ModelCriterion.AIC,
        inducerUnits = "M", reporterUnits = "M",
        figPath = None, quiet=True, plot = True):
    '''
    Select the best model interpreting given inducer/reporter characterization data.
    Notice that fitting of alternative models is paralled. Thus, selectModel() cannot be run multiprocessed!
    '''

    global gInducer, gReporter, gReporterStd, gModelSolver, gModelCriterion
    gInducer, gReporter, gReporterStd = inducer, reporter, reporterStd
    gModelSolver, gModelCriterion = modelSolver, modelCriterion

    # Generate model candicates (alternative mechanistic hypotheses)
    modelSet = modelBase.genModelSet(modelSet)

    # parallel model fitting
    metas = []
    minIC = float("inf")
    count = 0    # number of model being solved

    cpuCount = int(mp.cpu_count())    # Num of cpu cores
    pool = mp.Pool(cpuCount)    # A pool of processes
    for meta, duration in pool.imap_unordered(__fitModel, modelSet):
        count += 1
        metas.append(meta)

        # Print model info
        if not quiet:
            print("========================================")
            print("# %d Model"%count)
            printModel(meta)
            print("Time elapsed:%.2f\n"%duration)

        # Best model
        if meta.IC < minIC:
            minIC = meta.IC
            bestModelMeta = meta
    pool.close()
    pool.join()

    # Plotting
    thetaDict = {key: val for key, val in zip(bestModelMeta.thetaKey, bestModelMeta.thetaVal)}
    _, bestModel, _, _ = modelBase.genModel(bestModelMeta.modelType, bestModelMeta.modelSpecs)
    if plot:
        plotUtils.plotModel2D(inducer, reporter, reporterStd, thetaDict, bestModel,
                inducerTags[0], reporterTag, inducerUnits, reporterUnits, figPath)

    return bestModelMeta, metas

def printModel(meta):
    print('=============================================')
    print('%s\nType = %s'%(meta.expression, meta.modelType.name))
    print('Specs = ' + ', '.join([spec.name for spec in meta.modelSpecs]))
    print('Parameters:\n' +\
            ', '.join([key + ' = ' + str(val) for key, val in zip(meta.thetaKey, meta.thetaVal)]))
    print('Residue = %.4f\nIC = %.4f'%(meta.residue, meta.IC))
    stdout.flush()
