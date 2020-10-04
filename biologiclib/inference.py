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
#from biologiclib.modelBase import jacBase
from sympy import symbols, lambdify, diff, pretty
import multiprocessing as mp
import pymc3 as pm
import logging
logger = logging.getLogger("pymc3")
logger.setLevel(logging.ERROR)    # hack to hide NUTS initiation message

# ModelSolver defines the possible procedures to estimate the parameters of a single mechanistic model
ModelSolver = Enum("ModelSolver", ("Nelder_Mead", "N_M", "BFGS", "COBYLA", "SLSQP", "MCMC"))

# ModelCriterion defines the standard to evaluate a given model{expression and the parameters}, basically Infromation Criteria
ModelCriterion = Enum("ModelCriterion", ("AIC"))

class Solution:
    '''
    To Store the metadata of a modol solution
    '''

    def __init__(self, modelType, modelSpecs, expression, thetaKey, thetaVal, residue, IC, thetaStd=None):
        self.modelType = modelType
        self.modelSpecs = modelSpecs
        self.expression = expression
        self.thetaKey = thetaKey
        self.thetaVal = thetaVal
        self.thetaStd = thetaStd
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
        if len(inducer) != len(reporter):
            raise Exception("Usage: genJac(inducer, reporter, reporterStd, expression, thetaList), where inducer and reporter should be of the same dimension. However, inducer(%d) and reporter(%d)"%(len(inducer), len(reporter)))
    except TypeError:
        raise Exception("Usage: genJac(inducer, reporter, reporterStd, expression, thetaList), where inducer, reporter, and reporterStd should be iterable")

    # Claim symbols
    thetaSymbols = symbols(thetaList)
    # Original function lambdify
    PstFunc = lambdify((symbols('A'), thetaSymbols), expression, "numpy")
    # Symbolic derivitives
    key = str(expression)
    #if key in jacBase:
    #    jacExp = jacBase[key]    # shortcut if initialized
    jacExp = [diff(expression, x) for x in thetaSymbols]
    # Derivitives lambdify
    jacFunc = [lambdify((symbols('A'), thetaSymbols), exp, "numpy") for exp in jacExp]

    # Avoid 0 in log
    # This has caused big troubles when minimizing
    inducerSize = len(inducer)
    inducer = np.array(inducer)
    inducer = inducer.reshape(inducerSize, -1)    # convert 1D vector to 2D
    flatInducer = inducer.reshape(1, -1)
    mval = float('inf')
    for x in flatInducer[0]:
        if abs(x) > 1E-128 and abs(x) < mval:
            mval = abs(x)
    row, col = inducer.shape
    for i in range(row):
        for j in range(col):
            if abs(inducer[i][j]) <= 1E-128:
                inducer[i][j] = mval / 1E2

    # Construct the function
    def jacobian(theta):
        jac = []
        for i in range(len(thetaList)):
            # derivitive of sse is a linear combinition
            dSSE = 2 * sum(
                [(PstFunc(A[0], theta) - x) * jacFunc[i](A[0], theta)\
                    for A, x in zip(inducer, reporter)]
            )
            jac.append(dSSE)
        return np.array(jac)    # Some scipy optimization funcs, e.g. BFGS, insist on using np.array jacobian

    # Return
    return jacobian

def estimatePara(sse, X0, jacobian = None, constraints = None, method = ModelSolver.Nelder_Mead):
    # Check input
    if type(method) != ModelSolver:
        raise Exception('Usage: paraEstimator(sse, X0, method), where method should be ModelSolver type Enum.')

    # Always first use Nelder_Mead to estimate initial theta
    NMres = minimize(sse, X0, method = "Nelder-Mead", options = {
        "maxiter": 50,
        "disp": False
    })
    newX0 = NMres.x

    # Alternative solvers
    # Nelder_Mead only
    if method == ModelSolver.Nelder_Mead or method == ModelSolver.N_M:
        res = NMres

    # COBYLA
    elif method == ModelSolver.COBYLA:
        res = minimize(sse, newX0, constraints = constraints, method = "COBYLA", options = {
            "tol": 1E-4,
            "maxiter": 50,
            "disp": False
        })

    # BFGS is a Quasi-Newton method, fast and rather accurate
    elif method == ModelSolver.BFGS:
        res = minimize(sse, newX0, jac = jacobian, method = "BFGS", options = {
            "maxiter": 100,
            "gtol": 1E-4,
            "disp": False
        })

    # SLSQP accepts constraints
    elif method == ModelSolver.SLSQP:
        res = minimize(sse, newX0, constraints = constraints, jac = jacobian, method = "SLSQP", options = {
            "maxiter": 100,
            "ftol": 1E-4,
            "disp": False
        })

    # Return
    return res.x, res.fun

def paraPosterior(inducer, reporter, dp, thetaKeys, eqnFunc):
    bayesianModel = pm.Model()
    with bayesianModel:
        thetaR = []
        for para0, key in zip(dp, thetaKeys):
            if key[0] != 'K':
                paraR = pm.Uniform(key, lower=0, upper=4*abs(para0))
            else:
                print('Solving log(K) instead')
                # K can be too small. instead, we calculate the dist of log(K)
                magnitude = abs(np.log(abs(para0)))
                paraR = np.exp(pm.Uniform(key, lower=-magnitude - 4, upper=magnitude + 4))
            thetaR.append(paraR)
        # The expected value for reporter
        mu = eqnFunc(np.squeeze(np.array(inducer)), thetaR)
        sigma = pm.HalfNormal('sigma', 1000)
        g = 10000
        #g = sigma * (abs(mu)**2 + 1)
        # TODO: duo-input
        # observed reporter
        obsReporter = pm.Normal('obsReporter', mu=mu, sigma=g, observed=reporter)
        # MAP
        #mapFit = pm.find_MAP(model=bayesianModel)
        # Sampling
        trace = pm.sample(draws=5000, tune=2000, cores=2, target_accept=0.95)
        # parameter means
        summary = pm.summary(trace)
        inferredTheta = [summary.loc[key, 'mean'] for key in thetaKeys]
        correctedTheta = []
        for key, val in zip(thetaKeys, inferredTheta):
            if key[0] == 'K':
                correctedTheta.append(np.exp(val))
            else:
                correctedTheta.append(val)
        # parameter stdev
        thetaStd = [summary.loc[key, 'sd'] for key in thetaKeys]
        # residue
        reporterHat = np.squeeze(eqnFunc(np.array(inducer), correctedTheta))
        residue = sum((yHat - y)**2 for yHat, y in zip(reporterHat, reporter))
    return correctedTheta, thetaStd, residue, trace    # K is reversed for the ease of plotting

def infoCriteria(sse, theta, reporter, method = ModelCriterion.AIC):
    if method == ModelCriterion.AIC:
        k, n = len(theta), len(reporter)
        AIC = 2 * k + n * np.log(sse / n)
        # AIC correction when small sample size
        #AICc = AIC + 2 * k * (k + 1) / (n - k - 1)
        return AIC
    else:
        return sse

def __fitModelWrapper(paras):
    '''
    Split parameters for __fitModel()
    '''
    return __fitModel(*paras)

def __fitModel(modelType, modelSpecs,
                inducer, reporter, reporterStd,
                modelSolver = ModelSolver["Nelder_Mead"], modelCriterion = ModelCriterion["AIC"]):
    '''
    Paralleled version of model fitting
    '''

    # Is this repression or activation model?
    repression = False
    if ModelSpec.Repression in modelSpecs and ModelSpec.Inducer_Repression not in modelSpecs:
        repression = True
    elif ModelSpec.Repression not in modelSpecs and ModelSpec.Inducer_Repression in modelSpecs:
        repression = True

    startTime = time.time()
    if modelSolver != ModelSolver.MCMC:
        # Generate the model
        (model, eqnStr, eqnSym), thetaList, constraints = modelBase.genModel(modelType, modelSpecs)

        # Get SSE
        sse = genSSE(inducer, reporter, reporterStd, model, thetaList)
        # Also get jacobian of SSE
        if modelSolver == ModelSolver.Nelder_Mead or modelSolver == ModelSolver.N_M:
            jacobian = None
        else:
            jacobian = genJac(inducer, reporter, reporterStd, eqnSym, thetaList)

        # Generate default parameters
        dp = modelBase.defaultPara(thetaList, inducer, reporter, repression = repression)
        # Parameterization
        inferredTheta, residue = estimatePara(
                sse, dp, jacobian, constraints, method = modelSolver
        )

    # MCMC
    else:
        # TODO: the initial value is a bit messy
        # generate sse for Nelder_Mead
        (model, eqnStr, eqnSym), thetaList, constraints = modelBase.genModel(modelType, modelSpecs)
        sse = genSSE(inducer, reporter, reporterStd, model, thetaList)
        # default parameters
        dp = modelBase.defaultPara(thetaList, inducer, reporter, repression = repression)
        # better initial theta
        NMres = minimize(sse, dp, method = "Nelder-Mead", options = {
            "maxiter": 50,
            "disp": False
        })
        newX0 = NMres.x
        (eqnFunc, eqnStr, eqnSym), thetaList, constraints = modelBase.genEquation(modelType, modelSpecs)
        inferredTheta, thetaStd, residue, trace = paraPosterior(inducer, reporter, newX0, thetaList, eqnFunc)

    duration = time.time() - startTime

    # Calculate Infomation Criteria
    IC = infoCriteria(residue, inferredTheta, reporter, modelCriterion)

    # Compose model meta
    meta = Solution(
            modelType, modelSpecs,
            eqnSym,
            thetaList, inferredTheta,
            residue, IC)
    # if MCMC, also adds trace and thetaStd
    if modelSolver == ModelSolver.MCMC:
        meta.thetaStd = thetaStd
        meta.trace = trace

    return meta, duration

def selectModel(inducer, reporter, reporterStd = None,
        modelSolver = ModelSolver.N_M, modelSet = ModelSet.Simple_Inducible_Promoter, modelCriterion = ModelCriterion.AIC,
        quiet=True, parallel = True):
    '''
    Select the best model interpreting given inducer/reporter characterization data.
    Notice that fitting of alternative models is paralled. Thus, selectModel() cannot be run multiprocessed!
    '''

    # Generate model candicates (alternative mechanistic hypotheses)
    modelSet = modelBase.genModelSet(modelSet)

    metas = []
    minIC = float("inf")
    count = 0    # number of model being solved

    if parallel and modelSolver != ModelSolver.MCMC:    # mcmc conflicts with parallel computation
        # parallel model fitting
        cpuCount = min(len(modelSet), int(mp.cpu_count()))    # Num of cpu cores
        pool = mp.Pool(cpuCount)    # A pool of processes

        for meta, duration in pool.imap_unordered(__fitModelWrapper,
                                                [(*mk, inducer, reporter, None, modelSolver, modelCriterion) for mk in modelSet]):
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

    else:
        for mk in modelSet:
            meta, duration = __fitModel(*mk, inducer, reporter, None, modelSolver, modelCriterion)
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

    return bestModelMeta, metas

def printModel(meta):
    print('%s\nType = %s'%(pretty(meta.expression, use_unicode=False), meta.modelType.name))
    print('Specs = ' + ', '.join([spec.name for spec in meta.modelSpecs]))
    if meta.thetaStd != None:
        print('Parameters:\n' +\
                ', '.join(['%s = %.4e; std: %.4e'%(key, val, std) for key, val, std in zip(meta.thetaKey, meta.thetaVal, meta.thetaStd)]))
    else:
        print('Parameters:\n' +\
                ', '.join(['%s = %.4e'%(key, val) for key, val in zip(meta.thetaKey, meta.thetaVal)]))
    print('Residue = %.4e\nIC = %.4f'%(meta.residue, meta.IC))
    stdout.flush()
