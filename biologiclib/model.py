# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.30

model.py describes the models to infer parameters and to select optimal mechanistic model. For actual mechanistic models please refer to modelBase.ph
'''

import warnings
import numpy as np
from enum import Enum
from scipy.optimize import minimize
from functools import reduce
from biologiclib import ioUtils, modelBase, plotUtils
from biologiclib.modelBase import ModelType, ModelSpec, ModelSet

# ModelSolver defines the possible procedures to estimate the parameters of a single mechanistic model
ModelSolver = Enum("ModelSolver", ("Nelder_Mead", "N_M"))

# ModelCriterion defines the standard to evaluate a given model{expression and the parameters}, basically Infromation Criteria
ModelCriterion = Enum("ModelCriterion", ("AIC"))

# thetaList is necessary for buidling sse func. However as thetaList should be considered as a part of the model, no additional info. is needed.
def genSSE(data, model, thetaList):
    # Check data type
    if not("input" in data.keys() and "output" in data.keys() and "std" in data.keys()):
        raise Exception('Fields "input", "output", and "std" are required in data')

    def SSE(X):
        # Check # parameter
        if len(X) != len(thetaList):
            raise Exception('Incorrect number of parameters. Expected %d but %d parameters are given'%(len(thetaList), len(X)))
        # Composing theta **Note that the order of X is important**
        theta = {key: value for key, value in zip(thetaList, X)}
        # the total SSE is the sum of the SSE of each experiemnt
        sseTotal = 0
        # the input should be tranpose, to a vertical vector
        mu = model(np.array(data["input"]), theta)
        for row in data["output"]:
            sse = sum([(m - x) ** 2 for m, x in zip(mu, row)])
            sseTotal += sse
        return sseTotal

    return SSE

def paraEstimator(sse, X0, method = ModelSolver.Nelder_Mead):
    # Check input
    if type(method) != ModelSolver:
        raise Exception('Usage: paraEstimator(sse, X0, method), where method should be ModelSolver type Enum.')
    # Nelder-Mead, do not need bounds nor gradients
    if method == ModelSolver.Nelder_Mead:
        res = minimize(sse, X0, method = "Nelder-Mead", options = {
            "xatol": 1E-6, "disp": False
        })
        return res.x, res.fun
    else:
        warnings.warn("The required solver %s has not been implemented yet!"%method.__name__)
        return -1, -1

def infoCriteria(sse, theta, data, method = ModelCriterion.AIC):
    if method == ModelCriterion.AIC:
        k, n = len(theta), len(data[0]) * len(data)
        AIC = 2 * k + n * np.log(sse / n)
        # AIC correction when small sample size
        #AICc = AIC + 2 * k * (k + 1) / (n - k - 1)
        return AIC
    else:
        return sse

def selectModel(config, quiet=False):
    # If the ontology of codes are well-defined, then the higher the logic, the less comments will be required
    data = ioUtils.readInput(config)
    modelSet = modelBase.genModelSet(ModelSet[config["modelSet"]])
    minIC, bestModel, bestModelMeta = float('inf'), None, None
    for (i, modelKey) in enumerate(modelSet):
        exp, model, thetaList = modelBase.genModel(*modelKey)

        if not quiet:
            print('#%d model calculating...'%(i+1))
            print('%s\nType = %s'%(exp, modelKey[0].name))
            print('Specs = ' + ', '.join([spec.name for spec in modelKey[1]]))

        sse = genSSE(data, model, thetaList)
        inferredTheta, residue = paraEstimator(sse, modelBase.defaultPara(thetaList),
                ModelSolver[config["methods"]["paraEstimator"]])
        inferredThetaDict = {key: val for key, val in zip(thetaList, inferredTheta)}
        IC = infoCriteria(residue, inferredTheta, data["output"],
                ModelCriterion[config["methods"]["infoCriteria"]])

        if not quiet:
            print('Parameters:\n' +\
                    ', '.join([key + ' = ' + str(val) for key, val in inferredThetaDict.items()]))
            print('Residue = %.4f\nIC = %.4f'%(residue, IC))
            print('=============================================')

        if IC < minIC:
            minIC = IC
            bestModel = model
            bestModelMeta = (exp, modelKey, inferredThetaDict, residue, IC)

    # Print model selection
    print('Model Choice:')
    print('%s\nType = %s'%(bestModelMeta[0], bestModelMeta[1][0].name))
    print('Specs = ' + ', '.join([spec.name for spec in bestModelMeta[1][1]]))
    print('Parameters:\n' +\
            ', '.join([key + ' = ' + str(val) for key, val in bestModelMeta[2].items()]))
    print('Residue = %.4f\nIC = %.4f'%bestModelMeta[3:])

    # Plotting
    figPath = config["figPath"] if "figPath" in config.keys() else None
    # Handling input: expand and transpose
    # TODO: these codes are arwful, try to rebuild
    X = []
    for i in range(len(data["output"])):
        X += data["input"]
    plotUtils.plotModel2D(X,
            reduce(lambda x, y: x + y, data["output"]),
            reduce(lambda x, y: x + y, data["std"]),
            bestModelMeta[2],
            bestModel,
            config["tags"]["inputTags"][0], config["tags"]["outputTag"],
            config["units"]["inputUnits"], config["units"]["outputUnits"],
            figPath)

    return bestModelMeta
