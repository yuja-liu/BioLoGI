# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.29

Input & output utilities
'''

import json
import csv
from biologiclib.modelBase import ModelType
from sympy import symbols
# customed version of simplesbml, to add units
from biologiclib import simplesbml

# Read input configurations, e.g. units and tags
def readConfig(filePath):
    configFile = open(filePath, 'r')
    try:
        config = json.loads(configFile.read())
        configFile.close()
    except Error:
        raise Exception("Error while parsing json file. Please check the grammar.")
    
    return config

# Tags are used to check the input file
def readSplitInput(filePath, inducerTags, replicateTags, reporterTag):
    # Input file should be tab-separated
    fi = open(filePath, 'r')
    lines = csv.reader(fi, dialect="excel-tab")

    # Check columns
    # #input = N, #experiment = M
    # It is expected to have N columns of input, followed by M columns of output and M columns of deviations
    expectedColumnNum = len(inducerTags) + len(replicateTags) * 2
    firstLine = lines.__next__()
    if len(firstLine) != expectedColumnNum:
        raise Exception('Unexpected input layout. Please check the number of columns')
    expectedInputNum, expectedExpNum = len(inducerTags), len(replicateTags)
    readInputTags, readExpTags, readStdTags = firstLine[:expectedInputNum], firstLine[expectedInputNum:expectedInputNum + expectedExpNum],\
            firstLine[expectedInputNum + expectedExpNum:expectedInputNum + 2 * expectedExpNum]
    for tag1, tag2 in zip(readExpTags, readStdTags):
        if tag1 != tag2:
            raise Exception('Unexpected input layout. "output" columns and "deviation" columns should have the same headers and in the same order')
    for key in inducerTags:
        if key not in readInputTags:
            raise Exception('Unexpected input layout. "input" columns should agree with config file.')
    for key in replicateTags:
        if key not in readExpTags:
            raise Exception('Unexpected input layput. "output" or "deviation" columns should agree with config file.')

    # Read data (to array)
    flatLines = []
    for line in lines:
        flatLines.append(line)
    fi.close()    # close the file, since lines depends on it
    try:
        tmp = [[float(line[i]) for line in flatLines]for i in range(len(firstLine))]    #transpose
    except ValueError:
        raise Exception('The input file requires all numeric values, which is not the case.')
    dataInput, dataOutput, dataStd = tmp[:expectedInputNum], tmp[expectedInputNum:expectedInputNum + expectedExpNum],\
            tmp[expectedInputNum + expectedExpNum:expectedInputNum + 2 * expectedExpNum]
    # Transpose back the input
    dataInput = [[row[i] for row in dataInput] for i in range(len(dataInput[0]))]
    # reshape
    inducer, receiver, receiverStd = [], [], []
    for rRow, sRow in zip(dataOutput, dataStd):
        inducer.append(dataInput)
        receiver.append(rRow)
        receiverStd.append(sRow)

    return inducer, receiver, receiverStd

def readInput(filePath, inducerTags, replicateTags, reporterTag):
    '''
    inducer, receiver, receiverStd = readInput(filePath, inducerTags, replicateTags, reporterTag)
    Read characterization data from file, returning a single vector of input (inducer), output (reporter) and the standard deviation of output.
    '''
    inSplit, outSplit, stdSplit = readSplitInput(filePath, inducerTags, replicateTags, reporterTag)
    
    # concatenate
    inducer, reporter, reporterStd = [], [], []
    for i, o, s in zip(inSplit, outSplit, stdSplit):
        inducer += i
        reporter += o
        reporterStd += s
    return inducer, reporter, reporterStd

def printSBML(modelType, eqn, theta, inducerTags, reporterTag,\
        inducerUnits, reporterUnits, outputPath, **kargs):
    '''
    printSBML(modelType, eqn, theta, inducerTags, reporterTag, inducerUnits, reporterUnit, outputPath, **kargs)
    Output a parameterized model to SBML
    '''

    model = simplesbml.sbmlModel()
    # assign compartment
    try:
        chassisType
    except NameError:
        chassisType = 'prokaryotes'    # default
    volumeMap = {
            "prokaryotes": 1E-15,
            "eukaryotes": 1E-12,
            "cell_free": 1E-6
    }    # units in L
    model.addCompartment(volumeMap[chassisType], "comp1")

    # register inducer & reporter to species
    # remove empty characters in tags
    emptyChar = [' ', '\t', '\n', '\r']
    for char in emptyChar:
        reporterTag = reporterTag.replace(char, '_')
    for i in range(len(inducerTags)):
        for char in emptyChar:
            inducerTags[i] = inducerTags[i].replace(char, '_')
    # if flanked by brackets, the unit of species will be automatically assigned to M instead of mole
    # also, the actual id is striped of brackets
    # check out https://simplesbml.readthedocs.io/en/latest/#simplesbml.sbmlModel.addSpecies
    for it in inducerTags:    # all initialize to 0
        model.addSpecies('[' + it + ']', 0.0, comp="comp1")
    model.addSpecies('[' + reporterTag + ']', 0.0, comp="comp1")

    # add parameters
    unitsMap = {
            "alpha": "molar_per_second", "alpha_1": "molar_per_second", "alpha_2": "molar_per_second", "alpha_3": "molar_per_second", "alpha_4": "molar_per_second",
            "b": "molar_per_second",
            "K": "molar", "K_I": "molar", "K_1": "molar", "K_2": "molar",
            "n": "dimensionless", "n_1": "dimensionless", "n_2": "dimensionless",
            "beta": "per_second"
    }    # notice though in Hill, alpha is dimensionless, in ODE it is of per second
    for key, val in theta.items():
        model.addParameter(key, val, unitsMap[key])
    # beta is set to 1.0 as default
    model.addParameter('beta', 1.0, unitsMap['beta'])

    # add generation rate rules
    degradTerm = '-beta*' + reporterTag
    if modelType == ModelType.Constant:
        genTerm = 'alpha'
    elif modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        A = symbols('A')
        inducer = symbols(inducerTags[0])
        genTerm = str(eqn.subs(A, inducer))
    elif modelType == ModelType.Duo_Input_Node:
        A_1, A_2 = symbols('A_1 A_2')
        inducer1, inducer2 = symbols(inducerTags)
        genTerm = str(eqn.subs((A_1, inducer1), (A_2, inducer2)))
    rateRule = genTerm + degradTerm
    rateRule = rateRule.replace('**', '^')    # a hack for exponential problem
    model.addRateRule(reporterTag, rateRule)

    # write to file
    fo = open(outputPath, 'w')
    fo.write(model.toSBML())
