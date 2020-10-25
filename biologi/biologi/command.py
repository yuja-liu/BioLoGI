#!/usr/bin/env python

# -*- coding: utf-8 -*- 

'''
Author: Yujia Liu rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.30

main.py is the entrance script. It handles the command line parameters, and returns the best model as well as its parameter(s).
'''

import sys, getopt
from biologi import ioUtils, inference, modelBase, plotUtils
from biologi.modelBase import ModelSet
from biologi.inference import ModelSolver, ModelCriterion
import time
import pymc3 as pm

helpInfo = '''
BioLogic

A command line tool to infer the optimal model, parameters, and mechanistic hypotheses for a biological *node*, given characterization data.

Nomenclatures: "Reporter" refers to the output of a circuit, e.g., a fluorescense protein or a gene whose expression being regulated.
"Inducer" is the term for input signal, which is often a small molecule, but can also be transcription factors.
Take a Lac operon for an example. The "inducer" can be IPTG or LacI, depends on the hypotheses, while "reporter" can be *lacZYA*, or GFP, if we have a recombinated strain here.

Usage:

  biologic

    -c, --config <file> : provide a file of necessary information.
    -h                  : show this message
    -q                  : quiet mode
    -u, --units <units> : Specify the units for both the inducer and the receiver
    --inducer-units     : Assign the units of the inducer
    --reporter-units    : Assign the units of the reporter
    --inducer-tags      : Labels for the inducers. Multiple values separated by commas are encouraged, e.g., "LacI, IPTG"
    --replicate-tags    : Names for each experiment
    --reporter-tag      : The label for the reporter, e.g., GFP
    --figure            : Designate the output path for the figure
    -f, --file          : Designate the path for characterization data
    -c, --config        : Use a configuration file to specify details of the model, and/or the solver
    -q, --quiet         : Not showing alternetive models. Verbose is default.
    -h, --help          : Show this text
'''

def main():
    # parse the arguments
    opts, args = getopt.getopt(sys.argv[1:], "-h-c:-q", ["help", "config=", "quiet"])
    configPath, quiet = "", False
    for optKey, optVal in opts:
        if optKey in ("-h", "--help"):
            print(helpInfo)
            exit()
        elif optKey in ("-c", "--config"):
            configPath = optVal
        elif optKey in ("-q", "--quiet"):
            quiet = True
    if configPath == "":
        print(helpInfo)
        exit()
    else:
        config = ioUtils.readConfig(configPath)

        # Parse the config
        if "input" not in config.keys():
            raise Exception("[Tentative] the input field is required in configuration")
        inducerUnits = config.get('input', {}).get('inducerUnits', None)
        reporterUnits = config.get('input', {}).get('reporterUnits', None)
        inducerTags = config.get('input', {}).get('inducerTags', None)
        replicateTags = config.get('input', {}).get('replicateTags', None)
        reporterTag = config.get('input', {}).get('reporterTag', None)
        chassisType = config.get('input', {}).get('chassisType', None)
        filePath = config.get('input', {}).get('file', None)
        modelSolver = ModelSolver[
                config.get('solver', {}).get('paraEstimator', None)]
        modelCriterion = ModelCriterion[
                config.get('solver', {}).get('infoCriteria', None)]
        modelSet = ModelSet[
                config.get('solver', {}).get('modelSet', None)]
        figurePath = config.get('output', {}).get('figureFile', None)
        sbmlPath = config.get('output', {}).get('sbmlFile', None)
        sbolPath = config.get('output', {}).get('sbolFile', None)

        # inducer/reporter matrix from file
        inSplit, outSplit, stdSplit = ioUtils.readSplitInput(filePath, inducerTags, replicateTags, reporterTag)
        inducer, reporter, reporterStd = [], [], []
        for i, o, s in zip(inSplit, outSplit, stdSplit):
            inducer += i
            reporter += o
            reporterStd += s

        startTime = time.time()
        #modelBase.initializeModelBase()
        bestModel, _ = inference.selectModel(
                inducer, reporter, reporterStd,
                modelSolver, modelSet, modelCriterion,
                quiet=False, parallel=True)
        print("Total Run Time:", time.time() - startTime)    # Note this includes plotting elapse

        # print model selection
        print('========================================')
        print('Model Choice:')
        inference.printModel(bestModel)

        # Plotting
        thetaDict = {key: val for key, val in zip(bestModel.thetaKey, bestModel.thetaVal)}
        modelFunc = modelBase.genModel(bestModel.modelType, bestModel.modelSpecs)[0][0]
        if modelSolver == ModelSolver.MCMC:    # plot MCMC trace
            pm.traceplot(bestModel.trace)
        plotUtils.plotModel2D(
                inSplit, outSplit, stdSplit, 
                thetaDict, modelFunc,
                inputTag = inducerTags[0], 
                outputTag = reporterTag, 
                inputUnits = inducerUnits, 
                outputUnits = reporterUnits, 
                logScale = True,
                save_fig = figurePath)
        print('========================================')
        print('Successfully plotted model fitting')

        # output to sbml
        if sbmlPath != None:
            ioUtils.printSBML(
                    bestModel.modelType,
                    bestModel.expression,
                    {key: val for key, val in zip(bestModel.thetaKey, bestModel.thetaVal)},
                    inducerTags,
                    reporterTag,
                    inducerUnits,
                    reporterUnits,
                    sbmlPath,
                    chassisType=chassisType
            )
            print('========================================')
            print('Successfully generated SBML')

        # show as SBOL visuals
        plotUtils.plotDNACircuit(
                bestModel.modelType,
                bestModel.modelSpecs,
                inducerTags,
                reporterTag,
                sbolPath
        )
        print('========================================')
        print('Successfully generated SBOL visuals')

# command line entrance
if __name__ == "__main__":
    main()
