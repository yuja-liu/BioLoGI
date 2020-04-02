# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.29

Input & output utilities
'''

import json
import csv

# Read input configurations, e.g. units and tags
def readConfig(filePath):
    configFile = open(filePath, 'r')
    try:
        config = json.loads(configFile.read())
        configFile.close()
    except Error:
        raise Exception("Error while parsing json file. Please check the grammar.")
    # Check required fields
    for key in ["units", "tags", "dataPath", "methods"]:
        if key not in config.keys():
            raise Exception('The fields "units", "tags", "methods", and "dataPath" are required in config file')
    if not("inputUnits" in config["units"].keys() and "outputUnits" in config["units"].keys()):
        raise Exception('The fields "inputUnits" and "outputUnits" are required in "units" section.')
    if not("inputTags" in config["tags"].keys() and "outputTag" in config["tags"].keys()\
            and "experimentTags" in config["tags"].keys()):
        raise Exception('The fields "inputTags", "outputTag" and "experimentTags" are required in "tags" section.')
    if type(config["tags"]["inputTags"]) != list or type(config["tags"]["outputTag"]) != str\
            or type(config["tags"]["experimentTags"]) != list:
        raise Exception('The field "inputTags" or "experimentTags" requires a list, while the field "outputTag" require a string.')
    
    return config

def readInput(config):
    # Input file should be tab-separated
    fi = open(config["dataPath"], 'r')
    lines = csv.reader(fi, dialect="excel-tab")

    # Check columns
    # #input = N, #experiment = M
    # It is expected to have N columns of input, followed by M columns of output and M columns of deviations
    expectedColumnNum = len(config["tags"]["inputTags"]) + len(config["tags"]["experimentTags"]) * 2
    firstLine = lines.__next__()
    if len(firstLine) != expectedColumnNum:
        raise Exception('Unexpected input layout. Please check the number of columns')
    expectedInputNum, expectedExpNum = len(config["tags"]["inputTags"]), len(config["tags"]["experimentTags"])
    readInputTags, readExpTags, readStdTags = firstLine[:expectedInputNum], firstLine[expectedInputNum:expectedInputNum + expectedExpNum],\
            firstLine[expectedInputNum + expectedExpNum:expectedInputNum + 2 * expectedExpNum]
    for tag1, tag2 in zip(readExpTags, readStdTags):
        if tag1 != tag2:
            raise Exception('Unexpected input layout. "output" columns and "deviation" columns should have the same headers and in the same order')
    for key in config["tags"]["inputTags"]:
        if key not in readInputTags:
            raise Exception('Unexpected input layout. "input" columns should agree with config file.')
    for key in config["tags"]["experimentTags"]:
        if key not in readExpTags:
            raise Exception('Unexpected input layput. "output" or "deviation" columns should agree with config file.')

    # Read data (to array)
    flatLines = []
    for line in lines:
        flatLines.append(line)
    try:
        tmp = [[float(line[i]) for line in flatLines]for i in range(len(firstLine))]    #transpose
    except ValueError:
        raise Exception('The input file requires all numeric values, which is not the case.')
    dataInput, dataOutput, dataStd = tmp[:expectedInputNum], tmp[expectedInputNum:expectedInputNum + expectedExpNum],\
            tmp[expectedInputNum + expectedExpNum:expectedInputNum + 2 * expectedExpNum]
    # Transpose back the input
    dataInput = [[row[i] for row in dataInput] for i in range(len(dataInput[0]))]

    return {"input": dataInput, "output": dataOutput, "std": dataStd}
