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
    
    return config

# Tags are used to check the input file
def readInput(filePath, inducerTags, replicateTags, reporterTag):
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
        inducer += dataInput
        receiver += rRow
        receiverStd += sRow

    return inducer, receiver, receiverStd
