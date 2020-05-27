# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.04.02

Ploting utilities
'''

import dnaplotlib as dpl
from matplotlib import pyplot as plt
import numpy as np
from biologiclib.modelBase import ModelType, ModelSpec
import matplotlib
matplotlib.use('TkAgg')    # dnaplotlib changes the backend to 'Agg'

def plotModel2D(X, mu, std, theta, model, inputTag = "Inducer", outputTag = "FP", inputUnits = "M", outputUnits = "AU", save_fig = None):
    # Check input
    try:
        if len(X) != len(std):
            raise Exception("Usage: plotModel(X, std, mu, theta, model). The dimension of X (%d), std (%d), and mu (%d) should be equal"\
                    %(len(X), len(std), len(mu)))
        if len(X) != len(mu):
            raise Exception("Usage: plotModel(X, std, mu, theta, model). The dimension of X (%d), std (%d), and mu (%d) should be equal"\
                    %(len(X), len(std), len(mu)))
    except ValueError:
        raise Exception("Usage: plotModel(X, std, mu, theta, model), where X, std and mu should be iterators")

    try:
        len(X[0])    # a hack to raise error if 1-dimensional
        X1D = [x[0] for x in X]    # Extract the 1st dimension
    except TypeError:
        X1D = X
    # 10% extrapolation on both ends, on a log scale
    # Handling 0 input
    lower = min(X1D)
    if lower < 1E-12:    # consider as read zero
        Xsorted = sorted(X1D)
        for x in Xsorted:
            if x >= 1E-12:
                Xmin = np.log(x)
                break
    else:
        Xmin = np.log(lower)
    Xmax = np.log(max(X1D))
    try:
        plotBounds = (Xmin - 0.2 * (Xmax - Xmin), Xmax + 0.2 * (Xmax - Xmin))
    except NameError:    # Xmin does not exist, which is caused by all-zero inducer
        raise Exception("All-zero inducer input is not accepted")
    plotRange = np.exp(np.arange(*plotBounds, 1E-2))
    if lower < 1E-12:    # give back the 0 at front
        plotRange = np.insert(plotRange, 0, 0.0)
    Y = model(plotRange.reshape(-1, 1), theta)

    fig, ax = plt.subplots(figsize=(8, 6))
    #ax.plot(X1D, mu, color="blue", marker='^', linestyle=None)
    ax.plot(plotRange, Y, 'k-')
    ax.errorbar(X1D, mu, yerr=std, fmt='o', color="blue")
    ax.set_xlabel("[%s]/%s"%(inputTag, inputUnits))
    ax.set_ylabel("[%s]/%s"%(outputTag, outputUnits))
    ax.set_xscale("symlog", linthreshx = np.exp(plotBounds[0]))
    if save_fig is not None:
        # Then save as a file
        plt.savefig(save_fig, dpi=300, format="png")
    else:
        plt.show()

def plotDNACircuit(modelType, modelSpecs, inducerTags, reporterTag, figPath=None):
    '''
    plotDNACircuit(modelType, modelSpecs, inducerTags, reporterTag, figPath=None)
    Renders the inferred model of the biological node as genetic circuits.
    The plot conforms SBOL visual standards
    '''

    # color scheme
    palatte = ['#F79F1F', '#A3CB38', '#1289A7', '#D980FA', '#B53471']
    np.random.shuffle(palatte)

    # components
    design = []
    lw = 1.0
    cur = 0
    if modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        inPart = {'type': 'CDS', 'name': 'inPart', 'fwd': True,
                'opts': {'linewidth': lw, 'color': palatte[cur], 'edge_color': 'black', 'x_extent': 24,
                    'label': inducerTags[0], 'label_style': 'italic', 'label_color': 'black', 'label_size': 15,    # label_size controls font size
                    'label_x_offset': -3, 'label_y_offset': 0
                }
        }
        cur += 1
        inPromoter = {'type': 'Promoter', 'name': 'pIn', 'fwd': True,
                'opts': {'linewidth': lw, 'color': 'black', 'label': ''}}
        design += [inPromoter, inPart]
    elif modelType == ModelType.Duo_Input_Node:
        inPart1 = {'type': 'CDS', 'name': 'inPart1', 'fwd': True,
                'opts': {'linewidth': lw, 'color': palatte[cur], 'edge_color': 'black', 'x_extent': 24,
                    'label': inducerTags[0], 'label_style': 'italic', 'label_color': 'black', 'label_size': 15,
                    'label_x_offset': -3, 'label_y_offset': 0
                }
        }
        cur += 1
        inPromoter1 = {'type': 'Promoter', 'name': 'pIn1', 'fwd': True,
                'opts': {'linewidth': lw, 'color': 'black', 'label': ''}}
        inPart2 = {'type': 'CDS', 'name': 'inPart2', 'fwd': True,
                'opts': {'linewidth': lw, 'color': palatte[cur], 'edge_color': 'black', 'x_extent': 24,
                    'label': inducerTags[1], 'label_style': 'italic', 'label_color': 'black', 'label_size': 15,
                    'label_x_offset': -3, 'label_y_offset': 0
                }
        }
        cur += 1
        inPromoter2 = {'type': 'Promoter', 'name': 'pIn2', 'fwd': True,
                'opts': {'linewidth': lw, 'color': 'black', 'label': ''}}
        design += [inPromoter1, inPart1, inPromoter2, inPart2]
    outPart = {'type': 'CDS', 'name': 'outPart', 'fwd': True,
            'opts': {'linewidth': lw, 'color': palatte[cur], 'edge_color': 'black', 'x_extent': 24,
                'label': reporterTag, 'label_style': 'italic', 'label_color': 'black', 'label_size': 15,
                'label_x_offset': -3, 'label_y_offset': 0
            }
    }
    cur += 1
    outPromoter = {'type': 'Promoter', 'name': 'pOut', 'fwd': True,
            'opts': {'linewidth': lw, 'color': 'black', 'label': ''}}
    design += [outPromoter, outPart]

    # regulations
    ac = 25    # arc height
    if modelType == ModelType.Constant:
        regs = []
    elif modelType == ModelType.Inducible or modelType == ModeType.Single_Input_Node:
        isActivated = ModelSpec.Activation in modelSpecs and ModelSpec.Inducer_Repression not in modelSpecs
        reg = {'type': 'Activation' if isActivated else 'Repression', 'from_part': inPart, 'to_part': outPromoter,
                'opts': {'color': 'black', 'linewidth': lw, 'linestyle': '-' if ModelSpec.Inducer not in modelSpecs else '--', 'arc_height': ac}}
        regs = [reg]
    elif modelType == ModelType.Duo_Input_Node:
        truthTable = {
                'AND': (True, True),
                'OR': (True, True),
                'IMPLY12': (True, False),
                'IMPLY21': (False, True),
                'NAND': (False, False),
                'NOR': (False, False)
        }
        isActivated = truthTable[modelSpecs[0].name]
        reg1 = {'type': 'Activation' if isActivated[0] else 'Repression', 'from_part': inPart1, 'to_part': outPromoter,
                'opts': {'color': 'black', 'linewidth': lw, 'linestyle': '-', 'arc_height': ac}}
        reg2 = {'type': 'Activation' if isActivated[1] else 'Repression', 'from_part': inPart2, 'to_part': outPromoter,
                'opts': {'color': 'black', 'linewidth': lw, 'linestyle': '-', 'arc_height': ac}}
        regs = [reg1, reg2]

    # plot
    fig, ax = plt.subplots(figsize=(6, 3))
    dr = dpl.DNARenderer()
    start, end = dr.renderDNA(ax, design, dr.SBOL_part_renderers(),
        regs=regs, reg_renderers=dr.std_reg_renderers())
    ax.set_xlim([start, end])
    ax.set_ylim([-10, 30])
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    if figPath != None:
        plt.savefig(figPath, dpi=300, format='png')
    else:
        plt.show()
