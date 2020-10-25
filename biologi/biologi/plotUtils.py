# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.04.02

Ploting utilities
'''

import dnaplotlib as dpl
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
from biologi.modelBase import ModelType, ModelSpec, genModel
from biologi.inference import Solution
import matplotlib
try:
    matplotlib.use('TkAgg')    # dnaplotlib changes the backend to 'Agg'
except ImportError:    # headless platform
    pass

def plotHelper3D(inducer, reporter, modelMeta, inducer_name = ('', ''), reporter_name = '', 
        reporterStd = None, ax = None, logScale = True, **plotKW):
    if type(modelMeta) != Solution:
        raise Exception("The argument modelMeta should be of type inference.Solution. Abort")

    # convert input to np.array
    inducer = np.array(inducer)
    reporter = np.array(reporter)
    # convert inducer to row vector
    inducer = inducer.transpose()
    # get plot axis if needed
    if ax == None:
        fig = plt.figure(figsize = (8, 6))
        ax = fig.gca(projection='3d')
    # scatter
    if 's' not in plotKW:
        plotKW['s'] = 10
    ax.scatter(*inducer, reporter, **plotKW)
    # rotate to let x=0, y=0 pointing out
    ax.view_init(30, -135)
    ax.set_xlabel(inducer_name[0])
    ax.set_ylabel(inducer_name[1])
    ax.set_zlabel(reporter_name)

    # surface
    # get plot range
    def get_plot_range(X):
        MARGIN = 0.1
        x_range = max(X) - min(X)
        x_min = max(0, min(X) - x_range * MARGIN)
        x_max = max(X) + x_range * MARGIN
        return x_min, x_max
    grid_x = np.linspace(*get_plot_range(inducer[0]), 50)
    grid_y = np.linspace(*get_plot_range(inducer[1]), 50)
    grid = np.meshgrid(grid_x, grid_y)
    grid_rearrange =\
            np.vstack((grid[0].reshape(1, -1, order = 'C'),
            grid[1].reshape(1, -1, order = 'C'))).transpose()
    # dictionize theta
    theta = {key: val for key, val in zip(modelMeta.thetaKey, modelMeta.thetaVal)}
    # get a function for the model
    modelFunc = genModel(modelMeta.modelType, modelMeta.modelSpecs)[0][0]
    mu = modelFunc(grid_rearrange, theta)
    mu = mu.reshape(grid[0].shape, order = 'C')
    ax.plot_surface(*grid, mu,
            cmap = cm.coolwarm,
            linewidth = 0,
            antialiased = True,
            alpha = 0.5)
    # axis log scale
    if logScale:
        ax.set_xscale('symlog', linthresh = 1E-3)
        ax.set_yscale('symlog', linthresh = 1E-3)
    plt.show()

def plotHelper(inducer, reporter, modelMeta, inducer_name = '', reporter_name = '', reporterStd = None, ax = None, logScale = True, **plotKW):
    # modelMeta should be of type Solution
    if type(modelMeta) != Solution:
        raise Exception("The argument modelMeta should be of type inference.Solution. Abort")

    # reshape inducer & reporter
    inducer, reporter = np.array(inducer), np.array(reporter)
    inducer = inducer[:, 0]    # extract only 1 dimension
    inducer, reporter = inducer.reshape(1, -1, 1), reporter.reshape(1, -1)
    # also give a zero std deviation
    std = np.zeros(reporter.shape)
    # dictionize theta
    theta = {key: val for key, val in zip(modelMeta.thetaKey, modelMeta.thetaVal)}
    # get a function for the model
    modelFunc = genModel(modelMeta.modelType, modelMeta.modelSpecs)[0][0]

    plotModel2D(inducer, reporter, std,
            theta,
            modelFunc,
            inputTag = inducer_name,
            outputTag = reporter_name,
            inputUnits = 'AU',
            outputUnits = 'AU',
            ax = ax,
            logScale = logScale,
            **plotKW)

def plotModel2D(X, Y, std, theta, model, 
        inputTag = "Inducer", outputTag = "FP", inputUnits = "M", outputUnits = "AU", 
        logScale = True, save_fig = None, ax=None, **plotKW):
    # TODO: the datastructure required by this functions is very different from others.
    # A reconstruction is recommended
    # Check input
    try:
        if len(X) != len(std):
            raise Exception("Usage: plotModel(X, std, Y, theta, model). The dimension of X (%d), std (%d), and Y (%d) should be equal"\
                    %(len(X), len(std), len(Y)))
        if len(X) != len(Y):
            raise Exception("Usage: plotModel(X, std, Y, theta, model). The dimension of X (%d), std (%d), and Y (%d) should be equal"\
                    %(len(X), len(std), len(Y)))
    except ValueError:
        raise Exception("Usage: plotModel(X, std, Y, theta, model), where X, std and Y should be iterators")

    # color scheme
    palette = ['#F79F1F', '#A3CB38', '#1289A7', '#D980FA', '#B53471', '#EE5A24', '#009432', '#0652DD', '#9980FA', '#833471']
    np.random.shuffle(np.array(palette))

    nonZeros = []
    X = np.squeeze(np.array(X), axis=2)
    try:    # check inducer 2D
        len(X[0])
    except TypeError:
        raise Exception("X should be 2D array")
    zeroFlag = False
    for row in X:
        for x in row:
            if x > 1E-12:
                nonZeros.append(np.log(x))
            else:
                zeroFlag = True
    Xmin, Xmax = min(nonZeros), max(nonZeros)

    if logScale:
        try:
            plotBounds = (Xmin - 0.2 * (Xmax - Xmin), Xmax + 0.2 * (Xmax - Xmin))
        except NameError:    # Xmin does not exist, which is caused by all-zero inducer
            raise Exception("All-zero inducer input is not accepted")
        plotRange = np.exp(np.arange(*plotBounds, 1E-2))
    else:
        plotRange = np.arange(max(np.exp(Xmin) - 0.5, 0.0), np.exp(Xmax) + 0.5, 1E-2)
    if zeroFlag:    # give back the 0 at front
        plotRange = np.insert(plotRange, 0, 0.0)
    mu = model(plotRange.reshape(-1, 1), theta)

    axFlag = True
    if ax == None:
        fig, ax = plt.subplots(figsize=(6, 4))
        axFlag = False
    ax.plot(plotRange, mu, 'g-')
    i= 0
    if 'fmt' not in plotKW:
        plotKW['fmt'] = 'o'
    if 'alpha' not in plotKW:
        plotKW['alpha'] = 0.7
    for vx, vy, vstd in zip(X, Y, std):
        if 'color' not in plotKW:
            plotKW['color'] = palette[i]
        ax.errorbar(vx, vy, yerr=vstd, **plotKW)
        i += 1
    ax.set_xlabel("[%s]/%s"%(inputTag, inputUnits))
    ax.set_ylabel("[%s]/%s"%(outputTag, outputUnits))
    if logScale:
        #ax.set_xscale("symlog", linthresh = np.exp(plotBounds[0] + 1))
        ax.set_xscale("log")
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    # remove the 0 tick label
    xlim = ax.get_xlim()
    ax.set_xticks(list(ax.get_xticks())[1:])
    ax.set_xlim(xlim)
    # set y tick label to scientific notation
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if save_fig is not None:
        # Then save as a file
        plt.savefig(save_fig, dpi=300, format="png")
    elif axFlag == False:
        plt.show()

def plotDNACircuit(modelType, modelSpecs, inducerTags, reporterTag, figPath=None):
    '''
    plotDNACircuit(modelType, modelSpecs, inducerTags, reporterTag, figPath=None)
    Renders the inferred model of the biological node as genetic circuits.
    The plot conforms SBOL visual standards
    '''

    # color scheme
    palette = ['#F79F1F', '#A3CB38', '#1289A7', '#D980FA', '#B53471']
    np.random.shuffle(palette)

    # components
    design = []
    lw = 1.0
    cur = 0
    if modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        inPart = {'type': 'CDS', 'name': 'inPart', 'fwd': True,
                'opts': {'linewidth': lw, 'color': palette[cur], 'edge_color': 'black', 'x_extent': 24,
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
                'opts': {'linewidth': lw, 'color': palette[cur], 'edge_color': 'black', 'x_extent': 24,
                    'label': inducerTags[0], 'label_style': 'italic', 'label_color': 'black', 'label_size': 15,
                    'label_x_offset': -3, 'label_y_offset': 0
                }
        }
        cur += 1
        inPromoter1 = {'type': 'Promoter', 'name': 'pIn1', 'fwd': True,
                'opts': {'linewidth': lw, 'color': 'black', 'label': ''}}
        inPart2 = {'type': 'CDS', 'name': 'inPart2', 'fwd': True,
                'opts': {'linewidth': lw, 'color': palette[cur], 'edge_color': 'black', 'x_extent': 24,
                    'label': inducerTags[1], 'label_style': 'italic', 'label_color': 'black', 'label_size': 15,
                    'label_x_offset': -3, 'label_y_offset': 0
                }
        }
        cur += 1
        inPromoter2 = {'type': 'Promoter', 'name': 'pIn2', 'fwd': True,
                'opts': {'linewidth': lw, 'color': 'black', 'label': ''}}
        design += [inPromoter1, inPart1, inPromoter2, inPart2]
    outPart = {'type': 'CDS', 'name': 'outPart', 'fwd': True,
            'opts': {'linewidth': lw, 'color': palette[cur], 'edge_color': 'black', 'x_extent': 24,
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
