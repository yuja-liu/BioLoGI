# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.04.02

Ploting utilities
'''

from matplotlib import pyplot as plt
import numpy as np

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
