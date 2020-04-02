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

    X1D = [x[0] for x in X]    # Extract the 1st dimension
    Xmin, Xmax = min(X1D), max(X1D)
    # 10% extrapolation on both ends
    plotBounds = (Xmin - 0.1 * (Xmax - Xmin), Xmax + 0.1 * (Xmax - Xmin))
    plotRange = np.arange(*plotBounds, 1E-5)
    Y = model(plotRange, theta)

    fig, ax = plt.subplots(figsize=(8, 6))
    #ax.plot(X1D, mu, color="blue", marker='^', linestyle=None)
    ax.plot(plotRange, Y, 'k-')
    ax.errorbar(X1D, mu, yerr=std, fmt='o', color="blue")
    ax.set_xlabel("[%s]/%s"%(inputTag, inputUnits))
    ax.set_ylabel("[%s]/%s"%(outputTag, outputUnits))
    ax.set_xscale("log")
    if save_fig is not None:
        # Then save as a file
        plt.savefig(save_fig, dpi=300, format="png")
    else:
        plt.show()
