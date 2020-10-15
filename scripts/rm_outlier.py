import numpy as np
from sklearn.neighbors import LocalOutlierFactor

def remove_outliers_2d(X, Y):
    X, Y = np.array(X), np.array(Y)
    D = np.hstack([X.reshape(-1, 1), Y.reshape(-1, 1)])
    lof = LocalOutlierFactor(n_neighbors = 100, contamination = 0.1)
    is_outlier = lof.fit_predict(D) == -1
    new_X, new_Y = X[~is_outlier], Y[~is_outlier]
    return new_X, new_Y

def remove_outliers_3d(X, Y, Z):
    X, Y, Z = np.array(X), np.array(Y), np.array(Z)
    D = np.hstack([X.reshape(-1, 1), Y.reshape(-1, 1), Z.reshape(-1, 1)])
    lof = LocalOutlierFactor(n_neighbors = 100, contamination = 0.1)
    is_outlier = lof.fit_predict(D) == -1
    new_X, new_Y, new_Z = X[~is_outlier], Y[~is_outlier], Z[~is_outlier]
    return new_X, new_Y, new_Z