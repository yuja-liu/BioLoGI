# /usr/bin/env python

import pandas as pd

def dataframe_check(mat):
    if type(mat) != pd.core.frame.DataFrame:
        raise Exception("argument mat must be a pandas DataFrame")
        
def normalize(mat):
    dataframe_check(mat)
    mat = mat.sub(mat.mean(axis = 1), axis = 0)
    mat = mat.div((mat * mat).sum(axis = 1)**(1 / 2), axis = 0)
    return mat

def pearson_corr_mat(mat):
    dataframe_check(mat)
    mat_norm = normalize(mat)
    return mat_norm.dot(mat_norm.transpose())

def spearman_corr_mat(mat):
    dataframe_check(mat)
    mat_rank = mat.rank(axis = 1)
    return pearson_corr_mat(mat_rank)