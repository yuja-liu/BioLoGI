#!/usr/bin/env python

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import time
from sklearn.neighbors import LocalOutlierFactor
from scipy.stats import pearsonr
import json
from multiprocessing import Pool, cpu_count, Array, RawArray

import sys
sys.path.append("../")
from biologiclib.modelBase import *
from biologiclib.plotUtils import *
from biologiclib.inference import *
from rm_outlier import *

def extract_results(best_model, all_models):
    specs = [spec.name for spec in best_model.modelSpecs]
    AIC = best_model.IC
    for model in all_models:
        if ModelSpec.Linear in model.modelSpecs:
            linear_model = model
    try:
        Delta_AIC = best_model.IC - linear_model.IC
    except NameError:
        raise Exception('Linear model not found. Abort')
    theta = {key: val for key, val in zip(best_model.thetaKey, best_model.thetaVal)}
    return specs, AIC, Delta_AIC, theta

# load expression matrix
h9_tpm = pd.read_csv("../data/GSE75748/h9_imputed.tsv", sep = '\t')
shared_h9_tpm = Array('d', h9_tpm.shape[0] * h9_tpm.shape[1])
shared_h9_tpm = np.frombuffer(shared_h9_tpm.get_obj())
shared_h9_tpm[:] = h9_tpm.to_numpy().ravel()
shared_h9_tpm = shared_h9_tpm.reshape(h9_tpm.shape)

def func(idx):
    count_inducer_all = np.copy(shared_h9_tpm[idx[0]])
    count_reporter_all = np.copy(shared_h9_tpm[idx[1]])
    # remove outliers
    count_inducer, count_reporter = remove_outliers_2d(count_inducer_all, count_reporter_all)
    # positive or negative correlation
    r, _ = pearsonr(count_inducer, count_reporter)
    if r > 0:
        model_set = ModelSet.Activation_System
    else:
        model_set = ModelSet.Repression_System
    # restore exponential
    count_inducer, count_reporter = np.exp(count_inducer) - 1, np.exp(count_reporter) - 1

    return_models = []
    # fitting
    best_model_1, all_models_1 = selectModel(count_inducer.reshape(-1, 1), count_reporter,
                            modelSolver = ModelSolver.SLSQP,
                            modelSet = model_set,
                            parallel = False)
    if ModelSpec.Linear not in best_model_1.modelSpecs\
    and len(best_model_1.modelSpecs) != 0:    # then non-linear model
        return_models.append(tuple([idx[0], idx[1]] +\
                              list(extract_results(best_model_1, all_models_1))))
    # fitting, switch inducer and reporter
    best_model_2, all_models_2 = selectModel(count_reporter.reshape(-1, 1), count_inducer,
                            modelSolver = ModelSolver.SLSQP,
                            modelSet = ModelSet.Activation_System,
                            parallel = False)
    if ModelSpec.Linear not in best_model_2.modelSpecs\
    and len(best_model_2.modelSpecs) != 0:    # then non-linear model
        return_models.append(tuple([idx[0], idx[1]] +\
                              list(extract_results(best_model_2, all_models_2))))
    return return_models

if __name__ == "__main__":
    # load filtered gene pairs
    filtered_indices = pd.read_csv("../data/GSE75748/filtered_indices.tsv", sep = '\t')
    NUM = 10
    PROCESSES = cpu_count()
    
    with Pool(PROCESSES) as pool:
        start_time = time.time()
        # Do NOT pass a series into apply_async(), the index will go wrong
        async_results = [pool.apply_async(func, (row,)) for row in filtered_indices.head(NUM).to_numpy()]
        detected_logi_raw = [res.get(timeout = 60) for res in async_results]
        elapsed = time.time() - start_time
        print("time elapse:", elapsed)
    
    detected_logi = []
    for res in detected_logi_raw:
        detected_logi += res
    df_detected_logi = pd.DataFrame(detected_logi,
            columns = ('name_inducer', 'name_reporter', 'specs',
                'AIC', 'Delta_AIC', 'theta'))
    # substitute gene names, cat keywords
    df_detected_logi[['name_inducer', 'name_reporter']] =\
            df_detected_logi[['name_inducer', 'name_reporter']].applymap(
            lambda i: h9_tpm.index.values[i])
    df_detected_logi['specs'] = df_detected_logi['specs'].map(
            lambda s: ', '.join(s))
    print(df_detected_logi.shape[0], "logic gates detected")
    df_detected_logi.to_csv("../data/GSE75748/detected_logi.tsv", sep = "\t", index = False)
