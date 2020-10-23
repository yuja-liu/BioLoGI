#!/usr/bin/env python

'''
BioLoGI AND gate/dual-input gate test
Using GEO GSE75748 h9 hESC data
picked RNA:
    input: ENSG00000235910.1 (APOA1-AS), ENSG00000141431.9 (ASXL3)
    output: ENSG00000136883.12 (KIF12)

Author: Yujia Liu rainl199922@gmail.com
2020.10.14
'''

import numpy as np
import pandas as pd
import time
import sys
sys.path.append("../")
from biologiclib.modelBase import *
from biologiclib.plotUtils import *
from biologiclib.inference import *
sys.path.append("../scripts")
from rm_outlier import *

# load expression matrix
h9_exp_mat = pd.read_csv("../data/GSE75748/h9_imputed.tsv",
                         sep = "\t", index_col = 0)

# define input and output
id_inducer = ("ENSG00000235910.1", "ENSG00000141431.9")
id_reporter = "ENSG00000136883.12"
symbol_inducer = ("APOA1-AS", "ASXL3")
symbol_reporter = "KIF12"
count_x = [[], []]
count_x[0], count_x[1] =\
        h9_exp_mat.loc[id_inducer[0]].to_numpy(),\
        h9_exp_mat.loc[id_inducer[1]].to_numpy()
count_y = h9_exp_mat.loc[id_reporter].to_numpy()
count_x[0], count_x[1], count_y =\
        remove_outliers_3d(*count_x, count_y)

# fitting
count_x = np.array(count_x).transpose()    # convert to column vectors
start_time = time.time()
best_model, all_models = selectModel(count_x,
                            count_y,
                            modelSolver = ModelSolver.SLSQP,
                            modelSet = ModelSet.Dual_Activation,
                            parallel = False,
                            quiet = False)
elapsed = time.time() - start_time
best_model_specs = [spec.name for spec in best_model.modelSpecs]
print("Best model keywords:", best_model_specs)
print("Best model AIC:", best_model.IC)
print("Best model parameters:",
      {key: val for key, val in zip(best_model.thetaKey, best_model.thetaVal)})
print("Time elapse:", elapsed)

# plot
plotHelper3D(count_x,
        count_y,
        best_model,
        inducer_name = symbol_inducer,
        reporter_name = symbol_reporter)
