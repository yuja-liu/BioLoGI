# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.08

modelBase.py is the knowledgebase of bio-models. It includes specifications/definitions such as controlled vocabularies, equations and default bio-model sets, etc. Users are allowed to compose biomodels based on the knowledgebase (and also make extensions to it).
'''

from enum import Enum
from sympy import symbols, lambdify, pretty
import warnings
import numpy as np

# ModelType includes basic types. Constant: constitutive expression; Inducible or Single_Input_Node: inducible promoter system; Due_Input_Node: promoter w/ 2 TFs, or so called "logic gates"
ModelType = Enum('ModelType', ('Constant', 'Inducible', 'Duo_Input_Node', 'Single_Input_Node'))

# ModelSpec includes model specifications, the combinations of which builds the model
ModelSpec = Enum('ModelSpec', ('Linear', 'Michaelis_Menten', 'M_M', 'Quadratic', 'Dimerized', 'Hill', 'Activation', 'Repression', 'Basal_expression', 'No_basal_expression',
    'Inducer', 'Inducer_Michaelis_Menten', 'Inducer_Quadratic', 'Inducer_Activation', 'Inducer_Repression'))

# model set constants
ModelSet = Enum("ModelSet", ("Simple_Inducible_Promoter", "Inducible_Promoter_with_Inducer", "Activation_System", "Repression_System", "All", "Minimum"))

# CVSet is the set of all controlled vocabulary
CVSet = {ModelType, ModelSpec, ModelSet}

def printCV() :
    for CVGroup in CVSet :
        print('# ' + CVGroup.__name__)
        nameList = [item.name for item in CVGroup]
        print(', '.join(nameList))

def __findCV(cv) :
    for CVGroup in CVSet :
        if cv in [item.name for item in CVGroup] :
            return CVGroup.__name__
    return None

def whichCV(cv) :
    tmp = __findCV(cv)
    print(tmp if tmp is not None else 'CV not found')

def isCV(cv) :
    print(__findCV(cv) is not None)

def genModel(modelType, modelSpecs = (), plain_print=True) :
    ########################################
    # Check variable types. modelType: enum, modelSpecs: list/tuple
    ########################################
    if type(modelSpecs) != list and type(modelSpecs) != tuple :
        raise Exception('Usage: genModel(modelType, modelSpecs), modelSpecs should be either list or tuple')
    for key in modelSpecs:
        if type(key) != ModelSpec:
            raise Exception('Usage: genModel(modelType, modelSpecs), modelSepcs should be a list/tuple of ModelSpec enums')
    if type(modelType) != ModelType :
        raise Exception('Usage: genModel(modelType, modelSpecs), modelType should be a ModelType Enum')

    ########################################
    # modelType = Constant(constitutive expression)
    ########################################
    if modelType == ModelType.Constant :
        if len(modelSpecs) != 0:
            raise Exception('Usage: genModel(modelType, modelSpecs), where if modelType is "Constant", no modelSpecs should be specified')
        alpha = symbols('alpha')
        Pst = alpha    # Pst is the steady-state [P]

    ########################################
    # modelType = Inducible(single-input node, i.e. inducible promoter system)
    ########################################
    if modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        # Some modelSpecs cannot co-exist
        # "Activation" cannot co-exist w/ "Repression"
        if ModelSpec.Activation in modelSpecs and ModelSpec.Repression in modelSpecs :
            raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Activation" and "Repression" cannot coexist')
        if ModelSpec.Inducer_Activation in modelSpecs and ModelSpec.Inducer_Repression in modelSpecs:
            raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Inducer_Activation" and "Inducer_Repression" cannot coexist')
        # "Michaelis_Menten" cannot co-exist w/ "Quadratic" or "Hill", though "Hill" does not specify Hill coefficient
        if int(ModelSpec.Michaelis_Menten in modelSpecs or ModelSpec.M_M in modelSpecs) + \
                int(ModelSpec.Quadratic in modelSpecs or ModelSpec.Dimerized in modelSpecs) + int(ModelSpec.Hill in modelSpecs) >= 2 :
            raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Michaelis_Menten" (or "M_M"), "Quadratic" (or "Dimerized"), or "Hill" cannot coexist')
        if ModelSpec.Inducer_Michaelis_Menten in modelSpecs and ModelSpec.Inducer_Quadratic in modelSpecs:
            # Inducer "Hill coefficient" is constrained to 1 or 2
            raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Inducer_Michaelis_Menten" and "Inducer_Quadratic" cannot coexist')
        # "Basal_expression" cannot co-exit w/ "No_basal_expression"
        if ModelSpec.Basal_expression in modelSpecs and ModelSpec.No_basal_expression in modelSpecs :
            raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Basal_expression" and "No_basal_expression" cannot coexist')
        # "Inducer_+" keyword depends on "Inducer"
        for keyword in (ModelSpec.Inducer_Activation, ModelSpec.Inducer_Repression, ModelSpec.Inducer_Michaelis_Menten, ModelSpec.Inducer_Quadratic):
            if keyword in modelSpecs and ModelSpec.Inducer not in modelSpecs:
                raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "%s" depends on keyword "Inducer"'%keyword.name)

        ########################################
        # Hill Function
        ########################################
        # A, activator or repressor; alpha, expression(txn + tln) rate; beta, degradation rate, now obsolete; b, leaky expression rate; K, dissociation constant; n, Hill coefficient
        A, alpha, b, K, n = symbols('A alpha b K n')

        Pst = (alpha * A**n + b * K**n) / (A**n + K**n)
        if ModelSpec.Michaelis_Menten in modelSpecs or ModelSpec.M_M in modelSpecs:
            Pst = Pst.subs(n, 1)
        elif ModelSpec.Quadratic in modelSpecs or ModelSpec.Dimerized in modelSpecs:
            Pst = Pst.subs(n, 2)
        if ModelSpec.No_basal_expression in modelSpecs:
            if ModelSpec.Activation in modelSpecs:
                Pst = Pst.subs(b, 0)
            elif ModelSpec.Repression in modelSpecs:
                Pst = Pst.subs(alpha, 0)
            else:
                Pst = Pst.subs(b, 0)

        ########################################
        # Negative control: linear model
        ########################################
        if ModelSpec.Linear in modelSpecs:
            Pst = alpha * A + b

        ########################################
        # "Secondary inducer"
        ########################################
        if ModelSpec.Inducer in modelSpecs:
            n_I, K_I, a_I, b_I = symbols('n_I K_I a_I b_I')
            # General model for inducer binding
            A_aster = ((a_I*A**n_I + b_I*K_I**n_I) / (A**n_I + K_I**n_I))
            if ModelSpec.Inducer_Activation in modelSpecs: 
                A_aster = A_aster.subs([(a_I, 1), (b_I, 0)])
            elif ModelSpec.Inducer_Repression in modelSpecs:
                A_aster = A_aster.subs([(a_I, 0), (b_I, 1)])
            else:   #default
                A_aster = A_aster.subs([(a_I, 1), (b_I, 0)])
            if ModelSpec.Inducer_Michaelis_Menten in modelSpecs:
                A_aster = A_aster.subs(n_I, 1)
            elif ModelSpec.Inducer_Quadratic in modelSpecs:
                A_aster = A_aster.subs(n_I, 2)
            else:    #default
                A_aster = A_aster.subs(n_I, 1)
            Pst = Pst.subs(A, A_aster)

    ########################################
    # Pretty print expressions
    ########################################
    if plain_print:
        returnExp = pretty(Pst, use_unicode=False)
    else:
        returnExp = Pst

    ########################################
    # Parameter injection
    ########################################
    model, constraints, paraList = __modelParaInterpreter(modelType, modelSpecs, Pst)
    return returnExp, model, constraints, paraList

# Generate default parameters: useful for initials
def defaultPara(thetaList, inducer, reporter, repression = False):
    default = {}
    __supplementPara(default, thetaList, True)

    # get default parameter according to input data
    # TODO: the relationship between alpha & b fits only activation!!
    if 'alpha' in thetaList:    # Actually, alpha should always in theta
        if 'b' in thetaList and not repression:
            default['alpha'] = np.median(reporter)
            default['b'] = default['alpha'] / 10
        elif 'b' in thetaList and repression:
            default['b'] = np.median(reporter)
            default['alpha'] = default['b'] / 10
    else:    # repression
        default['b'] = np.median(reporter)
    if 'I' in thetaList:
        default['K_I'] = np.median(inducer)
    else:
        default['K'] = np.median(inducer)
    theta = [default[key] for key in thetaList]
    return theta

# Supplement missing parameters
def __supplementPara(theta, thetaList, quiet=False):
    defaultTheta = {'alpha': 1.0, 'beta': 1.0, 'b': 0.1, 'K': 1E-3, 'n': 1.0, 'K_I': 1E-3}

    if type(theta) != dict :
        raise Exception('Usage: model(x, theta), where theta should be a dict including e.g. "alpha", the transcription & translation strength')
    for key in thetaList:
        if key not in theta.keys():
            theta[key] = defaultTheta[key]
            if not quiet:
                warnings.warn('The parameter %s is not set, set to default value %f'%(key, theta[key]))

def __modelParaInterpreter(modelType, modelSpecs, modelExp) :
    # modelType = Constant
    if modelType == ModelType.Constant:
        thetaList = ['alpha']
        def model(X, theta = {}):
            __supplementPara(theta, thetaList)
            val = modelExp.evalf(subs=theta)
            try:
                return [val] * len(X)
            except ValueError:
                return val

    # modelType = Inducible
    elif modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        # Generate the required theta list
        if ModelSpec.Linear in modelSpecs:
            thetaList = ["alpha", "b"]
        else:
            thetaList = ['alpha', 'b', 'K', 'n']
        for key in (ModelSpec.Michaelis_Menten, ModelSpec.M_M, ModelSpec.Quadratic, ModelSpec.Dimerized):
            if key in modelSpecs:
                thetaList.remove('n')
        if ModelSpec.No_basal_expression in modelSpecs:
            if ModelSpec.Activation in modelSpecs:
                thetaList.remove('b')
            elif ModelSpec.Repression in modelSpecs:
                thetaList.remove('alpha')
            else:    #default
                thetaList.remove('b')
        if ModelSpec.Inducer in modelSpecs:
            thetaList.append('K_I')    # If "Inducer", activator/repressor concentration is considered as "hidden", i.e. parameter

        # Returning reporter = f(inducer; theta)
        def model(X, theta = {}):

            # theta must be dict. If not, errors would be raised inside __supplementPara
            __supplementPara(theta, thetaList)

            # X is obligated to be iterable
            try:
                len(X)
            except ValueError:
                raise Exception("Usage: model(X, theta)), where X should be iterable, e.g., a list. Furthermore, X should be a N*k matrix, where k is the dimension of inducer, N is the sample size.")
            
            # "Invalid value walls"
#            #########################################
#            # Check for negative parameters
#            for key, val in theta.items():
#                if val < 0:
#                    #warnings.warn("The parameter %s should be non-negative, but here is set to %f. Returning INF"%(key, val))
#                    return X * 1E5
#
#            # Activation/Repression is controled by parameter conditions
#            if ModelSpec.No_basal_expression not in modelSpecs:
#                # If there is no basal expression, than the system is purely activation or repression
#                if ModelSpec.Activation in modelSpecs and theta['alpha'] < theta['b']:
#                    warnings.warn('"Activation" keyword implies parameter alpha >= b, which is not the case. Returning INF')
#                    return X * 1E5    # a hack for vector input
#                elif ModelSpec.Repression in modelSpecs and theta['alpha'] > theta['b']:
#                    warnings.warn('"Repression" keyword implies parameter alpha <= b, which is not the case. Returning INF')
#                    return X * 1E5
#            #########################################

            # Evaluation
            newModelExp = modelExp.subs(theta)
            modelFunc = lambdify(symbols('A'), newModelExp, 'numpy')
            Y = []
            for x in X:
                # inducer is 1 dimensional
                try:
                    if len(x) != 1:
                        raise Exception("Each inducer should be 1 dimensional vector")
                    x = x[0]
                except:
                    pass
                Y.append(modelFunc(x))
            return Y

    ########################################
    # Generate constraints
    ########################################
    constraints = []
    nonNegPara = ['alpha', 'b', 'K', 'n', 'K_I']
    for para in nonNegPara:
        if para in thetaList:
            idx = thetaList.index(para)
            constraint = {
                "type": "ineq",
                "fun": lambda x: x[idx]
            }
            constraints.append(constraint)
    # for Activation or Repression
    if ModelSpec.Activation in modelSpecs and ModelSpec.No_basal_expression not in modelSpecs:
        constraint = {
            "type": "ineq",
            "fun": lambda x: x[thetaList.index('alpha')] - x[thetaList.index('b')]
        }
        constraints.append(constraint)
    elif ModelSpec.Repression in modelSpecs and ModelSpec.No_basal_expression not in modelSpecs:
        constraint = {
            "type": "ineq",
            "fun": lambda x: x[thetaList.index('b')] - x[thetaList.index('alpha')]
        }
        constraints.append(constraint)

    ########################################
    # Return
    ########################################
    return model, constraints, thetaList

# The func gives the common sets of mechanistic models
# Each model is a pair of (modelType, ModelSpecs)
def genModelSet(modelSet):
    # Check parameter
    if type(modelSet) != ModelSet:
        raise Exception("Usage: genModelSet(modelSet), where modelSet should be a ModelSet Enum")
    
    models = []
    if modelSet == ModelSet.All:
        # List all feasible models, a total of 61
        models.append((ModelType.Constant, ()))
        models.append((ModelType.Inducible, (ModelSpec.Linear,)))
        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic, ModelSpec.Hill):
            for j in (ModelSpec.Activation, ModelSpec.Repression):
                for k in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                    models.append((ModelType.Inducible, (i, j, k)))
                    for p in (ModelSpec.Inducer_Michaelis_Menten, ModelSpec.Inducer_Quadratic):
                        for q in (ModelSpec.Inducer_Activation, ModelSpec.Inducer_Repression):
                            models.append((ModelType.Inducible, (i, j, k, ModelSpec.Inducer, p, q)))

    elif modelSet == ModelSet.Simple_Inducible_Promoter:
        # Netative control
        models.append((ModelType.Constant, ()))
        models.append((ModelType.Inducible, (ModelSpec.Linear,)))

        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic, ModelSpec.Hill):
            for j in (ModelSpec.Activation, ModelSpec.Repression):
                for k in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                    models.append((ModelType.Inducible, (i, j, k)))

    elif modelSet == ModelSet.Inducible_Promoter_with_Inducer:
        # Netative control
        models.append((ModelType.Constant, ()))
        models.append((ModelType.Inducible, (ModelSpec.Linear,)))

        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic, ModelSpec.Hill):
            for j in (ModelSpec.Activation, ModelSpec.Reperssion):
                for k in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                    for p in (ModelSpec.Inducer_Michaelis_Menten, ModelSpec.Inducer_Quadratic):
                        for q in (ModelSpec.Inducer_Activation, ModelSpec.Inducer_Repression):
                            models.append((ModelType.Inducible, (i, j, k, ModelSpec.Inducer, p, q)))

    elif modelSet == ModelSet.Activation_System:
        # Netative control
        models.append((ModelType.Constant, ()))
        models.append((ModelType.Inducible, (ModelSpec.Linear,)))
        # Simple Activation
        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic, ModelSpec.Hill):
            for j in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                models.append((ModelType.Inducible, (i, j, ModelSpec.Activation)))
        # Activation with inducer
        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic):
            for j in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                for k in (ModelSpec.Inducer_Michaelis_Menten,):    # For simplification, only keep Michaelis Menten for inducer Hill
                    for p in ((ModelSpec.Activation, ModelSpec.Inducer_Activation),
                            (ModelSpec.Repression, ModelSpec.Inducer_Repression)):
                        models.append((ModelType.Inducible, (i, j, k, *p, ModelSpec.Inducer)))
    elif modelSet == ModelSet.Repression_System:
        # Netative control
        models.append((ModelType.Constant, ()))
        models.append((ModelType.Inducible, (ModelSpec.Linear,)))
        # Simple repression
        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic, ModelSpec.Hill):
            for j in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                models.append((ModelType.Inducible, (i, j, ModelSpec.Repression)))
        # Repression with inducer
        for i in (ModelSpec.Michaelis_Menten, ModelSpec.Quadratic):
            for j in (ModelSpec.Basal_expression, ModelSpec.No_basal_expression):
                for k in (ModelSpec.Inducer_Michaelis_Menten,):
                    for p in ((ModelSpec.Activation, ModelSpec.Inducer_Repression),
                            (ModelSpec.Repression, ModelSpec.Inducer_Activation)):
                        models.append((ModelType.Inducible, (i, j, k, *p, ModelSpec.Inducer)))

    # A small set for testing
    elif modelSet == ModelSet.Minimum:
        models += [(ModelType.Constant, ()),
                (ModelType.Inducible, (ModelSpec.Linear,)),
                (ModelType.Inducible, (ModelSpec.Michaelis_Menten, ModelSpec.Activation, ModelSpec.Basal_expression)),
                (ModelType.Inducible, (ModelSpec.Michaelis_Menten, ModelSpec.Repression, ModelSpec.Basal_expression)),
                (ModelType.Inducible, (ModelSpec.Quadratic, ModelSpec.Activation, ModelSpec.Basal_expression)),
                (ModelType.Inducible, (ModelSpec.Quadratic, ModelSpec.Repression, ModelSpec.Basal_expression)),
                (ModelType.Inducible, (ModelSpec.Hill, ModelSpec.Activation, ModelSpec.Basal_expression)),
                (ModelType.Inducible, (ModelSpec.Hill, ModelSpec.Repression, ModelSpec.Basal_expression))]

    return models
