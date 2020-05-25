# -*- coding: utf-8 -*-

'''
Author: Yujia Liu    rainl199922@gmail.com
Xi'an Jiaotong University, Biomedical Informatics & Genomics Center
2020.03.08

modelBase.py is the knowledgebase of bio-models. It includes specifications/definitions such as controlled vocabularies, equations and default bio-model sets, etc. Users are allowed to compose biomodels based on the knowledgebase (and also make extensions to it).
'''

from enum import Enum
from sympy import symbols, lambdify, pretty, diff, Mul
import warnings
import numpy as np

# ModelType includes basic types. Constant: constitutive expression; Inducible or Single_Input_Node: inducible promoter system; Due_Input_Node: promoter w/ 2 TFs, or so called "logic gates"
ModelType = Enum('ModelType', ('Constant', 'Inducible', 'Duo_Input_Node', 'Single_Input_Node'))

# ModelSpec includes model specifications, the combinations of which builds the model
ModelSpec = Enum('ModelSpec', (
    'Linear', 'Michaelis_Menten', 'M_M', 'Quadratic', 'Dimerized', 'Hill',
    'Activation', 'Repression', 'Basal_expression', 'No_basal_expression',
    'Inducer', 'Inducer_Michaelis_Menten', 'Inducer_Quadratic', 'Inducer_Activation', 'Inducer_Repression',
    'AND', 'OR', 'IMPLY12', 'IMPLY21', 'NOR', 'NAND', 'Duo_Activation', 'Duo_Repression'))

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

# Store model equations and jacobians to boost the process
eqnBase = {}
jacBase = {}

def initializeModelBase():
    '''
    initializeModelBase() constructs all functional models in the model-base,
    specifically, those used in genModel() and genJacobian().
    When fitting models in large scale, call this function first to accelerate the process.
    '''
    
    modelSetAll = genModelSet(ModelSet["All"])
    for modelType, modelSpecs in modelSetAll:
        eqn, thetaKeys, _ = genEquation(modelType, modelSpecs)
        key = (modelType, frozenset(modelSpecs))
        eqnBase[key] = eqn
        jac = [diff(eqn, x) for x in thetaKeys]
        jacBase[str(eqn)] = jac

def genEquation(modelType, modelSpecs=()):
    '''
    (equation_func, equation_str), theta_keys, constraints = genEquation(model_type, model_specifications=())
    
    Generate the mathematical expression of the model corresponding to given model-type and model-specifications.
    Other functions, e.g. genModel() wraps the equation and equips it with a parameter handler
    '''

    # TODO: Parameter type inspection

    # Shortcut model equation if is initialized
    key = (modelType, frozenset(modelSpecs))
    if key in eqnBase:
        return eqnBase[key]
    
    # Mathematical expression
    P_st = None
    # names of theta
    thetaKeys = []
    # inequity or equity constraints
    constraints = []

    # Generate model (and corresponding theta list & constraints)
    # Constitutive expression
    if modelType == ModelType.Constant:
        alpha, A = symbols('alpha A')
        # for each model, we define P_st (steady-state protein level), thetaKeys (the name of theta),
        # and constraints (equity and inequity)
        P_st = alpha + Mul(A, 0, evaluate=False)    # A hack to preserve the dimension if input is a vector
        thetaKeys = ['alpha']
        constraint = {
            "type": "ineq",
            "fun": lambda x: x[0]
        }
        constraints.append(constraint)

    # Inducible system (single-input node)
    # For all inducible models, 'A' is the independent variable
    elif modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        # Linear (null) model
        if ModelSpec.Linear in modelSpecs:
            A, alpha, b = symbols('A alpha b')
            P_st = alpha * A + b
            thetaKeys = ["alpha", "b"]
            # no constraints

        # Non-linear models
        else:
            A, alpha, b, K, n = symbols('A alpha b K n')
            P_st = (alpha * A**n + b * K**n) / (A**n + K**n)
            thetaKeys = ['alpha', 'b', 'K', 'n']
            # Michaelian (n=1)
            if ModelSpec.Michaelis_Menten in modelSpecs or ModelSpec.M_M in modelSpecs:
                P_st = P_st.subs(n, 1)
                thetaKeys.remove('n')

            # Dimerization (n=2, cooperative)
            elif ModelSpec.Quadratic in modelSpecs or ModelSpec.Dimerized in modelSpecs:
                P_st = P_st.subs(n, 2)
                thetaKeys.remove('n')

            # Activation
            if ModelSpec.Activation in modelSpecs and ModelSpec.No_basal_expression not in modelSpecs:
                alphaIdx, bIdx = thetaKeys.index('alpha'), thetaKeys.index('b')
                constraint = {
                    "type": "ineq",
                    "fun": lambda x: x[alphaIdx] - x[bIdx]
                }
                constraints.append(constraint)

            # Repression
            elif ModelSpec.Repression in modelSpecs and ModelSpec.No_basal_expression not in modelSpecs:
                alphaIdx, bIdx = thetaKeys.index('alpha'), thetaKeys.index('b')
                constraint = {
                    "type": "ineq",
                    "fun": lambda x: x[bIdx] - x[alphaIdx]
                }
                constraints.append(constraint)

            # No basal expression
            if ModelSpec.No_basal_expression in modelSpecs:
                if ModelSpec.Activation in modelSpecs:
                    P_st = P_st.subs(b, 0)
                    thetaKeys.remove('b')
                elif ModelSpec.Repression in modelSpecs:
                    P_st = P_st.subs(alpha, 0)
                    thetaKeys.remove('alpha')
                else:    # default: activation
                    P_st = P_st.subs(b, 0)
                    thetaKeys.remove('b')

            # Secondary inducer
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
                P_st = P_st.subs(A, A_aster)
                thetaKeys.append('K_I')

            # General constraints for all inducible systems
            for key in ['alpha', 'b', 'K', 'K_I', 'n']:
                if key in thetaKeys:
                    idx = thetaKeys.index(key)
                    constraint = {
                            "type": "ineq",
                            "fun": lambda x: x[idx]
                    }
                    constraints.append(constraint)

    elif modelType == ModelType.Duo_Input_Node:
        # The general 2-input model
        A_1, A_2, K_1, K_2, n_1, n_2, alpha_1, alpha_2, alpha_3, alpha_4 = \
                symbols('A_1 A_2 K_1 K_2 n_1 n_2 alpha_1 alpha_2 alpha_3 alpha_4')
        P_st = (alpha_1 * (A_1/K_1) ** n_1 * (A_2/K_2) ** n_2 + alpha_2 * (A_1/K_1) ** n_1 + alpha_3 * (A_2/K_2) ** n_2 + alpha_4) /\
                ((A_1/K_1) ** n_1 * (A_2/K_2) ** n_2 + (A_1/K_1) ** n_1 + (A_2/K_2) ** n_2 + 1)
        thetaKeys = ['K_1', 'K_2', 'n_1', 'n_2', 'alpha_1', 'alpha_2', 'alpha_3', 'alpha_4']
        # AND gate
        if ModelSpec.AND in modelSpecs:
            P_st = P_st.subs([(alpha_2, 0), (alpha_3, 0), (alpha_4, 0)])
            thetaKeys = ['K_1', 'K_2', 'n_1', 'n_2', 'alpha_1']

        # OR gate
        elif ModelSpec.OR in modelSpecs:
            P_st = P_st.subs(alpha_4, 0)
            thetaKeys.remove('alpha_4')

        # 1 IMPLY 2 gate
        elif ModelSpec.IMPLY12 in modelSpecs:
            P_st = P_st.subs(alpha_3, 0)
            thetaKeys.remove('alpha_3')

        # 2 IMPLY 1 gate
        elif ModelSpec.IMPLY21 in modelSpecs:
            P_st = P_st.subs(alpha_2, 0)
            thetaKeys.remove('alpha_2')

        # NOR gate
        elif ModelSpec.NOR in modelSpecs:
            P_st = P_st.subs([(alpha_1, 0), (alpha_2, 0), (alpha_3, 0)])
            thetaKeys = ['K_1', 'K_2', 'n_1', 'n_2', 'alpha_4']

        # NAND gate
        elif ModelSpec.NAND in modelSpecs:
            P_st = P_st.subs(alpha_1, 0)
            thetaKeys.remove('alpha_1')

        # General activation
        if ModelSpec.Duo_Activation in modelSpecs:
            i1, i2, i3, i4 = thetaKeys.index('alpha_1'), thetaKeys.index('alpha_2'),\
                    thetaKeys.index('alpha_3'), thetaKeys.index('alpha_4')
            con1 = {
                    "type": "ineq",
                    "fun": lambda x: x[i1] - x[i2]
            }
            con2 = {
                    "type": "ineq",
                    "fun": lambda x: x[i1] - x[i3]
            }
            con3 = {
                    "type": "ineq",
                    "fun": lambda x: x[i2] - x[i4]
            }
            con4 = {
                    "type": "ineq",
                    "fun": lambda x: x[i3] - x[i4]
            }
            constraints += [con1, con2, con3, con4]
        
        # General Repression
        elif ModelSpec.Duo_Repression in modelSpecs:
            i1, i2, i3, i4 = thetaKeys.index('alpha_1'), thetaKeys.index('alpha_2'),\
                    thetaKeys.index('alpha_3'), thetaKeys.index('alpha_4')
            con1 = {
                    "type": "ineq",
                    "fun": lambda x: x[i4] - x[i2]
            }
            con2 = {
                    "type": "ineq",
                    "fun": lambda x: x[i4] - x[i3]
            }
            con3 = {
                    "type": "ineq",
                    "fun": lambda x: x[i2] - x[i1]
            }
            con4 = {
                    "type": "ineq",
                    "fun": lambda x: x[i3] - x[i1]
            }
            constraints += [con1, con2, con3, con4]

        # General constraints
        # TODO: do alpha_i have to be non-negative?
        for key in ['K_1', 'K_2', 'n_1', 'n_2']:
            if key in thetaKeys:
                idx = thetaKeys.index(key)
                constraint = {
                        "type": "ineq",
                        "fun": lambda x: x[idx]
                }
                constraints.append(constraint)

    eqnFunc, eqnStr = __lambdifyEqn(P_st, thetaKeys, modelType)
    return (eqnFunc, eqnStr, P_st), thetaKeys, constraints

def __lambdifyEqn(eqn, thetaKeys, modelType):
    '''
    Lambdify symbolic expression of model. ModelType is required to determine the dimension of inputs
    '''

    # Pretty print expressions
    eqnStr = pretty(eqn, use_unicode=False)

    # lambdify equaition
    A, A_1, A_2 = symbols('A A_1 A_2')
    if modelType == ModelType.Constant or\
            modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        eqnFunc = lambdify((A, symbols(thetaKeys)), eqn, 'numpy')
    elif modelType == ModelType.Duo_Input_Node:
        eqnFunc = lambdify(((A_1, A_2), symbols(thetaKeys)), eqn, 'numpy')
    
    return eqnFunc, eqnStr

def genModel(modelType, modelSpecs = ()) :
    '''
    model, constraints, thetaKeys, eqnStr, eqnSym = genModel(modelType, modelSpecs)
    '''

    # Check variable types. modelType: enum, modelSpecs: list/tuple
    if type(modelSpecs) != list and type(modelSpecs) != tuple :
        raise Exception('Usage: genModel(modelType, modelSpecs), modelSpecs should be either list or tuple')
    for key in modelSpecs:
        if type(key) != ModelSpec:
            raise Exception('Usage: genModel(modelType, modelSpecs), modelSepcs should be a list/tuple of ModelSpec enums')
    if type(modelType) != ModelType :
        raise Exception('Usage: genModel(modelType, modelSpecs), modelType should be a ModelType Enum')

    # Check validation and coexistence of modelType & modelSpecs
    __modelSpecsInspector(modelType, modelSpecs)

    # Generate the mathematical expression of models
    (eqnFunc, eqnStr, eqnSym), thetaKeys, constraints = genEquation(modelType, modelSpecs)

    # Parameter injection
    model = __modelParaInterpreter(modelType, eqnFunc, thetaKeys)

    return model, constraints, thetaKeys, eqnStr, eqnSym

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

defaultTheta = {
        'alpha': 1.0, 'alpha_1': 1.0, 'alpha_2': 0.2, 'alpha_3': 0.2, 'alpha_4': 0.1,
        'beta': 1.0, 'b': 0.1,
        'K': 1E-3, 'K_I': 1E-3, 'K_1': 1E-3, 'K_2': 1E-3,
        'n': 1.0, 'n_1': 1.0, 'n_2': 1.0
}

# Supplement missing parameters
def __supplementPara(theta, thetaKeys, quiet=False):
    for key in thetaKeys:
        if key not in theta.keys():
            theta[key] = defaultTheta[key]
            if not quiet:
                warnings.warn('The parameter %s is not set, set to default value %f'%(key, theta[key]))
    # this func directly modifies the original theta

def __modelParaInterpreter(modelType, modelFunc, thetaKeys):
    def model(X, theta={}):
        # Check input
        if modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
            numInput = 1
        elif modelType == ModelType.Duo_Input_Node:
            numInput = 2
        elif modelType == ModelType.Constant:
            numInput = 0
        try:
            inputDim = len(X[0])
        except TypeError:
            raise Exception('Usage: model(X, theta), where X should be k * N 2D matrix, k is input dimension, N is sample size')
        if numInput > 0 and numInput != inputDim:
            raise Exception('Usage: model(X, theta), where X should be k * N 2D matrix, k is input dimension, N is sample size')

        # supplement missing parameters, if any
        __supplementPara(theta, thetaKeys)
        serializedTheta = [theta[key] for key in thetaKeys]

        X = np.array(X, dtype='float64')    # convert input to np array
        reporter = np.array([modelFunc(A, serializedTheta) for A in X])

        return np.squeeze(reporter)    # eliminate unwanted dimensions

    return model

def __modelSpecsInspector(modelType, modelSpecs):
    '''
    __modelSpecsInspector() checks the validation of given model-type and model-specifications.
    Some model-specifications and/or model-types cannot coexist.
    '''

    if modelType == ModelType.Constant :
        if len(modelSpecs) != 0:
            raise Exception('Usage: genModel(modelType, modelSpecs), where if modelType is "Constant", no modelSpecs should be specified')


    elif modelType == ModelType.Inducible or modelType == ModelType.Single_Input_Node:
        # 'Linear' cannot co-exist w/ any except 'Activation' or 'Repression
        if ModelSpec.Linear in modelSpecs:
            if len(modelSpecs) > 2:
                raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Linear" cannot co-exist with any other keywords other than "Repression" or "Activation')
            elif len(modelSpecs) == 2:
                tmpModelSpecs = list(modelSpecs).remove(ModelSpec.Linear)
                if not (tmpModelSpecs[0] == ModelSpec.Activation or\
                        tmpModelSpecs[0] == ModelSpec.Repression):
                    raise Exception('Usage: genModel(modelType, modelSpecs), in modelSpecs, "Linear" cannot co-exist with any other keywords other than "Repression" or "Activation')

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

    elif modelType == ModelType.Duo_Input_Node:
        # Inducible/Constant keyword are not compatible with duo input
        duoInputKeys = [ModelSpec.AND, ModelSpec.OR, ModelSpec.IMPLY12,
                ModelSpec.IMPLY21, ModelSpec.NOR, ModelSpec.NAND,
                ModelSpec.Duo_Activation, ModelSpec.Duo_Repression]
        for keyword in modelSpecs:
            if keyword not in duoInputKeys:
                raise Exception('Usage: genModel(modelType, modelSpecs), \
                        where model-spec "%s" is not compatible with model-type "Duo_Input_Node"'%keyword.name)
        # For duo-input-allowed specifications, non can co-exist
        if len(modelSpecs) > 1:
            raise Exception('Usage: genModel(modelType, modelSpecs), where model-spec "%s" and "%s" cannot coexist'\
                    %(modelSpecs[0].name, modelSpecs[1].name))

def genModelSet(modelSet):
    '''
    The func gives the common sets of mechanistic models
    Each model is a pair of (modelType, ModelSpecs)
    '''

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
