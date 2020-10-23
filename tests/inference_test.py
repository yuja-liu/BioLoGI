import unittest

import sys
sys.path.append('..')

from biologiclib import inference, ioUtils, modelBase
from biologiclib.modelBase import ModelType, ModelSpec, ModelSet
from biologiclib.inference import ModelSolver
import numpy as np
import time

class ModelMethods(unittest.TestCase):

    def test_readInput(self):
        inducer, reporter, reporterStd = ioUtils.readInput("../data/BBa_K2791004.tsv", ["Psicose"], ["Exp1", "Exp2"], "GFP")
        correct_inducer = [
                [0.00001], [0.00010], [0.00050], [0.00100], [0.00200], [0.00500], [0.01000], [0.05000],
                [0.00001], [0.00010], [0.00050], [0.00100], [0.00200], [0.00500], [0.01000], [0.05000]]
        correct_reporter = [
                1.243524, 1.243080, 1.504122, 1.144187, 1.307069, 1.944845, 1.624899, 3.257670,
                1.422357, 1.382356, 1.737699, 1.259750, 1.412099, 1.998595, 1.821425, 3.213432]
        correct_reporterStd = [
                0.043419, 0.055680, 0.061010, 0.037402, 0.023440, 0.024382, 0.082444, 0.051219,
                0.046885, 0.040743, 0.070671, 0.034629, 0.027485, 0.035141, 0.059574, 0.057311]
        for first, second in zip(inducer, correct_inducer):
            self.assertAlmostEqual(first[0], second[0], places=6)
        for first, second in zip(reporter, correct_reporter):
            self.assertAlmostEqual(first, second, places=6)
        for first, second in zip(reporterStd, correct_reporterStd):
            self.assertAlmostEqual(first, second, places=6)

    def test_genModel(self):

        # generate model
        specs = [ModelSpec.Activation, ModelSpec.Hill, ModelSpec.Basal_expression]
        expression, model, constraints, paraList = modelBase.genModel(ModelType.Inducible, specs, plain_print=False)
        correct_paraList = ('alpha', 'n', 'b', 'K', 'beta')
        self.assertCountEqual(paraList, correct_paraList)

        # model evaluation
        theta = {
            'alpha': 1.5,
            'n': 0.8,
            'b': 0.05,
            'K': 2E-3,
            'beta': 0.1
        }
        inducer = np.array([1E-4, 1E-3, 1E-2])
        reporter = model(inducer, theta)
        correct_reporter=[1.7098, 5.7898, 11.8641]
        for first, second in zip(reporter, correct_reporter):
            self.assertAlmostEqual(first, second, places=4)
        
    def test_genSSE(self):

        # first, generate the model
        specs = [ModelSpec.Activation, ModelSpec.Hill, ModelSpec.Basal_expression]
        expression, model, _, _ = modelBase.genModel(ModelType.Inducible, specs, plain_print=False)

        # Generate sse
        reporter = [1.0, 5.0, 11.0]
        inducer = [[1E-4], [1E-3], [1E-2]]
        paraList = ('alpha', 'n', 'b', 'beta', 'K')
        # prototype: genSSE(inducer, reporter, reporterStd, model, thetaList)
        sse = inference.genSSE(inducer, reporter, None, model, paraList)

        theta = (1.5, 0.8, 0.05, 0.1, 2E-3)
        sseVal = sse(theta)

        mu = [1.7098, 5.7898, 11.8641]
        correct_sse = sum([(x - y) ** 2 for x, y in zip(mu, reporter)])

        self.assertAlmostEqual(sseVal, correct_sse, places=3)

    def test_genJac(self):
        # Preperations
        specs = [ModelSpec.Activation, ModelSpec.Hill, ModelSpec.Basal_expression]
        exp, _, _, _ = modelBase.genModel(ModelType.Inducible, specs, plain_print=False)
        reporter = [1.0, 5.0, 11.0]
        std = [0.1, 0.1, 0.1]
        inducer = [[1E-4], [1E-3], [1E-2]]
        thetaList = ('alpha', 'n', 'b', 'beta', 'K')

        # generate jac
        jacobian = inference.genJac(inducer, reporter, std, exp, thetaList)
        initialTheta = (1.5, 0.8, 0.05, 0.1, 2E-3)
        #print(jacobian(initialTheta))

    def test_estimatePara(self):

        # Generate model and sse
        specs = [ModelSpec.Activation, ModelSpec.Hill, ModelSpec.Basal_expression]
        exp, model, constraints, _ = modelBase.genModel(ModelType.Inducible, specs, plain_print=False)
        print('Len of constraints', len(constraints))
        reporter = [1.7098, 5.7898, 11.8641]
        std = [0.1, 0.1, 0.1]
        inducer = [1E-4, 1E-3, 1E-2]
        thetaList = ('alpha', 'n', 'b', 'beta', 'K')
        print(thetaList)
        sse = inference.genSSE(inducer, reporter, None, model, thetaList)

        # Generate jacobian
        jacobian = inference.genJac(inducer, reporter, std, exp, thetaList)

        # Optimization
        initialTheta = (1.0, 1.0, 0.01, 0.8, 1E-3)
        # Try different methods
        for method in (ModelSolver.N_M, ModelSolver.BFGS, ModelSolver.SLSQP):
            print('\n' + method.name)
            startTime = time.time()
            theta, residue = inference.estimatePara(sse, initialTheta, jacobian, constraints, method = method)
            duration = time.time() - startTime
            print("Theta:", theta)
            print("Residue:", residue)
            print("Time elapse: %.3f"%duration)
            # The real theta:
            correct_theta = (1.5, 0.8, 0.05, 0.1, 2E-3)
            #for first, second in zip(theta, correct_theta):
            #    self.assertAlmostEqual(first, second, places=2)


if __name__ == '__main__':
    unittest.main()
