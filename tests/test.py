import ioUtils
import model
from modelBase import ModelSpec, ModelType, ModelSet
from modelBase import defaultPara, genModel, genModelSet
from model import ModelSolver, ModelCriterion

print("test readConfig")
config = ioUtils.readConfig("config.json")
print(config)
print("test readInput")
data = ioUtils.readInput(config)
print(data)
print("test genModelSet")
modelSet = genModelSet(ModelSet["Simple_Inducible_Promoter"])
print(modelSet)
print("test genSSE")
exp, paraList, m = genModel(
        ModelType.Inducible, (
            ModelSpec.M_M, ModelSpec.No_basal_expression
        )
    )
sse = model.genSSE(data, paraList, m)
print("test paraEstimator")
default = defaultPara(paraList)
print(paraList, default)
para, res = model.paraEstimator(sse,
        default, ModelSolver[config["methods"]["paraEstimator"]])
print("exp:")
print(exp)
print("paraList:", paraList)
print("sse:", res)
print("parameters:", para)
