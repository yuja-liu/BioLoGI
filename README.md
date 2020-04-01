# BioLogic

![Biologic LOGO](src/logo.png)

Biologic is a Python package and a command line tool which infers the mechanistic model of a single biological node. Given characterization data, e.g., from a microplate reader, Biologic can automatically estimate model parameters and select the optimal model, based on mechanistic hypotheses.
"Biological nodes" are commonly regarded as *logic gates* in synthetic biology, but *biological nodes* is a more generalized term, since the nodes may have more than 2 inputs and may behave far from discretely.

## Quick Start

The Quick Start uses a [biological part](http://parts.igem.org/wiki/index.php?title=Part:BBa_K2791004) from the 2018 XJTU iGEM project. The part features a D-psicose inducible promoter.

```sh
./biologic -c config.json
```

It is easy to use the command line. It is recommended to store configurations in a file. Though, in the future, configurations may be passed directly to the command line as options.
The config file should have the fields as followed:

```json
{
    "units": {
        "inputUnits": "M",
        "outputUnits": "AU"
    },
    "tags": {
        "inputTags": ["Psicose"],
        "outputTag": "GFP",
        "experimentTags": ["Exp1", "Exp2"]
    },
    "dataPath": "data/BBa_K2791004.json",
    "modelSet": "Activation_System",
    "methods": {
        "paraEstimator": "Nelder_Mead",
        "infoCriteria": "AIC"
    }
}
```

The "units" and "tags" offers additional information to the input data. "dataPath" specifies the path to such input file. "modelSet" asks the user to choose a keyword which refers to the set of mechanistic hypotheses as model candidates.
The "method" field assigns the algorithms for parameter estimation or model selection.
Details of available keyword please refer full [Documentations](docs/modelBase.md).

The results would be like:

```
#1 model calculating...
A*alpha + K*b
-------------
 beta*(A + K)
Type = Inducible
Specs = Michaelis_Menten, Basal_expression, Activation
Parameters:
alpha = 1.7956828304895099e-06, beta = 6.715611346991595e-07, b = 7.634137876835703e-07, K = 0.003402318454950439
Residue = 2.2724
IC = -19.5914

... ...

#30 model calculating...
           n
          K *b
------------------------
     /                n\
     |     /       2 \ |
     | n   |  A*K_I  | |
beta*|K  + |---------| |
     |     | 2      2| |
     \     \I  + K_I / /
Type = Inducible
Specs = Hill, No_basal_expression, Inducer_Quadratic, Repression, Inducer_Repression, Inducer
Parameters:
beta = 0.02170608718198557, b = 0.22919855389527066, K = 3.3504046038232686e-09, n = 0.15253056333394976, A = 0.0014373552873615742, K_I = 0.0013820565707686108
Residue = 0.8132
IC = -26.3361
=============================================
Choice:
 2          2
A *alpha + K *b
---------------
      / 2    2\
 beta*\A  + K /
Type = Inducible
Specs = Quadratic, Basal_expression, Activation
Parameters:
alpha = 2.173222877354994, beta = 0.6233347166737514, b = 0.8823534145312986, K = 0.018853500592474218
IC = -38.206461
```

## Liscense

GPLv3
