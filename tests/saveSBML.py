#!/bin/python

import simplesbml
model = simplesbml.sbmlModel()
model.addCompartment(1e-14, comp_id='comp')
model.addSpecies('E', 5e-21, comp='comp')
model.addSpecies('S', 1e-20, comp='comp')
model.addSpecies('P', 0.0, comp='comp')
model.addSpecies('ES', 0.0, comp='comp')
model.addReaction(['E', 'S'], ['ES'], 'comp*(kon*E*S-koff*ES)', local_params={'koff': 0.2, 'kon': 1000000.0}, rxn_id='veq')
model.addReaction(['ES'], ['E', 'P'], 'comp*kcat*ES', local_params={'kcat': 0.1}, rxn_id='vcat')

# output to file
fo = open('mytest.sbml', 'w')
with fo:
    fo.write(model.toSBML())
