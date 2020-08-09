#!/usr/bin/env python
# coding: utf-8

import networkx as nx
from pyomo.environ import *
from pyomo.opt import SolverFactory


# Input data as AMPL format
# Input data as AMPL format
model = ConcreteModel()

association = nx.read_edgelist('./association.txt')
colocalization = nx.read_edgelist('./colocalization.txt')
direct_interaction = nx.read_edgelist('./direct_interaction.txt')
Physical_association = nx.read_edgelist('./physical_association.txt')

network_layer = ['association','colocalization','direct_interaction','Physical_association']

# construct betweeness
association_betweeness = nx.edge_betweenness(association,normalized=False)
colocalization_betweeness = nx.edge_betweenness(colocalization,normalized=False)
direct_interaction_betweeness = nx.edge_betweenness(direct_interaction,normalized=False)
Physical_association_betweeness = nx.edge_betweenness(Physical_association,normalized=False)

association_bet = dict(association_betweeness)
colocalization_bet = dict(colocalization_betweeness)
direct_interaction_bet = dict(direct_interaction_betweeness)
Physical_association_bet = dict(Physical_association_betweeness)

association_node = list(association.nodes)
colocalization_node = list(colocalization.nodes)
direct_interaction_node = list(direct_interaction.nodes)
Physical_association_node = list(Physical_association.nodes)

# create node
node = list(association.nodes)
node2 = list(colocalization.nodes)
node3 = list(direct_interaction.nodes)
node4 = list(Physical_association.nodes)
node.extend(node2)
node.extend(node3)
node.extend(node4)
node = list(set(node))
#----------------------------------------------------------
# construct li and d

# degree of nodes in each layer
d ={}
li={}
value_association = 0
for i in association.nodes:
    d[(i, network_layer[0])] = association.degree[i]
    value_association += association.degree[i]
li[network_layer[0]] = value_association

value_colocalization = 0
for i in colocalization.nodes:
    d[(i, network_layer[1])] = colocalization.degree[i]
    value_colocalization += colocalization.degree[i]
li[network_layer[1]] = value_colocalization

value_direct_interaction = 0
for i in direct_interaction.nodes:
    d[(i, network_layer[2])] = direct_interaction.degree[i]
    value_direct_interaction += direct_interaction.degree[i]
li[network_layer[2]] = value_direct_interaction

value_Physical_association = 0
for i in Physical_association.nodes:
    d[(i, network_layer[3])] = Physical_association.degree[i]
    value_Physical_association += Physical_association.degree[i]
li[network_layer[3]] = value_Physical_association

#------------------------------------------------------------
# Construct alpha
alpha = {}

for i in node:
    if (i,i) not in association.edges:
        alpha[(i,network_layer[0])] = 0
    else:
        alpha[(i,network_layer[0])] = association_betweeness[i,i]

# network_layer_2
for i in node:
    if (i,i) not in colocalization.edges:
        alpha[(i,network_layer[1])] = 0
    else:
        alpha[(i,network_layer[1])] = colocalization_betweeness[i,i]

# network_layer_3
for i in node:
    if (i,i) not in direct_interaction.edges:
        alpha[(i,network_layer[2])] = 0
    else:
        alpha[(i,network_layer[2])] = direct_interaction_betweeness[i,i]

#network_layer_4
for i in node:
    if (i,i) not in Physical_association.edges:
        alpha[(i,network_layer[3])] = 0
    else:
        alpha[(i,network_layer[3])] = Physical_association_betweeness[i,i]

# ------------------------------------------
# construct belta
# belta
belta = {}

# Layer 1
for i in node:
    for j in node:
        if (i, j) not in association_bet.keys():
            weight = 0
            belta[(i, j, network_layer[0])] = weight
        else:
            weight_1 = association_bet[(i, j)]
            belta[(i, j, network_layer[0])] = weight_1

# Layer_2
for i in node:
    for j in node:
        if (i, j) not in colocalization_bet.keys():
            weight = 0
            belta[(i, j, network_layer[1])] = weight
        else:
            weight_1 = colocalization_bet[(i, j)]
            belta[(i, j, network_layer[1])] = weight_1

# Layer 3
for i in node:
    for j in node:
        if (i, j) not in direct_interaction_bet.keys():
            weight = 0
            belta[(i, j, network_layer[2])] = weight
        else:
            weight_1 = direct_interaction_bet[(i, j)]
            belta[(i, j, network_layer[2])] = weight_1

# Layer 4
for i in node:
    for j in node:
        if (i, j) not in Physical_association_bet.keys():
            weight = 0
            belta[(i, j, network_layer[3])] = weight
        else:
            weight_1 = Physical_association_bet[(i, j)]
            belta[(i, j, network_layer[3])] = weight_1

print("Begin construct model")
# define set ----- Node Name
model.n = Set(initialize=node)
model.e = Set(initialize=node)
# define Layer Number
model.Layer = Set(initialize=network_layer)

# the number of models, for example: number = 50
model.m = RangeSet(1, 30)
# define Li
model.l = Param(model.Layer, within=NonNegativeReals, initialize=li)
# define Dni
model.d = Param(model.n, model.Layer, within=NonNegativeReals, initialize=d)
# define weight of loop node
model.alpha = Param(model.n, model.Layer, within=NonNegativeReals, initialize=alpha)
# define belta_1
model.belta = Param(model.n, model.e, model.Layer, within=NonNegativeReals, initialize=belta)

# define variables
# equal to 1 if node n is in module m
model.Y = Var(model.n, model.m, within=Binary)
# sum of the weighted degrees (strength) of nodes in module m in network slice i
model.D = Var(model.m, model.Layer, within=NonNegativeReals)
# sum of the weights of the edges that are in module m in network slice i
model.L = Var(model.m, model.Layer, within=NonNegativeReals)
#
model.Dd = Var(model.n, model.Layer, within=NonNegativeReals)


# Objective and Constraint
def objrule(model):
    Q = sum([sum([model.L[m, i] / model.l[i] - pow((model.D[m, i] / model.l[i] * 2), 2) for m in model.m]) for i in
             model.Layer])
    return Q / len(model.Layer)


model.value = Objective(rule=objrule, sense=maximize)

model.value.pprint()

# Constraint 1
model.constraint_1 = Constraint(expr=sum([model.Y[n, m] for n in model.n for m in model.m]) == len(model.n))

model.constraint_1.pprint()


# Constraint 2
def constraint_2(model, n):
    return sum([model.Y[n, m] for m in model.m]) == 1


model.constraint_2 = Constraint(model.n, rule=constraint_2)
model.constraint_2.pprint()


# Constraint 3
def constraint_3(model, m, l):
    return model.D[m, l] == sum([model.Dd[n, l] * model.Y[n, m] for n in model.n])


model.constraint_3 = Constraint(model.m, model.Layer, rule=constraint_3)

model.constraint_3.pprint()


# Constraint 4
def constraint_4(model, m, l):
    alpha_value = sum([model.alpha[n, l] * model.Y[n, m] for n in model.n])
    belta_value = sum([model.belta[n, e, l] * model.Y[n, m] * model.Y[e, m] for n in model.n for e in model.n])
    return model.L[m, l] == alpha_value + belta_value


model.constraint_4 = Constraint(model.m, model.Layer, rule=constraint_4)

model.constraint_4.pprint()

SolverFactory('mindtpy').solve(model, mip_solver='glpk', nlp_solver='ipopt')
model.Y.pprint()
print('Finish')
model.display('./model.txt')