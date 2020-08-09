import networkx as nx
import json

association = nx.read_edgelist('/Users/shuhezheng/Desktop/zz/association.txt')
colocalization = nx.read_edgelist('/Users/shuhezheng/Desktop/zz/colocalization.txt')
direct_interaction = nx.read_edgelist('/Users/shuhezheng/Desktop/zz/direct_interaction.txt')
Physical_association = nx.read_edgelist('/Users/shuhezheng/Desktop/zz/physical_association.txt')

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

#print(alpha[('290','association')])
print(alpha)

print(len(alpha))

print(len(belta))