print("> import tskit")
import tskit

print("> tree = tskit.load('./trees/tree-S0.3-R1_subset_250.trees')")
tree = tskit.load("./trees/tree-S0.3-R1.trees")

print("> tree.individual_locations:")
print(tree.individual_locations)

print("> len(tree.individual_locations)")
print(len(tree.individual_locations))

print("> tree.num_individuals")
print(tree.num_individuals)

print("> tree.individual_locations[1000]")
print(tree.individual_locations[1000])

print("> tree.indvidiuals()[1000]")
print(tree.individuals()[1000])

print("> tree.individuals()[1000].nodes")
print(tree.individuals()[1000].nodes)

print("> node1 = tree.individuals()[1000].nodes[0]")
node1 = tree.individuals()[1000].nodes[0]

print("> node2 = tree.individuals()[1000].nodes[1]")
node2 = tree.individuals()[1000].nodes[1]

print("> tree.nodes()[node1]")
print(tree.nodes()[node1])

print("> tree.nodes()[node2]")
print(tree.nodes()[node2])

print("> tree.nodes()[node2].individual")
print(tree.nodes()[node2].individual)

print("> tree.individuals()[1000].metadata['pedigree_id']")
print(tree.individuals()[1000].metadata['pedigree_id'])

print("> tree.num_samples")
print(tree.num_samples)

print("> tree.num_nodes")
print(tree.num_nodes)

print("> tree.num_individuals")
print(tree.num_individuals)

print("> tree.samples()[10]")
print(tree.samples()[10])

print("> tree.nodes()[tree.individuals()[10].nodes[0]")
print(tree.nodes()[tree.individuals()[10].nodes[0]])

print("> tree.nodes()[tree.individuals()[10].nodes[0].time")
print(tree.nodes()[tree.individuals()[10].nodes[0]].time)

print("> tree.samples()[36900]")
print(tree.samples()[36900])

print("> tree.nodes()[tree.individuals()[36900].nodes[0]")
print(tree.nodes()[tree.individuals()[36900].nodes[0]])

print("> tree.nodes()[tree.individuals()[36900].nodes[0].time")
print(tree.nodes()[tree.individuals()[36900].nodes[0]].time)

print("> tree.nodes()[tree.individuals()[36900].nodes[0].flags")
print(tree.nodes()[tree.individuals()[36900].nodes[0]].flags)

print("> tree.nodes()[tree.individuals()[40000].nodes[0]")
print(tree.nodes()[tree.individuals()[40000].nodes[0]])

print("> tree.nodes()[tree.individuals()[40000].nodes[0].time")
print(tree.nodes()[tree.individuals()[40000].nodes[0]].time)

print("> tree.nodes()[tree.individuals()[40000].nodes[0].flags")
print(tree.nodes()[tree.individuals()[40000].nodes[0]].flags)