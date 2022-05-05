import networkx as nx
import sys

def inorder_newick(G, root):
	'''
	Inorder recursive traversal and newick construction.
	Builds the children newick first, 
	then builds and adds parent newicks as siblings.
	Collapsed nodes are handled as siblings. 
	'''
	if len(list(G.neighbors(root))) != 0:
		#process children
		output = ''
		output += '('
		for child in G.neighbors(root):				
			output += inorder_newick(G, child)
			output += ","
		output = output[:-1]  # remove last comma
		output += ')'
		#process parent
		output += ','
		for seq in nodes[root]:					
			output += seq
			output += ','
		output = output[:-1]  # remove last comma
		#return '(' + output + '))'
		return '(' + output + ')'
	else:
		#base case for children
		return(','.join(nodes[root]))


def build_node_seqs_dict(nodes_file):
	nodes = {}
	with open(nodes_file, 'r') as f:
		for index, line in enumerate(f.readlines()):
			if index != 0:
				seq, node = line.strip().split(',')
				node = int(node)
				if node in nodes:
					nodes[node].append(seq)
				else:
					nodes[node] = [seq]
	return nodes


def build_edge_list(edges_file):
	edges = []
	with open(edges_file, 'r') as f:
		for index, line in enumerate(f.readlines()):
			if index != 0:
				src, trg, nmux = line.strip().split(',')
				edges.append((int(src), int(trg), int(nmux)))
	return edges





#What if instead of that we work directly on the newick string

#It's actually way more complicated
# def convert_nwk_binary(nwk):
# 	new = ''
# 	buffer = []
# 	pos = 0
# 	while pos < len(nwk):
# 		if nwk[pos] == '(':#Read child set
# 			buffer.append(nwk[pos]) #'('
# 			pos += 1
# 			while nwk[pos] != ')':
# 				seq_id = ''
# 				while nwk[pos] not in [',', ')']:
# 					seq_id += nwk[pos]
# 					pos += 1
# 				buffer.append(seq_id)
# 				if nwk[pos] == ',':
# 					buffer.append(nwk[pos]) # ','
# 					pos += 1
# 				#else, nwk[pos] = ')', so not incr pos => break the loop
# 			buffer.append(nwk[pos])
# 			#now buffer contains the sibling string between ( and )
# 			if buffer.count(',') < 2:
# 				new += ''.join(buffer)
# 			else:
# 				#[(a, b, c)] => [(a,b),c)]
# 				#[(a, b, c, d)] => [(a,b),(c,d)]
# 				#[[]]
# 				#read two chars at a time and group them
# 				mod_buf = []
# 				i = 0
# 				while i < len(buffer)
# 					c = buffer[i]
# 					if c in ['(', ',', ')']:
# 						mod_buf.append(c)
# 					else:
# 						while c not in ['(', ]

# 					i += 1


# 			pos += 1

# 		else:
# 			new += nwk[pos]
# 			pos += 1


#So let's work on the tree after all. First expand all collapsed nodes to children
#Then create hierarchy out of siblings

# def convert_binary(G):
# 	for n in G.nodes():
# 		if len(nodes[n]) > 1:
# 			#create a binary hierarchy out of internal nodes



# 		# seqs = []
# 		# if len(list(G.neighbors(n))) > 2:
# 		# 	for x in G.neighbors(n):
# 		# 		seqs.append(nodes[x])
# 			# pass
# '''
# 				if pos >= len(nwk): 
# 					raise('Unclosed parenthesis')

# '''

	


try:
	nodes_file = sys.argv[1]
	edges_file = sys.argv[2]
	seqs_file = sys.argv[3]
	nodes = build_node_seqs_dict(nodes_file)
	edges = build_edge_list(edges_file)

except:
	print('''
	nodes_file = sys.argv[1]
	edges_file = sys.argv[2]
	are missing. Running default test case.
		''')

	nodes = {
		0: ['a', 'b', 'c'],
		1: ['e'],
		2: ['d'],
		3: ['g','h'],
		4: ['x'],
	}

	edges = [
		(0, 1, 1),
		(0, 2, 1),
		(2, 4, 1),
		(4, 3, 1),
	]



# print('nodes', nodes)
# print('edges', edges)

G = nx.DiGraph()

for n in nodes.items():
	G.add_node(n[0])

for e in edges:
	G.add_edge(e[0], e[1])

#G_b = convert_binary(G)

# print('newick')
nwk = inorder_newick(G, 0) + ';'
print(nwk)

with open(sys.argv[4], 'w') as f:
	f.write(nwk)
