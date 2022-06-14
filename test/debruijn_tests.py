## Single sequence
from biodive import debruijn

# Create k-mers from seq
# seqs = ["ATCGTTGCGCGACCG"]
seqs = ["ATCGTTG","ATCGTTG"]
kmers = [debruijn.build_k_mer(seq, 4) for seq in seqs]
print(kmers)

# Make graph
G = debruijn.make_debruijn_graph(kmers)
print(f"Nodes:\n{G[0]}\n\nEdges:\n{G[1]}\n\nStarts:\n{G[2]}")

# Make node-edge
m = debruijn.make_node_edge_map(G[1])
print(m)

# Find start of path
start = G[2][0] if (len(G[2]) > 0) else G[0][0]
print(start)

# Find trails
t = debruijn.eulerian_trail(m, start)
print(t)

# assemble
debruijn.assemble_trail(t)

## Multiple sequences
from biodive import debruijn

# Create k-mers from seq
seqs = ["ATCGTTGCGCGACC", "TCGTTGCGCGACCGT", "GTTGCGCGACCGTA", "TGCGCGACCGTAA"]
kmers = [debruijn.build_k_mer(seq, 6) for seq in seqs]
print(kmers)

# Make graph
G = debruijn.make_debruijn_graph(kmers)
print(f"Nodes:\n{G[0]}\n\nEdges:\n{G[1]}\n\nStarts:\n{G[2]}")

# Make node-edge
m = debruijn.make_node_edge_map(G[1])
print(m)

# Find start of path
start = G[2][0] if (len(G[2])>0) else G[0][0]
print(start)

# Find trails
t = debruijn.eulerian_trail(m, start)
print(t)

# assemble
assembly = debruijn.assemble_trail(t)
assembly == "ATCGTTGCGCGACCGTAA"

## All in one
from biodive import debruijn
keys = ["ATCGTTGCGCGACCG", "TCGTTGCGCGACCGT", "GTTGCGCGACCGTA", "TGCGCGACCGTAA"]
assembly = debruijn.asmbl_keys(keys)
assembly == "ATCGTTGCGCGACCGTAA"

## Repetitive fails!
from biodive import debruijn
keys = ["TTTTTT", "TTTTTT", "TTTTTT", "TTTTTT"]
assembly = debruijn.asmbl_keys(keys)
assembly == "TTTTTTTTT"
