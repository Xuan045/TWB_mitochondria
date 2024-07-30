from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'font.weight': 'bold'})

# Define the Clades from the image
CZ = Clade(name="", clades=[Clade(name="Z"), Clade(name="C")])
M = Clade(name="M", clades=[CZ, Clade(name="D"), Clade(name="E"), Clade(name="G"), Clade(name="Q")])

JT = Clade(name="", clades=[Clade(name="T"), Clade(name="J")])
U = Clade(name="U", clades=[Clade(name="K")])
R0 = Clade(name="", clades=[Clade(name="HV"), Clade(name="H"), Clade(name="V")])
R = Clade(name="R", clades=[Clade(name="P"), Clade(name="B"), Clade(name="F"), JT, U, R0])
N = Clade(name="N", clades=[Clade(name="Y"), Clade(name="O"), Clade(name="S"), Clade(name="A"), Clade(name="I"), Clade(name="X"), Clade(name="W"), R])
L3 = Clade(name="L3", clades=[M, N])
L = Clade(name="L*", clades=[L3])

# Main root clades
root = Clade(clades=[L])

# Optional: color mapping and other aesthetic features
base_colors = {
    "M": "#1A759F", "C": "#1E6091", "Z": "#1A759F",
    "E": "#3479A3", "G": "#168AAD", "D": "#6BA1B3",
    "N": "#95D5B2", "Y": "#52B788", "A": "#52B788",
    "F": "#D1B3C4", "B": "#B392AC", "R": "#B392AC",
    "H": "#797D7F", "J": "#717D7E", "K": "#7E846B",
    "P": "grey", "S": "grey", "O": "grey", "I": "grey", "X": "grey", "W": "grey", "T": "grey", "V": "grey", "HV": "grey", "U": "grey", "Q": "grey",
    "CZ": "grey", "JT": "grey", "R0": "grey", "L3": "grey", "L*": "grey"
}
# Drawing the tree with the root at the top
fig, ax = plt.subplots(figsize=(6, 6))
Phylo.draw(root, axes=ax, do_show=False, label_colors=base_colors)

ax.axis('off')
plt.savefig('phylogenetic_tree.png', dpi=300, bbox_inches='tight')

plt.show()