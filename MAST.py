import Bio.Phylo
from io import StringIO


def MAST(tree1, tree2):
    # Get all leaf nodes from tree1 and tree2
    leaves1 = [node for node in tree1.find_clades() if not node.clades]
    leaves2 = [node for node in tree2.find_clades() if not node.clades]
    # Get number of leaf nodes in tree1 and tree2
    n = len(leaves1)
    m = len(leaves2)
    # Initialize a 2D matrix F with 0's
    F = [[0] * (m + 1) for i in range(n + 1)]

    # Fill the matrix F using dynamic programming
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # If the names of the leaf nodes match, increase the value in F
            if leaves1[i - 1].name == leaves2[j - 1].name:
                F[i][j] = F[i - 1][j - 1] + 1
            else:
                # If the names of the leaf nodes do not match, take the maximum value from the previous row or column
                F[i][j] = max(F[i - 1][j], F[i][j - 1])

    # Initialize variables i and j to be the last row and column of the matrix
    i, j = n, m
    mast = []
    # Traverse the matrix backwards to find the MAST
    while i > 0 and j > 0:
        if leaves1[i - 1].name == leaves2[j - 1].name:
            # If the names of the leaf nodes match, add the name to the MAST
            mast.append(leaves1[i - 1].name)
            i -= 1
            j -= 1
        elif F[i - 1][j] > F[i][j - 1]:
            # If the value in the previous row is greater, go to the previous row
            i -= 1
        else:
            # If the value in the previous column is greater, go to the previous column
            j -= 1

    # Return the MAST in reverse order
    return mast[::-1]


tree1 = Bio.Phylo.read(StringIO("((A,B),(C,D));"), "newick")
tree2 = Bio.Phylo.read(StringIO("((B,C),(A,D));"), "newick")

print(MAST(tree1, tree2))
