# This is a method to get the transpose of a matrix
# it has two methods one for the normal matrices that works with numbers called transpose
# and another for the matrices that has symbols (i.e.: x, y, z) called transpose_with_symbols
# the second one takes the original matrix and a dictionary
# the dictionary should contain a string that represent the symbol and the corresponding value for this symbol
# The corresponding value could be the same symbol (or any other symbol) but in sympy
# EX :
#       [[x, 2, 3], [4, y, 6]]
#       dic is {"x": x, "y": y}
# the return value will be: [[x, 4], [2, y], [3, 6]]





from sympy import Matrix, symbols

class Transpose:
    def __init__(self, matrix, substitutions=None):
        self.original = matrix
        self.substitutions = substitutions or {}
        self.symbolic_matrix = self._apply_substitution()

    def _apply_substitution(self):
        resolved = {self.substitutions[k] if isinstance(k, str) else k: v for k, v in self.substitutions.items()}
        return Matrix(self.original).subs(resolved)

    def transpose(self):
        transposed_A = []
        rows, cols = self.symbolic_matrix.shape
        for i in range(cols):
            transposed_A.append([])
            for j in range(rows):
                transposed_A[i].append(self.symbolic_matrix[j, i])
        return transposed_A