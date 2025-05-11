# This is a class to get the transpose of a matrix
# It could be used for a numerical matrix, a symbolic one, or a mix of them
# Has two parameters, the first one is the original matrix, the other is a dictionary of the symbols you want to substitute 
# You could substitute the symbols with numbers or keep them as symbols
# EX :
#       [[x, 2, 3], [4, y, 6]]
#       dic is {"x": x, "y": y}
# The return value will be: [[x, 4], [2, y], [3, 6]]




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
