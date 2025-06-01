# This is a method to get the transpose of a matrix
# It has two methods, one for the normal matrices that work with numbers called transpose
# and another for the matrices that have symbols (i.e., x, y, z) called transpose_with_symbols
# The second one takes the original matrix and a dictionary
# The dictionary should contain a string that represents the symbol and the corresponding value for this symbol
# The corresponding value could be the same symbol (or any other symbol), but in sympy
# EX :
#       [[x, 2, 3], [4, y, 6]]
#       dic is {"x": x, "y": y}
# The return value will be: [[x, 4], [2, y], [3, 6]]





from MatrixBase import MatrixBase

class Transpose(MatrixBase):
    def __init__(self, matrix, substitutions=None):
        super().__init__(matrix, substitutions)

    def transpose(self):
        transposed_A = []
        rows, cols = self.symbolic_matrix.shape
        for i in range(cols):
            transposed_A.append([])
            for j in range(rows):
                transposed_A[i].append(self.symbolic_matrix[j, i])
        return transposed_A
