from sympy import Matrix, symbols



class Product:
    def __init__(self, matrix1, matrix2, substitutions=None):
        self.substitutions = substitutions or {}
        self.symbolic_matrix_a = self._apply_substitution(matrix1)
        self.symbolic_matrix_b = self._apply_substitution(matrix2)
        
    def _apply_substitution(self, matrix):
        resolved = {symbols(k) if isinstance(k, str) else k: v for k, v in self.substitutions.items()}
        return Matrix(matrix).subs(resolved)
        
    def evaluate(self):
        rows_A, cols_A = self.symbolic_matrix_a.shape
        rows_B, cols_B = self.symbolic_matrix_b.shape

        if cols_A != rows_B:
            raise ValueError("Matrix dimensions are not compatible for multiplication")
    
        result = [[0 for _ in range(cols_B)] for _ in range(rows_A)]
    
        for i in range(rows_A):
            for j in range(cols_B):
                for k in range(cols_A):
                    result[i][j] += self.symbolic_matrix_a[i, k] * self.symbolic_matrix_b[k, j]
    
        return Matrix(result)