import Identity
from sympy import Matrix, symbols
from Product import Product


class Power:
    def __init__(self, matrix, n, substitutions=None):
        self.matrix = matrix
        self.substitutions = substitutions or {}
        self.n = n
        self.symbolic_matrix = self._apply_substitution(self.matrix)

    def _apply_substitution(self, matrix):
        resolved = {symbols(k) if isinstance(k, str) else k: v for k, v in self.substitutions.items()}
        return Matrix(matrix).subs(resolved)

    def power(self):
        rows, cols = self.symbolic_matrix.shape
        if rows != cols:
            raise ValueError("Error: Matrix is not square")
        if self.n == 0:
            return Identity.Identity(rows)
        if self.n < 1 or not isinstance(self.n, int) :
            raise ValueError("Error: n must be a positive integer")

        result = Identity.Identity(rows)
        base = self.symbolic_matrix

        exp = self.n
        while exp > 0:
            if exp & 1 == 1:
                result = Product(result, base).evaluate()
            base = Product(base, base).evaluate()
            exp >>= 1

        return result