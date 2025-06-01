from MatrixBase import MatrixBase
from Identity import Identity
from Product import Product


class Power(MatrixBase):
    def __init__(self, matrix, n, substitutions=None):
        super().__init__(matrix, substitutions)
        self.n = n

    def power(self):
        rows, cols = self.symbolic_matrix.shape
        if rows != cols:
            raise ValueError("Error: Matrix is not square")
        if self.n == 0:
            return Identity(rows)
        if self.n < 1 or not isinstance(self.n, int) :
            raise ValueError("Error: n must be a positive integer")

        result = Identity(rows)
        base = self.symbolic_matrix

        exp = self.n
        while exp > 0:
            if exp & 1 == 1:
                result = Product(result, base).evaluate()
            base = Product(base, base).evaluate()
            exp >>= 1

        return result
