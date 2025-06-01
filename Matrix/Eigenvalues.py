from MatrixBase import MatrixBase
from CharomaticPolynomial import CharomaticPolynomial
from sympy import Symbol, roots, simplify


class Eigenvalues(MatrixBase):
    def __init__(self, matrix, substitutions=None):
        super().__init__(matrix, substitutions)
        self.lambda_symbol = Symbol('λ')

    def eigenvalues(self):
        poly = CharomaticPolynomial(self.symbolic_matrix).berkowitz()
        root_dict = roots(poly, self.lambda_symbol, multiple=False)

        eigen_data = {}
        for root, multiplicity in root_dict.items():
            simplified_root = simplify(root)
            eigen_data[simplified_root] = multiplicity

        return eigen_data