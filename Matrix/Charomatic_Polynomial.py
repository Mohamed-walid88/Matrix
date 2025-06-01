import sympy
from sympy import Matrix, NonSquareMatrixError, Symbol, symbols, zeros


class Charomatic_Polynomial:
    def __init__(self, matrix, substitutions=None):
        self.original = sympy.Matrix(matrix)
        self.substitutions = substitutions or {}
        self.symbolic_matrix = self._apply_substitutions()
        self.one = sympy.Integer(1)
        self._new = lambda coeffs: Matrix(coeffs)
        self._row_swaps = 0


    def _apply_substitutions(self):
        resolved = {k if isinstance(k, str) else k: v
                    for k, v in self.substitutions.items()}
        return self.original.subs(resolved)

    def berkowitz(self):
        if not self.symbolic_matrix:
            return sympy.Integer(1)

        if not self.symbolic_matrix.is_square:
            raise NonSquareMatrixError()

        A, N = self.symbolic_matrix, self.symbolic_matrix.shape[0]
        transforms = [None] * (N - 1)

        for n in range(N, 1, -1):
            T, k = zeros(n + 1, n), n - 1

            R = -A[k, :k]
            C = A[:k, k]

            A, a = A[:k, :k], -A[k, k]

            items = [C]
            for i in range(n - 2):
                items.append(A * items[i])

            for i, B in enumerate(items):
                items[i] = (R * B)[0, 0]

            items = [self.one, a] + items

            for i in range(n):
                T[i:, i] = items[: (n - i + 1)]

            transforms[k - 1] = T

        polys = [self._new([self.one, -A[0, 0]])]

        for i, T in enumerate(transforms):
            polys.append(T * polys[i])

        last_coeffs = tuple(polys[-1])

        lam = Symbol('λ')
        degree = len(last_coeffs) - 1

        poly_expr = sum(last_coeffs[i] * lam ** (degree - i) for i in range(degree + 1))

        return sympy.simplify(poly_expr)