from CharomaticPolynomial import CharomaticPolynomial
from Identity import Identity
from MatrixBase import MatrixBase


class Determinant(MatrixBase):
    def __init__(self, matrix, substitutions=None):
        super().__init__(matrix, substitutions)
        self._row_swaps = 0

    def _LU_decomposition(self):
        A = self.symbolic_matrix
        n = A.rows
        L = Identity(n)
        U = A.copy()
        P = Identity(n)
        sign = 1

        for i in range(n):
            max_row = max(range(i, n), key=lambda k: abs(U[k, i]))
            if max_row != i:
                U.row_swap(i, max_row)
                P.row_swap(i, max_row)
                sign *= -1
                if i > 0:
                    L.row_swap(i, max_row)
                    L.col_swap(i, max_row)

            if U[i, i] == 0:
                return 0

            for j in range(i + 1, n):
                factor = U[j, i] / U[i, i]
                L[j, i] = factor
                U[j, i:] -= factor * U[i, i:]

        det = sign
        for i in range(n):
            det *= U[i, i]

        return det.simplify()

    def _expand_along_line(self, A, fixed_index, is_row):
        n = A.shape[0]
        det = 0
        for j in range(n):
            i, j_ = (fixed_index, j) if is_row else (j, fixed_index)
            elem = A[i, j_]
            if elem == 0:
                continue
            minor = A.minor_submatrix(i, j_)
            sign = (-1) ** (i + j_)
            det += sign * elem * self._Laplace_expansion(minor)
        return det

    def _Laplace_expansion(self, A):
        n = A.shape[0]
        if n == 1:
            return A[0, 0]

        best_count = -1
        use_row = True
        index = 0

        for i in range(n):
            row_zeros = sum(1 for j in range(n) if A[i, j] == 0)
            col_zeros = sum(1 for j in range(n) if A[j, i] == 0)

            if row_zeros > best_count:
                best_count = row_zeros
                use_row = True
                index = i
            if col_zeros > best_count:
                best_count = col_zeros
                use_row = False
                index = i

        return self._expand_along_line(A, index, use_row)

    def _Bareiss_algorithm(self):
        A = self.symbolic_matrix.copy()
        n = A.rows
        det = 1

        for k in range(n - 1):
            if A[k, k].is_zero:
                swap_row = None
                for i in range(k + 1, n):
                    if not A[i, k].is_zero:
                        swap_row = i
                        break

                if swap_row is not None:
                    A.row_swap(k, swap_row)
                    det *= -1
                else:
                    return self._Berkowitz_det()

            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    numerator = A[i, j] * A[k, k] - A[i, k] * A[k, j]
                    denominator = 1 if k == 0 else A[k - 1, k - 1]
                    A[i, j] = numerator / denominator

        return det * A[n - 1, n - 1].simplify()

    def _Berkowitz_det(self):
        poly = CharomaticPolynomial(self.symbolic_matrix)
        p = poly.berkowitz()

        return (-1)**(self.symbolic_matrix.shape[0]) * p.subs('λ', 0)

    def getDet(self):
        n = self.symbolic_matrix.shape[0]
        if n <= 4 and self.is_symbolic() and not self.is_numeric():
            return self._Laplace_expansion(self.symbolic_matrix)
        elif self.is_numeric() and not self.is_symbolic():
            return self._LU_decomposition()
        else:
            return self._Bareiss_algorithm()
