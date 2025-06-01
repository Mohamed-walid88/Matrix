from sympy import Matrix, Symbol


class MatrixBase:
    def __init__(self, matrix, substitutions=None):
        self.original = Matrix(matrix)
        self.substitutions = substitutions or {}
        self.symbolic_matrix = self._apply_substitutions()

    def _apply_substitutions(self):
        """Apply substitutions consistently handling symbols"""
        resolved = {}
        for k, v in self.substitutions.items():
            if isinstance(k, str):
                # Convert string keys to symbols
                resolved[Symbol(k)] = v
            else:
                resolved[k] = v
        return self.original.subs(resolved)

    def is_square(self):
        return self.symbolic_matrix.is_square

    def validate_square(self):
        if not self.is_square():
            raise ValueError("Matrix must be square for this operation")

    def is_symbolic(self):
        """Check if matrix contains symbolic elements"""
        return any(element.is_Symbol for element in self.symbolic_matrix)

    def is_numeric(self):
        """Check if matrix contains numbers"""
        return any(element.is_number for element in self.symbolic_matrix)