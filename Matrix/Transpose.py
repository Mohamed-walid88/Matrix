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

def transpose_with_symbols(A, dec):
    resolved_dec = {dec[k] if isinstance(k, str) else k: v for k, v in dec.items()}
    def_A = Matrix(A).subs(resolved_dec)

    transposed_A = []
    for i in range(def_A.shape[1]):
        transposed_A.append([])
        for j in range(def_A.shape[0]):
            transposed_A[i].append(def_A[j, i])
    return transposed_A

def transpose(A):
    transposed_A = []
    for i in range(len(A[0])):
        transposed_A.append([])
        for j in range(len(A)):
            transposed_A[i].append(A[j][i])
    return transposed_A