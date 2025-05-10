def product(a, b):
    if len(a[0]) != len(b):
        return "Error: Matrix dimensions are not compatible for multiplication."

    result = [[0 for _ in range(len(b[0]))] for _ in range(len(a))]

    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(a[0])):
                result[i][j] += a[i][k] * b[k][j]

    return result