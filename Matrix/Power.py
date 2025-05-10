import Product, Identity

def Power(A,n):
    if len(A) != len(A[0]):
        return "Error: Matrix is not square"
    if n < 1 :
        return "Error: n must be a positive integer"

    result = Identity.Identity(len(A))
    base = A

    while n > 0:
        if n & 1 == 1:
            result = Product.product(result,base)
        base = Product.product(base,base)
        n >>= 1

    return result

print(Power([[1,2],[3,4]],2))