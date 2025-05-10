def Identity(n):
    id = []
    for i in range(n):
        id.append([])
        for j in range(n):
            if i == j:
                id[i].append(1)
            else:
                id[i].append(0)
    return id