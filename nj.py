import Bio.Phylo as phylo
from io import StringIO

def menorvalor (matriz):
    menorvalor = float("inf")
    x,y= 0, 0
    for i in range(len(matriz)):
        for j in range(len(matriz[i])):
            if matriz[i][j] < menorvalor:
                menorvalor = matriz[i][j]
                x,y = i,j

    return x,y

def agrupamento_nomes (nomes, x,y):
    if y < x:
        x, y = y, x
    nomes[x] = "("+ nomes[x] + "," + nomes[y] + ")"
    del nomes[y]

def up_distancia (matriz, x,y):
    if y < x:
        x, y = y, x

    row =[]
    for i in range (0,x):
        row.append((matriz[x][i] + matriz[y][i])/2)
    matriz[x] = row
    for i in range(x+1, y):
            matriz[i][x] = (matriz[i][x]+matriz[y][i])/2
    for i in range(y+1, len(matriz)):
        matriz[i][x] = (matriz[i][x]+matriz[i][y])/2
        del matriz[i][y]
    del matriz[y]



def UPGMA(matriz, nomes):
    # Until all matriz have been joined...
    while len(nomes) > 1:
        # Locate lowest cell in the table
        x, y = menorvalor(matriz)

        # Join the table on the cell co-ordinates
        up_distancia(matriz, x, y)

        # Update the matriz accordingly
        agrupamento_nomes(nomes, x, y)
    # Return the final label
    return nomes[0]


def matriz_tsv(filename):
    with open(filename,"r") as f:
        f.readline()
        matriz=[]
        flag = 1
        for line in f:
            parts = line.rstrip().split('\t')

            for i in range(1):
                nums = parts[1+i: flag]
                matriz.append(nums)
            flag += 1

        for i in range(len(matriz)):
            for j in range(len(matriz[i])):
                matriz[i][j]=float(matriz[i][j])

    f.close()
    return matriz


def nomes_tsv (filename):

    with open(filename,"r") as f:
        f.readline()
        nomes=[]
        for line in f:
            parts = line.rstrip().split('\t')
            for i in range(1):
                nomes.append(parts[0])
    f.close()

    return nomes


tree= UPGMA(matriz_tsv("mytxtfile.txt"), nomes_tsv("mytxtfile.txt"))


phylo.draw(phylo.read(StringIO(tree), 'newick'))
