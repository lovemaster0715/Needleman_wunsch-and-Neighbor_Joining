import sys

def needleman_wunsch(s1, s2, match, mismatch, gap):

    # tamanhos de linha e coluna
    tamanho_coluna = (len(s1))+1
    tamanho_linha = (len(s2))+1


    # criação da matriz com inicialização 0
    matrix = [[0] * tamanho_coluna for i in range(tamanho_linha)]


    # instruções das direções do caminho
    direc = {}
    for i in range(1, tamanho_coluna):
        direc[(0, i)] = (0, i - 1)
    for i in range(i, tamanho_linha):
        direc[(i, 0)] = (i - 1, 0)

    # guardar o maior valor
    info_max_value = 3 * [-1]


    # função de máximo valor
    def max_value(i, j):

        v1 = matrix[i - 1][j - 1] + (match if s1[j - 1] == s2[i - 1] else mismatch)
        v2 = matrix[i - 1][j] + gap
        v3 = matrix[i][j - 1] + gap

        max_value = max([0, v1, v2, v3])

        # checa o maior valor
        if max_value == v1:
            direc[(i, j)] = (i - 1, j - 1)
        elif max_value == v2:
            direc[(i, j)] = (i - 1, j)
        else:
            direc[(i, j)] = (i, j - 1)

        # checa o maior valor da casa da matriz
        if info_max_value[0] <= max_value:
            info_max_value[0], info_max_value[1], info_max_value[2] = max_value, i, j
        return max_value


    # preenche a matriz
    for i in range (1, tamanho_linha):
        for j in range (1, tamanho_coluna):
            matrix [i][j] =  max_value(i,j)

    # alinhamentos
    s1_align, s2_align = '', ''

    # índices com o último valor
    i, j = tamanho_linha - 1, tamanho_coluna - 1

    while i>0 and j>0:

        k, l = direc[(i,j)]

        if (i-1) == k and (j-1) == l:  # diagonal
            s1_align += s1[l]
            s2_align += s2[k]
        elif (i-1) == k and j == l:  # vertical
            s1_align += '-'
            s2_align += s2[k]
        elif i == k and (j-1) == l:  # horizontal
            s1_align += s1[l]
            s2_align += '-'

        i, j = k, l
        if not i and not j:
            break

    def get_score(s1, s2, match, mismatch, gap):
        score, identidade = 0,0
        for i in range(len(s1)):
            if s1[i] == s2[1]:
                identidade +=1
            if s1[i] != '-' and s2[i] != '-':
                if s1[i] == s2[i]:
                    score += match
                else:
                    score += mismatch
            else:
                score += gap
        #print("A identidade é ", "{:.2f}".format(identidade * 100 / (min(len(s1), len(s2)))), "%")
        return score

    #s1_align, s2_align = s1_align[::-1], s2_align[::-1]

    #print(s1_align)
    #print(s2_align,"\n")

    #print(get_score(s1_align, s2_align, match, mismatch, gap))
    return (get_score(s1_align, s2_align, match, mismatch, gap))

from Bio import SeqIO
f = open("matriz_distancia.txt", "a")
filename= ("seq1.fasta")
seq_object = SeqIO.read(filename,"fasta")
seq1=seq_object.seq
'''    if i == 1:
        print(seq1[9482:9732],)
    elif i == 2:
        print(seq1[177903:178153],
    elif i == 3:
    print(seq1[3519919:3520169],
    elif i == 4:
        print(seq1[1490249:1490499],
    elif i == 5:
        print(seq1[731123:731373],)'''

filename2 = open("seq6.fasta", 'r')
seq_object2 = SeqIO.read(filename2,"fasta")
seq2=seq_object2.seq
f.write(str(needleman_wunsch(seq1[9482:9732], seq2[731123:731373], 1, -1, -2)))

    #elif i == 6:
        #print(seq1[682942:683192],"\n")
    #elif i == 7:
#        print(seq1[454434:454284],"\n")






