#!/usr/bin/python3



def cluster_atoms(coeff, vec_a, vec_b, n_layers, writefile):

    # Creation of vector matrix

    v_mat_all = []

    writefile.write('Calculating the vectors to all atoms of each neighbour set...\n\n')

    writefile.write('The following vectors were calculated:\n')

    for num_N in range(0, n_layers + 1):

        v_mat_test = []

        for N1 in range(0, len(coeff[num_N])):

            v_mat_test.append([coeff[num_N][N1][0] * vec_a[0] + coeff[num_N][N1][1] * vec_b[0],
                               coeff[num_N][N1][0] * vec_a[1] + coeff[num_N][N1][1] * vec_b[1]])



            #print(coeff[num_N][N1][0], coeff[num_N][N1][1], vec_a[0], vec_b[0], v_mat_test)

            writefile.write('Vector to atom ')
            writefile.write(str(N1 + 1))
            writefile.write(' of neighbour set ')
            writefile.write(str(num_N))
            writefile.write(' = ')
            writefile.write(str(v_mat_test[N1]))
            writefile.write('\n')

        writefile.write('\n')


        v_mat_all.append(v_mat_test)

    return v_mat_all


def cluster_set2_atoms(set1, shift, writefile):

    set2 = []

    for num1 in range(0, len(set1)):

        set2_test = []

        for num2 in range(0, len(set1[num1])):

            set2_test.append([set1[num1][num2][0] - s[0], set1[num1][num2][1] - shift[1]])

            writefile.write('Vector to atom ')
            writefile.write(str(num2 + 1))
            writefile.write(' of neighbour set ')
            writefile.write(str(num1))
            writefile.write('B = ')
            writefile.write(str(set2_test[num2]))
            writefile.write('\n')

        set2.append(set2_test)

        writefile.write('\n')

    return set2