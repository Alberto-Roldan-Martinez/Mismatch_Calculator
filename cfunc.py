#!/usr/bin/python3

def rounded_coeff(coeff_matrix):

    rounded_coeff_all = []

    for num1 in range(0, len(coeff_matrix)):

        rounded_coeff_temp = []

        for num2 in range(0, len(coeff_matrix[num1])):

            rounded_coeff_temp.append([round(coeff_matrix[num1][num2][0]), round(coeff_matrix[num1][num2][1])])



        rounded_coeff_all.append(rounded_coeff_temp)

    return rounded_coeff_all


def calc_coeff_neu(cluster_vec, surf_vec1, surf_vec2):
    """The cluster vector is described as a linear combination of the surface vectors """

    coeff = [(cluster_vec[0] * surf_vec2[1] - cluster_vec[1] * surf_vec2[0])
              / (surf_vec1[0] * surf_vec2[1] - surf_vec1[1] * surf_vec2[0]),
             (cluster_vec[0] * surf_vec1[1] - cluster_vec[1] * surf_vec1[0])
              / (surf_vec1[1] * surf_vec2[0] - surf_vec1[0] * surf_vec2[1])]

    return coeff


def calc_coeff_all(v_matrix, surf_v1, surf_v2, writefile):
    ci_all = []

    for num1 in range(0, len(v_matrix)):  # range -1 than for v_mat_all

        ci_test = []

        writefile.write('Neighbour set ')
        writefile.write(str(num1))
        writefile.write('\n\n')

        for num2 in range(0, len(v_matrix[num1])):
            coeff_temp = calc_coeff_neu(v_matrix[num1][num2], surf_v1, surf_v2)

            #print(coeff_temp)

            ci_test.append(coeff_temp)

            writefile.write('Atom ')
            writefile.write(str(num2 + 1))
            writefile.write(' = ')
            writefile.write(str(coeff_temp[0]))
            writefile.write('k + ')
            writefile.write(str(coeff_temp[1]))
            writefile.write('l')
            writefile.write('\n')

        writefile.write('\n')

        ci_all.append(ci_test)

    return ci_all