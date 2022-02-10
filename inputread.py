#!/usr/bin/python3

import math
import random
import re

def read_input(writefile):
    # Default input information (Give input in mismatch.in)

    length_a = 1.0
    length_b = 1.0
    angle_phi = 60.0
    length_k = 1.0
    length_l = 1.0
    angle_theta = 90.0
    r_cut = 0.95
    do_plot = 1
    do_cluster_plot = 0
    angle_alpha = 0
    alpha_min = 0
    alpha_rand = 0
    final_run = 0
    n_set2 = 1
    n_set2_l = 0.6
    n_set2_angle = 30
    normalisation = 0
    n_layers = 100

    with open("mismatch.in", "r") as infile:
        cont_line = []
        for line in infile:
            cont_line.append(line)

    length_a_read = 0
    length_b_read = 0
    angle_phi_read = 0
    length_k_read = 0
    length_l_read = 0
    angle_theta_read = 0
    n_set2_read = 0
    n_set2_angle_read = 0
    n_set2_l_read = 0
    alpha_min_read = 0
    alpha_rand_read = 0
    angle_alpha_read = 0
    n_layers_read = 0
    normalisation_read = 0
    do_cluster_plot_read = 0
    do_plot_read = 0
    r_cut_read = 0

    for nline in range(0, len(cont_line)):
        cont_temp = str(cont_line[nline])

        cont_check_a = cont_temp.find('LENGTH_A')

        if cont_check_a != -1:
            length_test_a = re.findall('\d+\.\d+', cont_temp)
            length_a = float(length_test_a[0])
            length_a_read = 1

        cont_check_b = cont_temp.find('LENGTH_B')

        if cont_check_b != -1:
            length_test_b = re.findall('\d+\.\d+', cont_temp)
            length_b = float(length_test_b[0])
            length_b_read = 1

        cont_check_k = cont_temp.find('LENGTH_K')

        if cont_check_k != -1:
            length_test_k = re.findall('\d+\.\d+', cont_temp)
            length_k = float(length_test_k[0])
            length_k_read = 1

        cont_check_l = cont_temp.find('LENGTH_L')

        if cont_check_l != -1:
            length_test_l = re.findall('\d+\.\d+', cont_temp)
            length_l = float(length_test_l[0])
            length_l_read = 1

        cont_check_phi = cont_temp.find('ANGLE_PHI')

        if cont_check_phi != -1:
            angle_test_phi = re.findall('\d+\.\d+', cont_temp)
            angle_phi = float(angle_test_phi[0])
            angle_phi_read = 1

        cont_check_theta = cont_temp.find('ANGLE_THETA')

        if cont_check_theta != -1:
            angle_test_theta = re.findall('\d+', cont_temp)
            angle_theta = float(angle_test_theta[0])
            angle_theta_read = 1

        cont_check_n_set2 = cont_temp.find('N_SET2_IF')

        if cont_check_n_set2 != -1:
            n_set2_test = re.findall('\d+', cont_temp)
            n_set2 = int(n_set2_test[1])
            n_set2_read = 1

        cont_check_n_set2_l = cont_temp.find('N_SET2_L')

        if cont_check_n_set2_l != -1:
            n_set2_l_test = re.findall('\d+\.\d+', cont_temp)
            n_set2_l = float(n_set2_l_test[0])
            n_set2_l_read = 1

        cont_check_n_set2_angle = cont_temp.find('N_SET_ANGLE')

        if cont_check_n_set2_angle != -1:
            n_set2_angle_test = re.findall('\d+\.\d+', cont_temp)
            n_set2_angle = float(n_set2_angle_test[0])
            n_set2_angle_read = 1

        cont_check_plot = cont_temp.find('PLOT_FUNC')

        if cont_check_plot != -1:
            do_plot_test = re.findall('\d+', cont_temp)
            do_plot = int(do_plot_test[0])
            do_plot_read = 1

        cont_check_cplot = cont_temp.find('PLOT_CLUSTER')

        if cont_check_cplot != -1:
            do_cluster_plot_test = re.findall('\d+', cont_temp)
            do_cluster_plot = int(do_cluster_plot_test[0])
            do_cluster_plot_read = 1

        cont_check_r_cut = cont_temp.find('R_CUT')

        if cont_check_r_cut != -1:
            r_cut_test = re.findall('\d+\.\d+', cont_temp)
            r_cut = float(r_cut_test[0])
            r_cut_read = 1

        cont_check_alpha_min = cont_temp.find('ALPHA_MIN')

        if cont_check_alpha_min != -1:
            alpha_min_test = re.findall('\d+', cont_temp)
            alpha_min = int(alpha_min_test[0])
            alpha_min_read = 1

        cont_check_alpha_rand = cont_temp.find('ALPHA_RAND')

        if cont_check_alpha_rand != -1:
            alpha_rand_test = re.findall('\d+', cont_temp)
            alpha_rand = int(alpha_rand_test[0])
            alpha_rand_read = 1

        cont_check_angle_alpha = cont_temp.find('ANGLE_ALPHA')

        if cont_check_angle_alpha != -1:
            angle_alpha_test = re.findall('\d+\.\d+', cont_temp)
            angle_alpha = float(angle_alpha_test[0])
            angle_alpha_read = 1

        cont_check_n_neighbours = cont_temp.find('N_NEIGHBOURS')

        if cont_check_n_neighbours != -1:
            n_neighbours_test = re.findall('\d+', cont_temp)
            n_layers = int(n_neighbours_test[0])
            n_layers_read = 1

        cont_check_d_norm = cont_temp.find('D_NORM')

        if cont_check_d_norm != -1:
            d_norm_test = re.findall('\d+', cont_temp)
            normalisation = int(d_norm_test[0])
            normalisation_read = 1

    if length_a == length_k and length_b == length_l or length_a == length_l and length_b == length_k:
        if angle_phi == angle_theta:
            print('Both materials are the same. There is no mismatch. Exiting..')
            exit()

    writefile.write('##CLUSTER LATTICE DEFINITION##\n')

    if length_a_read == 1:
        writefile.write('l_a = ')
        writefile.write(str(length_a))
        writefile.write('\n')

    else:
        writefile.write('l_a = ')
        writefile.write(str(length_a))
        writefile.write('(default value)\n')

    if length_b_read == 1:
        writefile.write('l_b = ')
        writefile.write(str(length_b))
        writefile.write('\n')

    else:
        writefile.write('l_b = ')
        writefile.write(str(length_b))
        writefile.write('(default value)\n')

    if angle_phi_read == 1:
        writefile.write('phi = ')
        writefile.write(str(angle_phi))
        writefile.write('\n')

    else:
        writefile.write('phi = ')
        writefile.write(str(angle_phi))
        writefile.write('(default value)\n')

    writefile.write('\n')
    writefile.write('##SURFACE LATTICE DEFINITION##\n')

    if length_k_read == 1:
        writefile.write('l_k = ')
        writefile.write(str(length_k))
        writefile.write('\n')

    else:
        writefile.write('l_k = ')
        writefile.write(str(length_k))
        writefile.write('(default value)\n')

    if length_l_read == 1:
        writefile.write('l_l = ')
        writefile.write(str(length_l))
        writefile.write('\n')

    else:
        writefile.write('l_l = ')
        writefile.write(str(length_l))
        writefile.write('(default value)\n')

    if angle_theta_read == 1:
        writefile.write('theta = ')
        writefile.write(str(angle_theta))
        writefile.write('\n')

    else:
        writefile.write('theta = ')
        writefile.write(str(angle_theta))
        writefile.write('(default value)\n')

    writefile.write('\n')
    writefile.write('##SITE SET2 DEFINITION##\n')

    if n_set2_read == 1 and n_set2 == 1:
        writefile.write('Second set of surface sites detected.\n')

        if n_set2_l_read == 1:
            writefile.write('l_s = ')
            writefile.write(str(n_set2_l))
            writefile.write('\n')

        else:
            writefile.write('l_s = ')
            writefile.write(str(n_set2_l))
            writefile.write('(default value)\n')

        if n_set2_angle_read == 1:
            writefile.write('kappa = ')
            writefile.write(str(n_set2_angle))
            writefile.write('\n')

        else:
            writefile.write('kappa = ')
            writefile.write(str(n_set2_angle))
            writefile.write('(default value)\n')

    elif n_set2_read == 1 and n_set2 == 0:
        writefile.write('Second set of surface sites deactivated in input.\n')

    else:
        writefile.write('No second set of surface sites detected. Only one set will be used (Default).\n')

    writefile.write('\n')
    writefile.write('##ROTATION SETTINGS##\n')

    if alpha_min_read == 1:

        if alpha_min == 0:
            writefile.write('No optimisation of rotation angle alpha.\n')

        elif alpha_min == 1:
            writefile.write('Using Nelder-Mead algorithm for optimisation of alpha.\n')

    elif alpha_min_read == 0 and alpha_min == 0:
        writefile.write('No optimisation of rotation angle alpha (Default).\n')

    elif alpha_min_read == 0 and alpha_min == 1:
        writefile.write('Using Nelder-Mead algorithm for optimisation of alpha (Default).\n')

    if alpha_rand_read == 1:

        if alpha_rand == 0:

            if angle_alpha_read == 1:
                writefile.write('alpha = ')
                writefile.write(str(angle_alpha))
                writefile.write('\n')

            else:
                writefile.write('alpha = ')
                writefile.write(str(angle_alpha))
                writefile.write('(default value)\n')

        elif alpha_rand == 1:
            writefile.write('Random initialisation of alpha.\n')
            angle_alpha = random.randint(1,89)
            writefile.write('alpha = ')
            writefile.write(str(angle_alpha))
            writefile.write('\n')

        elif alpha_rand == 2:
            writefile.write('alpha set to (theta - phi) / 2.\n')
            angle_alpha = (angle_theta - angle_phi) / 2
            writefile.write('alpha = ')
            writefile.write(str(angle_alpha))
            writefile.write('\n')

    writefile.write('\n')
    writefile.write('##PLOT SETTINGS##\n')

    if do_plot_read == 1:
        if do_plot == 1:
            writefile.write('Plot of linear fit will be shown.\n')

        else:
            writefile.write('Plot of linear fit will NOT be shown.\n')

    else:
        writefile.write('Plot of linear fit will NOT be shown (default).\n')

    if r_cut_read == 1:
        writefile.write('R2_cut = ')
        writefile.write(str(r_cut))
        writefile.write('\n')

    else:
        writefile.write('R2_cut = ')
        writefile.write(str(r_cut))
        writefile.write(' (default value)\n')

    if do_cluster_plot_read == 1:
        if do_cluster_plot == 1:
            writefile.write('Plot of cluster orientation will be shown.\n')

        else:
            writefile.write('Plot of cluster orientation will NOT be shown.\n')

    else:
        writefile.write('Plot of cluster orientation will NOT be shown (default).\n')

    writefile.write('\n')
    writefile.write('##MISC SETTINGS##\n')

    if n_layers_read == 1:
        writefile.write('number of neighbour sets = ')
        writefile.write(str(n_layers))
        writefile.write('\n')

    else:
        writefile.write('number of neighbour sets = ')
        writefile.write(str(n_layers))
        writefile.write(' (default value)\n')

    if normalisation_read == 1:
        if normalisation == 0:
            writefile.write('d_N and m_N will not be normalised.\n')
        elif normalisation == 1:
            writefile.write('d_N and m_N will be normalised to (a+b)/2.\n')
        elif normalisation == 2:
            writefile.write('d_N and m_N will be normalised to (k+l)/2.\n')
        elif normalisation == 3:
            writefile.write('d_N and m_N will be normalised to (a+b+k+l)/4.\n')

    else:
        writefile.write('d_N and m_N will not be normalised (default).\n')

    writefile.write('\n')

    return length_a, length_b, math.radians(angle_phi),length_k, length_l, math.radians(angle_theta), \
           n_set2, n_set2_l, math.radians(n_set2_angle), alpha_min, alpha_rand, math.radians(angle_alpha), \
           r_cut, do_plot, do_cluster_plot, n_layers, normalisation