#!/usr/bin/python3

import math
import datetime

from vecfunc import vector_rotation_neu
from vecfunc import surface_sites
from vecfunc import distance_btw_vec_neu

from atomfunc import cluster_atoms
from atomfunc import cluster_set2_atoms

from cfunc import rounded_coeff
from cfunc import calc_coeff_all

from plotfunc import linear_fit


def minimisation_function(f_alpha, f_comm, f_my_rank, f_n_proc, f_length_a, f_a, f_length_b, f_b, f_length_k, f_k,
                          f_length_l, f_l, f_n_set2, f_s, f_coefficients, f_unique_dist_list, f_normalisation,
                          f_n_layers, f_r_cut, writefile, stop):

    if f_my_rank == 0:
        for num in range(1, f_n_proc):
            f_comm.send(stop, dest=num)

    else:
        stop = f_comm.recv(source=0)

        writefile.write("process")
        writefile.write(str(f_my_rank))
        writefile.write("received information. stop =")
        writefile.write(str(stop))
        writefile.write("at")
        writefile.write(str(datetime.datetime.now()))
        writefile.write(".\n")


    if stop == 0:

        v_mat = None
        v_mat2 = None
        v_surf_mat = None
        v_surf_mat2 = None

        alpha = f_comm.bcast(f_alpha)

        if f_my_rank == 0:

            # rotate cell vectors

            arot = vector_rotation_neu(f_a, alpha)

            brot = vector_rotation_neu(f_b, alpha)

            # calculate cluster atom positions

            v_mat = cluster_atoms(f_coefficients, arot, brot, f_n_layers, writefile)

            if f_n_set2 == 1:  # If set 2 exists

                writefile.write(
                    'A second set of surface sites exist. Creating second cluster lattice shifted with -s...\n\n')

                v_mat2 = cluster_set2_atoms(v_mat, f_s, writefile)

            # Calculation of coefficient matrix ci from vector matrix v_mat
            writefile.write('Calculating coefficients c_i of cluster atom positions...')
            writefile.write('\n\n')

            ci = calc_coeff_all(v_mat, f_k, f_l, writefile)

            if f_n_set2 == 1:  # If set 2 exists
                writefile.write('Calculating coefficients c_i2 of shifted cluster atom positions...')
                writefile.write('\n\n')
                ci2 = calc_coeff_all(v_mat2, f_k, f_l, writefile)

            writefile.write('Calculating coefficients c_iR of occupied surface site positions...')
            writefile.write('\n\n')

            cir = rounded_coeff(ci)

            if f_n_set2 == 1:  # If set 2 exists
                writefile.write('Calculating coefficients c_iR2 of shifted occupied surface site positions...')
                writefile.write('\n\n')
                cir2 = rounded_coeff(ci2)

            writefile.write('Calculating surface site vectors...')
            writefile.write('\n\n')

            v_surf_mat = surface_sites(cir, f_k, f_l)

            if f_n_set2 == 1:  # If set 2 exists
                v_surf_mat2 = surface_sites(cir2, f_k, f_l)

        joblist = None

        if f_my_rank == 0:

            writefile.write('Evaluating distance between cluster atoms and nearest surface site...')
            writefile.write('\n\n')

            joblist = []

            elements = math.floor(len(v_mat) / f_n_proc)  # minimum number of elements for each process

            rest = len(v_mat) % f_n_proc  # number of remaining elements

            for num in range(0, f_n_proc):

                if rest > 0:

                    elements_temp = elements + 1

                    rest = rest - 1

                else:

                    elements_temp = elements

                joblist.append(elements_temp)

            print("rest", rest)

            debugfile = open("debug.txt", 'a')
            debugfile.write("joblist was created by process")
            debugfile.write(str(f_my_rank))
            debugfile.write(".\n")
            debugfile.write("There are ")
            debugfile.write(str(len(v_mat)))
            debugfile.write("jobs for ")
            debugfile.write(str(f_n_proc))
            debugfile.write("processors.\n The joblist has ")
            debugfile.write(str(joblist))
            debugfile.write("elements.\n")
            debugfile.close()

        joblist = f_comm.bcast(joblist)

        d_mat_proc = []

        startpoint = sum(joblist[:f_my_rank])

        endpoint = sum(joblist[:(f_my_rank + 1)])

        v_mat = f_comm.bcast(v_mat)
        v_surf_mat = f_comm.bcast(v_surf_mat)

        for num1 in range(startpoint, endpoint):

            d_mat_test = []

            for num2 in range(0, len(v_mat[num1])):
                d_mat_temp = distance_btw_vec_neu(v_mat[num1][num2], v_surf_mat[num1][num2])

                d_mat_test.append(d_mat_temp)

            d_mat_proc.append(d_mat_test)

        if f_my_rank != 0:
            f_comm.send(d_mat_proc, dest=0)

        elif f_my_rank == 0:

            d_mat_all = d_mat_proc

            for process in range(1, f_n_proc):

                writefile.write('Receiving from process ')
                writefile.write(str(process))
                writefile.write(' and appending to d_mat.\n')

                d_mat_proc = f_comm.recv(source=process)

                for num1 in range(0, len(d_mat_proc)):

                    writefile.write('Neighbour set ')
                    writefile.write(str(num1))
                    writefile.write('of transmission.\n\n')

                    for num2 in range(0, len(d_mat_proc[num1])):
                        writefile.write('Distance of Atom ')
                        writefile.write(str(num2 + 1))
                        writefile.write(' = ')
                        writefile.write(str(d_mat_proc[num1][num2]))
                        writefile.write('.\n')

                    writefile.write('\n')
                    writefile.write('Average distance in set = ')
                    writefile.write(str(sum(d_mat_proc[num1]) / len(d_mat_proc[num1])))
                    writefile.write('\n\n')

                d_mat_all.extend(d_mat_proc)

        v_mat2 = f_comm.bcast(v_mat2)
        v_surf_mat2 = f_comm.bcast(v_surf_mat2)

        d_mat2_proc = []

        if f_my_rank == 0:

            if f_n_set2 == 1:
                writefile.write('Evaluating distance between shifted cluster atoms and nearest surface site...')
                writefile.write('\n\n')

                writefile.write('WARNING: not correctly implemented yet!')

            if f_n_set2 == 1:

                for num1 in range(startpoint, endpoint):

                    d_mat2_test = []

                    for num2 in range(0, len(v_mat[num1])):
                        d_mat2_temp = distance_btw_vec_neu(v_mat2[num1][num2], v_surf_mat2[num1][num2])

                        d_mat2_test.append(d_mat2_temp)

                    d_mat2_proc.append(d_mat2_test)

                if f_my_rank != 0:
                    f_comm.send(d_mat2_proc, dest=0)

                elif f_my_rank == 0:

                    d_mat2_all = []

                    d_mat2_all.append(d_mat2_proc)

                    for process in range(1, f_n_proc):

                        writefile.write('Receiving from process ')
                        writefile.write(str(process))
                        writefile.write(' and appending to d_mat2.\n')

                        d_mat2_proc = f_comm.recv(source=process)

                        for num1 in range(0, len(d_mat2_proc)):

                            writefile.write('Neighbour set ')
                            writefile.write(str(num1))
                            writefile.write('of transmission.\n\n')

                            for num2 in range(0, len(d_mat2_proc[num1])):
                                writefile.write('Distance of Atom ')
                                writefile.write(str(num2 + 1))
                                writefile.write(' = ')
                                writefile.write(str(d_mat2_proc[num1][num2]))
                                writefile.write('.\n')

                            writefile.write('\n')
                            writefile.write('Average distance in set = ')
                            writefile.write(str(sum(d_mat2_proc[num1]) / len(d_mat2_proc[num1])))
                            writefile.write('\n\n')

                        d_mat2_all.append(d_mat2_proc)

        if f_my_rank == 0:

            if f_normalisation == 2:
                normalisation_factor = (f_length_k + f_length_l) / 2

            elif f_normalisation == 1:
                normalisation_factor = (f_length_a + f_length_b) / 2

            elif f_normalisation == 3:
                normalisation_factor = (f_length_a + f_length_b + f_length_k + f_length_l) / 4

            else:
                normalisation_factor = 1

            d_N_norm = []

            for num in range(0, len(f_unique_dist_list)):
                d_N_norm.append(f_unique_dist_list[num] / normalisation_factor)

            d_mat_norm = []

            for num1 in range(0, len(d_mat_all)):

                d_mat_norm_temp = []

                for num2 in range(0, len(d_mat_all[num1])):
                    print(num1, num2, d_mat_all[num1][num2], normalisation_factor)
                    d_mat_norm_temp.append(d_mat_all[num1][num2] / normalisation_factor)

                d_mat_norm.append(d_mat_norm_temp)

            writefile.write('Calculating average mismatch with r...\n\n')

            average_matrix = []
            entries = 0
            entry_sum = 0

            debugfile = open("debug.txt", 'a')
            debugfile.write(str(d_mat_norm))
            debugfile.close()

            for num in range(0, len(d_mat_norm)):
                entries = entries + len(d_mat_norm[num])

                entry_sum = entry_sum + sum(d_mat_norm[num])

                average_matrix.append(entry_sum / entries)

                writefile.write('Neighbour set ')
                writefile.write(str(num))
                writefile.write(' (')
                writefile.write(str(entries))
                writefile.write(' Atoms): Average mismatch with r = ')
                writefile.write(str(entry_sum / entries))
                writefile.write('\n')

            writefile.write('\n')

            av_mat_add = average_matrix

            fit_results_add = linear_fit(d_N_norm, f_coefficients, av_mat_add, f_r_cut, math.degrees(alpha),
                                         writefile)

            print('Initial slope: ', fit_results_add[0])
            print('Final value: ', av_mat_add[-1])

            return fit_results_add[0]

        elif f_my_rank != 0:

            return alpha, stop

    else:
        return f_alpha, stop