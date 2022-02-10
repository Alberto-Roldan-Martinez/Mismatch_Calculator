#!/usr/bin/python3

# Import function packages
from mpi4py import MPI
import math
import datetime
import scipy
from scipy import stats
from scipy import optimize

from vecfunc import surface_vectors
from vecfunc import cluster_vectors
from vecfunc import shift_vector
from vecfunc import vector_rotation_neu
#from vecfunc import vector_length_angle_neu
from vecfunc import surface_sites
from vecfunc import distance_btw_vec_neu

from atomfunc import cluster_atoms
from atomfunc import cluster_set2_atoms

from cfunc import rounded_coeff
#from cfunc import calc_coeff_neu
from cfunc import calc_coeff_all

from inputread import read_input

from plotfunc import structure_plot
from plotfunc import linear_fit
from plotfunc import fit_plot

from mmcalc import minimisation_function

#length_a = 1.0
#length_b = 1.0
#phi = 60
#length_k = 1.1
#length_l = 1.1
#theta = 90
#n_layers = 10000

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
n_proc = comm.Get_size()

if my_rank == 0:
    debugfile = open("debug.txt", 'w')
    debugfile.write("debug file created by process")
    debugfile.write(str(my_rank))
    debugfile.write("\n")
    debugfile.write("total number of processes is ")
    debugfile.write(str(n_proc))
    debugfile.write("\n")
    debugfile.close()

for proc in range(0, n_proc):
    debugfile = open("debug.txt", 'a')
    debugfile.write("this is process")
    debugfile.write(str(my_rank))
    debugfile.write("speaking")
    debugfile.write("\n")
    debugfile.close()

length_a = None
a = None
length_b = None
b = None
length_k = None
k = None
length_l = None
l = None
n_set2 = None
s = None
n_layers = None
alpha = None
alpha_min = None
coefficients = None
outfile = None
joblist = None
v_mat = None
v_mat2 = None
v_surf_mat = None
v_surf_mat2 = None
unique_dist_list = None
normalisation = None
r_cut = None

if my_rank == 0:
    outfile = open("mismatch_out.txt", 'w')

    print('Mismatch_MPI v1.00 by Julien Engel - 21/01/2020')
    outfile.write('Mismatch_MPI v1.00 by Julien Engel - 21/01/2020 \n \n')

    now = datetime.datetime.now()
    outfile.write('The current time is: ')
    outfile.write(str(now))
    outfile.write('\n \n')

    length_a, length_b, phi, length_k, length_l, theta, n_set2, n_set2_l, n_set2_angle, alpha_min, alpha_rand, \
    alpha, r_cut, do_plot, do_cluster_plot, n_layers, normalisation = read_input(outfile)

    outfile.write('calculating surface vectors \n')
    k, l = surface_vectors(length_k, length_l, theta, outfile)
    outfile.write('vector k:')
    outfile.write(str(k))
    outfile.write('\n')
    outfile.write('vector l:')
    outfile.write(str(l))
    outfile.write('\n\n')

    outfile.write('calculating cluster vectors \n')
    a, b = cluster_vectors(length_a, length_b, phi, outfile)
    outfile.write('vector a:')
    outfile.write(str(k))
    outfile.write('\n')
    outfile.write('vector b:')
    outfile.write(str(l))
    outfile.write('\n\n')

    if n_set2 == 1:

        s = shift_vector(n_set2_l, n_set2_angle, k)

    else:
        s = 0

# Distribute input values to all processes
length_a = comm.bcast(length_a)
length_b = comm.bcast(length_b)
length_k = comm.bcast(length_k)
length_l = comm.bcast(length_l)
a = comm.bcast(a)
b = comm.bcast(b)
k = comm.bcast(k)
l = comm.bcast(l)

alpha = comm.bcast(alpha)
alpha_min = comm.bcast(alpha_min)
n_set2 = comm.bcast(n_set2)
s = comm.bcast(s)
n_layers = comm.bcast(n_layers)
normalisation = comm.bcast(normalisation)
r_cut = comm.bcast(r_cut)

for proc in range(0, n_proc):
    debugfile = open("debug.txt", 'a')
    debugfile.write("variable k on process")
    debugfile.write(str(my_rank))
    debugfile.write("is")
    debugfile.write(str(k))
    debugfile.write("\n")
    debugfile.close()

n_coeff = 2 * round(math.sqrt(2 * n_layers))  # 10 + round(n_layers / 10)

if n_layers <= 2000:
    n_coeff = round(math.sqrt(6 * n_layers))

coeff_range = ((2 * n_coeff) + 1)

# Creation of the final lists of distances and coefficients on process 0
dist_list_total = []
coeff_list_total = []

if my_rank == 0:
    outfile.write('Calculating the distances from the origin of each neighbour set...\n\n')
    outfile.write('Evaluating coefficients up to ')
    outfile.write(str(n_coeff))
    outfile.write('. \n\n')

# Create distance and coefficient lists on each process
dist_list_proc = []
coeff_list_proc = []

if my_rank == 0:

    joblist = []

    for num in range(0, n_proc):
        joblist.append([])

    process = 0

    for n_a in range(-1 * n_coeff, n_coeff + 1):

        joblist[process].append(n_a)

        process = process + 1

        if process >= n_proc:
            process = 0

    debugfile = open("debug.txt", 'a')
    debugfile.write("joblist was created by process")
    debugfile.write(str(my_rank))
    debugfile.write("has the length")
    debugfile.write(str(len(joblist)))
    debugfile.write("and is as follows")
    debugfile.write(str(joblist))
    debugfile.write("\n")
    debugfile.close()

    for process in range(1, n_proc):

        debugfile = open("debug.txt", 'a')
        debugfile.write("process")
        debugfile.write(str(my_rank))
        debugfile.write("sends joblist")
        debugfile.write(str(joblist))
        debugfile.write("to process")
        debugfile.write(str(process))
        debugfile.write("\n")
        debugfile.close()

        comm.send(joblist[process], dest=process)

if my_rank != 0:
    queue = comm.recv(source=0)

else:
    queue = joblist[0]

for n_a in queue:

    for n_b in range(0, n_coeff + 1): # SYMMETRY!!

        dist = round(math.sqrt(pow((n_a * a[0] + n_b * b[0]), 2) + pow((n_a * a[1] + n_b * b[1]), 2)), 5)

        # print('lenght of ', n_a, 'a + ', n_b, 'b is ', dist)

        dist_list_proc.append(dist)
        coeff_list_proc.append([n_a, n_b])

        if n_b != -1 * n_b:
            dist_list_proc.append(dist)
            coeff_list_proc.append([-1 * n_a, -1 * n_b])


if my_rank != 0:

    list_package = [dist_list_proc, coeff_list_proc]

    debugfile = open("debug.txt", 'a')
    debugfile.write("process")
    debugfile.write(str(my_rank))
    debugfile.write("sends list package")
    debugfile.write(str(list_package))
    debugfile.write("to process")
    debugfile.write(str(0))
    debugfile.write("\n")
    debugfile.close()

    comm.send(list_package, dest=0)

elif my_rank == 0:

    for process in range(1, n_proc):

        list_package = comm.recv(source=process)

        debugfile = open("debug.txt", 'a')
        debugfile.write("process")
        debugfile.write(str(my_rank))
        debugfile.write("received list package")
        debugfile.write(str(list_package))
        debugfile.write("from process")
        debugfile.write(str(process))
        debugfile.write("\n")
        debugfile.close()

        dist_list_total.extend(list_package[0])
        coeff_list_total.extend(list_package[1])

    unique_dist_list = list(sorted(set(dist_list_total)))

    outfile.write('The following distances were found: \n')
    outfile.write('N   d_N\n')

    for num in range(0, len(unique_dist_list)):
        outfile.write(str(num))
        outfile.write('   ')
        outfile.write(str(unique_dist_list[num]))
        outfile.write('\n')
    outfile.write('\n')

    outfile.write('Calculating the coefficients of each neighbour set...\n\n')

    coefficients = []

    for N in range(0, n_layers + 1):

        coefficients_temp = []

        for num in range(len(dist_list_total)):

            if dist_list_total[num] == unique_dist_list[N]:
                coefficients_temp.append([coeff_list_total[num][0], coeff_list_total[num][1]])

        coefficients.append(coefficients_temp)

    outfile.write('The following coefficients were calculated: \n')

    for num1 in range(0, len(coefficients)):

        outfile.write('Neighbour set ')
        outfile.write(str(num1))
        outfile.write(' at distance ')
        outfile.write(str(unique_dist_list[num1]))
        outfile.write('\n\n')

        for num2 in range(0, len(coefficients[num1])):
            outfile.write('Atom ')
            outfile.write(str(num2 + 1))
            outfile.write(' = ')
            outfile.write(str(coefficients[num1][num2][0]))
            outfile.write('a + ')
            outfile.write(str(coefficients[num1][num2][1]))
            outfile.write('b')
            outfile.write('\n')

        outfile.write('\n')

    max_coeff = []

    for num in range(0, len(coefficients)):
        max_coeff.append(max(max(coefficients[num])))
        # print(coefficients[num], max(max(coefficients[num])))

    outfile.write('Maximum coefficient for ')
    outfile.write(str(n_layers))
    outfile.write(' neighbour sets is ')
    outfile.write(str(max(max_coeff[:(n_layers + 1)])))
    outfile.write('\n\n')

coefficients = comm.bcast(coefficients) # send coefficients to all processes
unique_dist_list = comm.bcast(unique_dist_list)


if alpha_min == 1:

    stop = 0

    if my_rank == 0:

        outfile.write('##Initiating optimisation procedure\n')
        outfile.write('\n')
        outfile.write('Nelder-Mead minimisation algorithm\n\n')
        result = scipy.optimize.minimize(minimisation_function, alpha, args=(comm, my_rank, n_proc, length_a, a,
                                                                             length_b, b, length_k, k, length_l, l,
                                                                             n_set2, s, coefficients, unique_dist_list,
                                                                             normalisation, n_layers, r_cut, outfile,
                                                                             stop),
                                         method='nelder-mead', options={'xtol': 1e-8, 'disp': True})
        outfile.write('##Optimisation procedure ended\n\n')
        outfile.write('Optimisation status ')
        outfile.write(str(result.message))
        outfile.write('\n\n')
        outfile.write('Number of iterations: ')
        outfile.write(str(result.nit))
        outfile.write('\n\n')
        outfile.write('##Final results of the mismatch calculation\n')
        outfile.write('\n')
        outfile.write('Mismatch between cluster interface (a,b,phi) and surface lattice (k,l,theta)\n')
        outfile.write('dmN/dr = ')
        outfile.write(str(result.fun))
        outfile.write('\n')
        outfile.write('Final orientation is rotated ')
        outfile.write(str(str(result.x).strip('[]')))
        outfile.write('°')
        outfile.write('\n\n')

        alpha = result.x[0]
        stop = 1

        for num in range(1, n_proc):
            comm.send(stop, dest=num)

    if my_rank != 0:

        temp_file = open("temp.txt", 'w')
        temp_file.close()

        while stop == 0:

            temp_file = open("temp.txt", 'a')

            #stop = comm.recv(source=0)

            if stop == 0:
                temp_file.write("process")
                temp_file.write(str(my_rank))
                temp_file.write("is starting minimisation function")
                temp_file.write(str(datetime.datetime.now()))
                temp_file.write("with stop = ")
                temp_file.write(str(stop))
                temp_file.write(".\n")

                alpha, stop = minimisation_function(alpha, comm, my_rank, n_proc, length_a, a, length_b, b, length_k, k,
                                              length_l, l, n_set2, s, coefficients, unique_dist_list, normalisation,
                                              n_layers, r_cut, temp_file, stop)


            temp_file.write("process")
            temp_file.write(str(my_rank))
            temp_file.write("is still here. stop = ")
            temp_file.write(str(stop))
            temp_file.write("at ")
            temp_file.write(str(datetime.datetime.now()))
            temp_file.write("with alpha = ")
            temp_file.write(str(alpha))
            temp_file.write("\n")
            temp_file.close()


# start here!

temp_file = open("temp.txt", 'a')
temp_file.write("process")
temp_file.write(str(my_rank))
temp_file.write("is going on at")
temp_file.write(str(datetime.datetime.now()))
temp_file.write(".\n")
temp_file.close()

if my_rank == 0:

    # rotate cell vectors

    arot = vector_rotation_neu(a, alpha)

    brot = vector_rotation_neu(b, alpha)

    # calculate cluster atom positions

    v_mat = cluster_atoms(coefficients, arot, brot, n_layers, outfile)

    if n_set2 == 1:  # If set 2 exists

        outfile.write('A second set of surface sites exist. Creating second cluster lattice shifted with -s...\n\n')

        v_mat2 = cluster_set2_atoms(v_mat, s, outfile)

    # Calculation of coefficient matrix ci from vector matrix v_mat
    outfile.write('Calculating coefficients c_i of cluster atom positions...')
    outfile.write('\n\n')

    ci = calc_coeff_all(v_mat, k, l, outfile)

    if n_set2 == 1:  # If set 2 exists
        outfile.write('Calculating coefficients c_i2 of shifted cluster atom positions...')
        outfile.write('\n\n')
        ci2 = calc_coeff_all(v_mat2, k, l, outfile)

    outfile.write('Calculating coefficients c_iR of occupied surface site positions...')
    outfile.write('\n\n')

    cir = rounded_coeff(ci)

    if n_set2 == 1:  # If set 2 exists
        outfile.write('Calculating coefficients c_iR2 of shifted occupied surface site positions...')
        outfile.write('\n\n')
        cir2 = rounded_coeff(ci2)

    outfile.write('Calculating surface site vectors...')
    outfile.write('\n\n')

    v_surf_mat = surface_sites(cir, k, l)

    if n_set2 == 1:  # If set 2 exists
        v_surf_mat2 = surface_sites(cir2, k, l)


joblist = None

if my_rank == 0:

    outfile.write('Evaluating distance between cluster atoms and nearest surface site...')
    outfile.write('\n\n')

    joblist = []

    process = 0

    elements = math.floor(len(v_mat) / n_proc) # minimum number of elements for each process

    rest = len(v_mat) % n_proc # number of remaining elements

    for num in range(0, n_proc):

        if rest > 0:

            elements_temp = elements + 1

            rest = rest - 1

        else:

            elements_temp = elements

        joblist.append(elements_temp)

    print("rest", rest)

    debugfile = open("debug.txt", 'a')
    debugfile.write("joblist was created by process")
    debugfile.write(str(my_rank))
    debugfile.write(".\n")
    debugfile.write("There are ")
    debugfile.write(str(len(v_mat)))
    debugfile.write("jobs for ")
    debugfile.write(str(n_proc))
    debugfile.write("processors.\n The joblist has ")
    debugfile.write(str(joblist))
    debugfile.write("elements.\n")
    debugfile.close()

joblist = comm.bcast(joblist)

d_mat_proc = []

startpoint = sum(joblist[:my_rank])

endpoint = sum(joblist[:(my_rank + 1)])

v_mat = comm.bcast(v_mat)
v_surf_mat = comm.bcast(v_surf_mat)

for num1 in range(startpoint, endpoint):

    d_mat_test = []

    for num2 in range(0, len(v_mat[num1])):
        d_mat_temp = distance_btw_vec_neu(v_mat[num1][num2], v_surf_mat[num1][num2])

        d_mat_test.append(d_mat_temp)

    d_mat_proc.append(d_mat_test)


if my_rank != 0:
    comm.send(d_mat_proc, dest=0)

elif my_rank == 0:

    #d_mat_all = []

    d_mat_all = d_mat_proc

    for process in range(1, n_proc):

        outfile.write('Receiving from process ')
        outfile.write(str(process))
        outfile.write(' and appending to d_mat.\n')

        d_mat_proc = comm.recv(source=process)

        for num1 in range(0, len(d_mat_proc)):

            outfile.write('Neighbour set ')
            outfile.write(str(num1))
            outfile.write('of transmission.\n\n')

            for num2 in range(0, len(d_mat_proc[num1])):
                outfile.write('Distance of Atom ')
                outfile.write(str(num2 + 1))
                outfile.write(' = ')
                outfile.write(str(d_mat_proc[num1][num2]))
                outfile.write('.\n')

            outfile.write('\n')
            outfile.write('Average distance in set = ')
            outfile.write(str(sum(d_mat_proc[num1]) / len(d_mat_proc[num1])))
            outfile.write('\n\n')

        d_mat_all.extend(d_mat_proc)

v_mat2 = comm.bcast(v_mat2)
v_surf_mat2 = comm.bcast(v_surf_mat2)

d_mat2_proc = []

if my_rank == 0:

    if n_set2 == 1:
        outfile.write('Evaluating distance between shifted cluster atoms and nearest surface site...')
        outfile.write('\n\n')

        outfile.write('WARNING: not correctly implemented yet!')

    if n_set2 == 1:

        for num1 in range(startpoint, endpoint):

            d_mat2_test = []

            for num2 in range(0, len(v_mat[num1])):
                d_mat2_temp = distance_btw_vec_neu(v_mat2[num1][num2], v_surf_mat2[num1][num2])

                d_mat2_test.append(d_mat2_temp)

            d_mat2_proc.append(d_mat2_test)

        if my_rank != 0:
            comm.send(d_mat2_proc, dest=0)

        elif my_rank == 0:

            d_mat2_all = []

            d_mat2_all.append(d_mat2_proc)

            for process in range(1, n_proc):

                outfile.write('Receiving from process ')
                outfile.write(str(process))
                outfile.write(' and appending to d_mat2.\n')

                d_mat2_proc = comm.recv(source=process)

                for num1 in range(0, len(d_mat2_proc)):

                    outfile.write('Neighbour set ')
                    outfile.write(str(num1))
                    outfile.write('of transmission.\n\n')

                    for num2 in range(0, len(d_mat2_proc[num1])):
                        outfile.write('Distance of Atom ')
                        outfile.write(str(num2 + 1))
                        outfile.write(' = ')
                        outfile.write(str(d_mat2_proc[num1][num2]))
                        outfile.write('.\n')

                    outfile.write('\n')
                    outfile.write('Average distance in set = ')
                    outfile.write(str(sum(d_mat2_proc[num1]) / len(d_mat2_proc[num1])))
                    outfile.write('\n\n')

                d_mat2_all.append(d_mat2_proc)



if my_rank == 0:

    structure_plot(v_mat, k, l, n_set2, s) #plotting grids

    if normalisation == 2:
        normalisation_factor = (length_k + length_l) / 2

    elif normalisation == 1:
        normalisation_factor = (length_a + length_b) / 2

    elif normalisation == 3:
        normalisation_factor = (length_a + length_b + length_k + length_l) / 4

    else:
        normalisation_factor = 1

    d_N_norm = []

    for num in range(0, len(unique_dist_list)):
        d_N_norm.append(unique_dist_list[num] / normalisation_factor)

    d_mat_norm = []

    print(d_mat_all[51])

    for num1 in range(0, len(d_mat_all)):

        d_mat_norm_temp = []

        for num2 in range(0, len(d_mat_all[num1])):
            print(num1, num2, d_mat_all[num1][num2], normalisation_factor)
            d_mat_norm_temp.append(d_mat_all[num1][num2] / normalisation_factor)

        d_mat_norm.append(d_mat_norm_temp)

    outfile.write('Calculating average mismatch with r...\n\n')

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

        outfile.write('Neighbour set ')
        outfile.write(str(num))
        outfile.write(' (')
        outfile.write(str(entries))
        outfile.write(' Atoms): Average mismatch with r = ')
        outfile.write(str(entry_sum / entries))
        outfile.write('\n')

    outfile.write('\n')

    av_mat_add = average_matrix

    fit_results_add = linear_fit(d_N_norm, coefficients, av_mat_add, r_cut, math.degrees(alpha), outfile)

    print("d_N_norm", d_N_norm)
    print("av_mat_add", av_mat_add)

    fit_plot(d_N_norm, av_mat_add, fit_results_add[0])

    print('Initial slope: ', fit_results_add[0])
    print('Final value: ', av_mat_add[-1])



if my_rank == 0:
    outfile.close()