#!/usr/bin/python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

def structure_plot(cluster_vector_matrix, surface_vec1, surface_vec2, second_set, shift_vector):
    # Plot structure
    coord_range = 5
    surface_points_x = []
    surface_points_y = []

    for num1 in range(-10 * coord_range, coord_range * 10):
        for num2 in range(-10 * coord_range, coord_range * 10):
            x_temp = num1 * surface_vec1[0] + num2 * surface_vec2[0]
            y_temp = num1 * surface_vec1[1] + num2 * surface_vec2[1]

            if abs(x_temp) < coord_range and abs(y_temp) < coord_range:
                surface_points_x.append(x_temp)
                surface_points_y.append(y_temp)

    if second_set == 1:
        surface_points_x2 = []
        surface_points_y2 = []

        for num1 in range(-10 * coord_range, coord_range * 10):
            for num2 in range(-10 * coord_range, coord_range * 10):
                x_temp = num1 * surface_vec1[0] + num2 * surface_vec2[0] + shift_vector[0]
                y_temp = num1 * surface_vec1[1] + num2 * surface_vec2[1] + shift_vector[1]

                if abs(x_temp) < coord_range and abs(y_temp) < coord_range:
                    surface_points_x2.append(x_temp)
                    surface_points_y2.append(y_temp)

    cluster_points_x = []
    cluster_points_y = []

    cluster_points_x.append(0)
    cluster_points_y.append(0)

    for num1 in range(0, len(cluster_vector_matrix)):
        for num2 in range(0, len(cluster_vector_matrix[num1])):
            x_temp = cluster_vector_matrix[num1][num2][0]
            y_temp = cluster_vector_matrix[num1][num2][1]

            if abs(x_temp) < coord_range and abs(y_temp) < coord_range:
                cluster_points_x.append(x_temp)
                cluster_points_y.append(y_temp)


    # plt.plot(surface_points_x, surface_points_y, 'kh')
    plt.axes().set_aspect('equal')
    plt.plot(surface_points_x, surface_points_y, color='black', linestyle='none', marker='o', markersize=10)
    if second_set == 1:
        plt.plot(surface_points_x2, surface_points_y2, color='#59656d', linestyle='none', marker='o', markersize=10)
    plt.plot(cluster_points_x, cluster_points_y, color='#f5bf03', linestyle='none', marker='o', markersize=7)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig('cluster.png')
    plt.clf()
    #plt.show()
    
def linear_fit(neighbour_distances, neighbour_coefficients, average_matrix, r_cut, rot_angle, writefile):
    # Calculate linear fit

    plotfile = open("mismatch_plot_add.out", "w")

    writefile.write('Performing linear fit of average set mismatch...\n\n')
    writefile.write('Input data:\n')
    writefile.write('r    d_N\n')
    for num in range(0, len(average_matrix)):
        writefile.write(str(neighbour_distances[num]))
        writefile.write('    ')
        writefile.write(str(average_matrix[num]))
        writefile.write('\n')

    writefile.write('\n')

    plotfile.write('N')
    plotfile.write(' ')
    plotfile.write('r_N')
    plotfile.write(' ')
    plotfile.write('m_N')
    plotfile.write(' ')
    plotfile.write('n_N')
    plotfile.write('\n')

    for num in range(0, len(average_matrix)):
        plotfile.write(str(num))
        plotfile.write(' ')
        plotfile.write(str(neighbour_distances[num]))
        plotfile.write(' ')
        plotfile.write(str(average_matrix[num]))
        plotfile.write(' ')
        plotfile.write(str(len(neighbour_coefficients[num])))
        plotfile.write('\n')

    np_neighbour_distances = numpy.array(neighbour_distances)
    np_average_matrix = numpy.array(average_matrix)

    c = 2
    r_value = 1.0

    divergence = -1

    while r_value >= 0.8 * r_cut and c <= len(average_matrix):

        x_val = np_neighbour_distances[:c][:, numpy.newaxis]
        y_val = np_average_matrix[:c]

        slope, resid, rank, s = numpy.linalg.lstsq(x_val, y_val, rcond=None)

        r_value = 1 - resid / (np_average_matrix[:c].size * np_average_matrix[:c].var())

        if r_value > r_cut:
            divergence = slope.item(0)
            # intercept_final = slope.item(1)
            r_value_final = r_value
            c_final = c

        c = c + 1

    if divergence == -1:
        writefile.write('Results of linear fit \n specified R^2 could not be reached')
        plotfile.write('Results of linear fit \n specified R^2 could not be reached')

    else:
        writefile.write('Results of linear fit: \nd_ini = ')
        writefile.write(str(divergence))
        writefile.write('\n')
        writefile.write('R^2 = ')
        writefile.write(str(r_value_final[0]))
        writefile.write('\n')
        writefile.write('last neighbour included: ')
        writefile.write(str(c_final))
        writefile.write('\n')
        writefile.write('rotation angle alpha: ')
        writefile.write(str(rot_angle))
        writefile.write('\n')


        plotfile.write('Results of linear fit \n slope = ')
        plotfile.write(str(divergence))
        plotfile.write('\n')
        plotfile.write('R^2 = ')
        plotfile.write(str(r_value_final))
        plotfile.write('\n')
        plotfile.write('last neighbour included: ')
        plotfile.write(str(c_final))
        plotfile.write('\n')
        plotfile.write('rotation angle alpha: ')
        plotfile.write(str(rot_angle))
        plotfile.write('\n')
        plotfile.write('maximum value = ')
        plotfile.write(str(max(average_matrix)))
        plotfile.write('\n')
        plotfile.write('average value = ')
        plotfile.write(str(sum(average_matrix) / len(average_matrix)))
        plotfile.write('\n')
        plotfile.write('last value = ')
        plotfile.write(str(average_matrix[-1]))
        plotfile.write('\n')
        
        plotfile.close()

    return [divergence, r_value_final, c_final]

def fit_plot(neighbour_distances, average_matrix1, divergence1):

    np_neighbour_distances = numpy.array(neighbour_distances)
    np_average_matrix1 = numpy.array(average_matrix1)

    print("np_neighbour", np_neighbour_distances)
    print("np_average", np_average_matrix1)
    print("max_np_av", max(np_average_matrix1))

    m_max = max(np_average_matrix1)
    m_av = sum(np_average_matrix1) / len(np_average_matrix1)

    # Plot results

    plt.xlabel("r_N")
    plt.ylabel("m_N")
    plt.xlim(0, neighbour_distances[len(average_matrix1) - 1])
    plt.ylim(0, max(average_matrix1) + 0.2)
    plt.plot(np_neighbour_distances[:len(average_matrix1)], np_average_matrix1, '-o')
    plt.plot(np_neighbour_distances, divergence1 * np_neighbour_distances, 'r-')
    plt.plot(np_neighbour_distances, m_max + 0 * np_neighbour_distances, 'c-')
    plt.plot(np_neighbour_distances, m_av + 0 * np_neighbour_distances, 'y-')
    plt.savefig('plot.png')
    plt.clf()

