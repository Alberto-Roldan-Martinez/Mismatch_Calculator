#!/usr/bin/python3

import math

def vector_length_angle_neu(vec1, vec2_length, vec_ang):
    """returns vector 2 with length vec2_length and angle vec_angle to vec1"""

    vec2 = [(vec2_length * (vec1[0] * math.cos(vec_ang) + vec1[1] * math.sin(vec_ang)))
             / math.sqrt(pow(vec1[0], 2) + pow(vec1[1], 2)),
            (vec2_length * (vec1[1] * math.cos(vec_ang) + vec1[0] * math.sin(vec_ang)))
             / math.sqrt(pow(vec1[0], 2) + pow(vec1[1], 2))]

    return vec2

def vector_rotation_neu(vec, rot_angle):
    """rotates vector vec with angle rot_angle"""

    # Define rotation matrix rot_mat (angles in rad!)
    rot_mat = ([math.cos(rot_angle), -1 * math.sin(rot_angle)], [math.sin(rot_angle), math.cos(rot_angle)])

    # Calculate components of rotated vector

    vec_rot = [rot_mat[0][0] * vec[0] + rot_mat[0][1] * vec[1], rot_mat[1][0] * vec[0] + rot_mat[1][1] * vec[1]]

    return vec_rot

def surface_vectors(length_k, length_l, angle_theta, writefile):

    writefile.write('Calculating surface vectors...\n\n')

    vec_k = [length_k, 0]

    #print('vector k =', vec_k)

    writefile.write('vector k = ')
    writefile.write(str(vec_k))
    writefile.write('\n\n')

    vec_l = vector_length_angle_neu(vec_k, length_l, angle_theta)

    #print('vector l =', vec_l)

    writefile.write('vector l = ')
    writefile.write(str(vec_l))
    writefile.write('\n\n')

    return vec_k, vec_l

def cluster_vectors(length_a, length_b, angle_phi, writefile):

    writefile.write('Calculating cluster vectors...\n\n')

    # Definition of vector a of the cluster lattice (a,b); Orientation parallel to the x axis.
    vec_a = [length_a, 0]

    #print('vector a =', vec_a)
    writefile.write('vector a = ')
    writefile.write(str(vec_a))
    writefile.write('\n\n')

    # Calculation of vector b of the cluster lattice (a,b).
    vec_b = vector_length_angle_neu(vec_a, length_b, angle_phi)

    #print('vector b =', vec_b)

    writefile.write('vector b = ')
    writefile.write(str(vec_b))
    writefile.write('\n\n')

    return vec_a, vec_b


def shift_vector(length_s, angle_s, vec_k):

    vec_s = vector_length_angle_neu(vec_k, length_s, angle_s)

    return vec_s


def surface_sites(coeff, surf_vec1, surf_vec2):

    surf_vec_all = []

    for num1 in range(0, len(coeff)):

        surf_vec_test = []

        for num2 in range(0, len(coeff[num1])):
            surf_vec_temp = [coeff[num1][num2][0] * surf_vec1[0] + coeff[num1][num2][1] * surf_vec2[0],
                               coeff[num1][num2][0] * surf_vec1[1] + coeff[num1][num2][1] * surf_vec2[1]]
            surf_vec_test.append(surf_vec_temp)
        surf_vec_all.append(surf_vec_test)

    return surf_vec_all


def distance_btw_vec_neu(vec1, vec2):
    """Calculates the distance between two vectors"""

    dist = math.sqrt(pow((vec1[0] - vec2[0]), 2) + pow((vec1[1] - vec2[1]), 2))

    return dist
